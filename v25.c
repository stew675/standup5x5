// A solution to the Parker 5x5 Unique Word Problem
//
// Author: Stew Forster (stew675@gmail.com)	Date: Aug 2022
//

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>
#include <stdatomic.h>
#include <immintrin.h>

#define USE_AVX2_SCAN

#include "utilities.h"

// ********************* SOLVER ALGORITHM ********************

static void
add_solution(register uint32_t *solution)
{
	register int i, pos = atomic_fetch_add(&num_sol, 1);
	register char *so = solutions + pos * 30;
	register const char *wd;

	for (i = 1; i < 6; i++) {
		wd = hash_lookup(solution[i], words);
		assert(wd != NULL);

		*so++ = *wd++; *so++ = *wd++; *so++ = *wd++; *so++ = *wd++;
		*so++ = *wd; *so++ = (i < 5) ? '\t' : '\n';
	}
} // add_solution

static inline uint32_t
vscan(register uint32_t mask, register uint32_t *set)
{
	__m256i vmask = _mm256_set1_epi32(mask);
	__m256i vkeys = _mm256_loadu_si256((__m256i *)set);
	__m256i vres = _mm256_cmpeq_epi32(_mm256_and_si256(vmask, vkeys), _mm256_setzero_si256());
	return (uint32_t)_mm256_movemask_epi8(vres);
} // vscan

// find_solutions() which is the busiest loop is kept
// as small and tight as possible for the most speed
void
find_solutions(register int depth, register struct frequency *f, register uint32_t mask,
		register int skipped, register uint32_t *solution, register uint32_t key)
{
	solution[depth] = key;
	if (depth == 5)
		return add_solution(solution);
	mask |= key;

	register struct frequency *e = frq + (min_search_depth + depth);
	for (; f < e; f++) {
		if (mask & f->m)
			continue;

		register struct tier *t = f->sets + !!(mask & f->tm1) + (2 * !!(mask & f->tm2));

		// Determine the values for set and end
		// The !! means we end up with only 0 or 1
		register int mf = !!(mask & f->tm3);
		register int ms = !!(mask & f->tm4);

		// A branchless calculation of end
		register uint32_t *end = t->s + (ms * t->toff3) + (!ms * t->l);

		// A branchless calculation of set
		ms &= !mf;
		register uint32_t *set = t->s + ((mf & !ms) * t->toff2) + (ms * t->toff1);

		for (; set < end; set += 8) {
			register uint32_t vresmask = vscan(mask, set);
			while (vresmask) {
				register int i = __builtin_ctz(vresmask);
				find_solutions(depth + 1, f + 1, mask, skipped, solution, set[i>>2]);
				vresmask ^= (0xFU << i);
			}
		}
		if (skipped)
			break;
		skipped = 1;
	}
} // find_solutions

// Thread driver
static void
solve_work()
{
	uint32_t solution[6] __attribute__((aligned(64)));
	register struct frequency *f = frq;
	register struct tier *t;
	register int32_t pos;

	// Solve starting with least frequent set
	t = f->sets;
	while ((pos = atomic_fetch_add(&f->pos, 1)) < t->l)
		find_solutions(1, f + 1, 0, 0, solution, t->s[pos]);

	// Solve after skipping least frequent set
	f++;
	t = f->sets;
	while ((pos = atomic_fetch_add(&f->pos, 1)) < t->l)
		find_solutions(1, f + 1, 0, 1, solution, t->s[pos]);

	atomic_fetch_add(&solvers_done, 1);
} // solve_work

void
solve()
{
	// Instruct any waiting worker-threads to start solving
	start_solvers();

	// The main thread also participates in finding solutions
	solve_work();

	// Wait for all solver threads to finish up
	while(solvers_done < nthreads)
		asm("nop");
} // solve


// ********************* MAIN SETUP AND OUTPUT ********************

int
main(int argc, char *argv[])
{
	struct timespec t1[1], t2[1], t3[1], t4[1], t5[1];
	char file[256];

	// Copy in the default file-name
	strcpy(file, "words_alpha.txt");

	nthreads = get_nthreads();

	if (argc > 1) {
		for (int i = 1; i < argc; i++) {
			if (!strncmp(argv[i], "-v", 2)) {
				write_metrics = 1;
				continue;
			}

			if (!strncmp(argv[i], "-f", 2)) {
				if ((i + 1) < argc) {
					strncpy(file, argv[i+1], 255);
					file[255] = '\0';
					i++;
					continue;
				}
			}

			if (!strncmp(argv[i], "-t", 2)) {
				if ((i + 1) < argc) {
					nthreads = atoi(argv[i+1]);
					i++;
					if (nthreads < 0)
						nthreads = 1;
					if (nthreads > MAX_THREADS)
						nthreads = MAX_THREADS;
					continue;
				}
			}

			printf("Usage: %s [-v] [-t num_threads] [-f filename]\n", argv[0]);
			exit(1);
		}
	}

	if (nthreads <= 0)
		nthreads = 1;
	if (nthreads > MAX_THREADS)
		nthreads = MAX_THREADS;

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t1);

	read_words(file);

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t2);

	setup_frequency_sets(8);

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t3);

	solve();

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t4);

	emit_solutions();

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t5);

	if (!write_metrics)
		exit(0);

	printf("\nFrequency Table:\n");
	for (int i = 0; i < 26; i++) {
		struct tier *t = frq[i].sets;
		char c = 'a' + __builtin_ctz(frq[i].m);
		printf("%c set_length=%4d  toff1=%4d, toff2=%4d, toff[3]=%4d\n",
			c, t->l, t->toff1, t->toff2, t->toff3);
	}
	printf("\n\n");

	printf("Num Unique Words  = %8d\n", nkeys);
	printf("Hash Collisions   = %8u\n", hash_collisions);
	printf("Number of threads = %8d\n", nthreads);

	printf("\nNUM SOLUTIONS = %d\n", num_sol);

	printf("\nTIMES TAKEN :\n");
	print_time_taken("Total", t1, t5);
	printf("\n");
	print_time_taken("File Load", t1, t2);
	print_time_taken("Frequency Set Build", t2, t3);
	print_time_taken("Main Algorithm", t3, t4);
	print_time_taken("Emit Results", t4, t5);

	exit(0);
} // main
