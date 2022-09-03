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

#ifdef __AVX512F__

#define USE_AVX2_SCAN

#define	MAX_SOLUTIONS	8192
#define	MAX_WORDS	8192
#define	MAX_THREADS	  24

static const char	*solution_filename = "solutions.txt";

// Worker thread state
static struct worker {
	char	*start;
	char	*end;
} workers[MAX_THREADS] __attribute__ ((aligned(64)));

// Set Pointers
struct tier {
	uint32_t	*s;	// Pointer to set
	uint32_t	l;	// Length of set
	uint32_t	toff1;	// Tiered Offset 1
	uint32_t	toff2;	// Tiered Offset 2
	uint32_t	toff3;	// Tiered Offset 3
};

// Character frequency recording
static struct frequency {
	struct tier	sets[2];	// Tier Sets
	uint32_t	m;		// Mask (1 << (c - 'a'))
	uint32_t	tm1;		// Tiered Mask 1
	uint32_t	tm2;		// Tiered Mask 2
	uint32_t	tm3;		// Tiered Mask 3

	uint32_t	pad[14];	// Pad to 64 bytes
	int32_t		f;		// Frequency
	atomic_int	pos;		// Position within a set
} frq[26] __attribute__ ((aligned(64)));

// Keep frequently modified atomic variables on their own CPU cache line
atomic_int 	num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done = 0;
atomic_int	solvers_done = 0;

// Put all general global variables on their own cache line
static int32_t	min_search_depth __attribute__ ((aligned(64))) = 0;
static int	write_metrics = 0;
static int	nthreads = 0;
static int	nkeys = 0;
static uint32_t hash_collisions = 0;

// We build the solutions directly as a character array to write out when done
static char	solutions[MAX_SOLUTIONS * 30] __attribute__ ((aligned(64)));

#include "utilities.h"

// ********************* SOLVER ALGORITHM ********************

static void
add_solution(uint32_t *solution)
{
	int i, pos = atomic_fetch_add(&num_sol, 1);
	char *so = solutions + pos * 30;
	const char *wd;

	for (i = 1; i < 6; i++) {
		wd = hash_lookup(solution[i], words);
		assert(wd != NULL);

		*so++ = *wd++; *so++ = *wd++; *so++ = *wd++; *so++ = *wd++;
		*so++ = *wd; *so++ = (i < 5) ? '\t' : '\n';
	}
} // add_solution

static inline uint16_t
vscan(uint32_t mask, uint32_t *set)
{
	__m512i vmask = _mm512_set1_epi32(mask);
	__m512i vkeys = _mm512_loadu_si512((__m256i *)set);
	return (uint16_t) _mm512_cmpeq_epi32_mask(_mm512_and_si512(vmask, vkeys), _mm512_setzero_si512());
} // vscan

// find_solutions() which is the busiest loop is kept
// as small and tight as possible for the most speed
void
find_solutions(int depth, struct frequency *f, uint32_t mask,
		int skipped, uint32_t *solution, uint32_t key)
{
	solution[depth] = key;
	if (depth == 5)
		return add_solution(solution);
	mask |= key;

	struct frequency *e = frq + (min_search_depth + depth);

	for (; f < e; f++) {
		if (mask & f->m)
			continue;

		struct tier *t = f->sets + !!(mask & f->tm1);

		// Determine the values for set and end
		// The !! means we end up with only 0 or 1
		int mf = !!(mask & f->tm2);
		int ms = !!(mask & f->tm3);

		// A branchless calculation of end
		uint32_t *end = t->s + (ms * t->toff3) + (!ms * t->l);

		// A branchless calculation of set
		ms &= !mf;
		uint32_t *set = t->s + ((mf & !ms) * t->toff2) + (ms * t->toff1);

		for (; set < end; set += 16) {
			uint16_t vresmask = vscan(mask, set);
			while (vresmask) {
				int i = __builtin_ctz(vresmask);
				find_solutions(depth + 1, f + 1, mask, skipped, solution, set[i]);
				vresmask ^= ((uint16_t)1 << i);
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
	struct frequency *f = frq;
	struct tier *t = f->sets;
	int32_t pos;

	// Solve starting with least frequent set
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

	setup_frequency_sets(16);

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t3);

	solve();

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t4);

	emit_solutions();

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t5);

	if (!write_metrics)
		exit(0);

	printf("\nFrequency Table:\n");
	for (int i = 0; i < 26; i++) {
		char c = 'a' + __builtin_ctz(frq[i].m);
		struct tier *t = frq[i].sets;
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
#else
int
main()
{
	fprintf(stderr, "AVX-512 not supported on this system\n");
	exit(1);
} // main
#endif
