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

#define	MAX_SOLUTIONS	8192
#define	MAX_WORDS	8192
#define	MAX_THREADS	  24

static const char	*solution_filename = "solutions.txt";

// Worker thread state
static struct worker {
	char     *start;
	char     *end;
} workers[MAX_THREADS] __attribute__ ((aligned(64)));

// Character frequency recording
static struct frequency {
	uint32_t  *s;		// Pointer to set
	struct frequency *e;	// Min Frequency Depth
	uint32_t   m;		// Mask (1 << (c - 'a'))
	int32_t    f;		// Frequency
	int32_t    l;		// Length of set
	atomic_int pos;		// Position within a set
} frq[26] __attribute__ ((aligned(64)));

// Keep frequently modified atomic variables on their own CPU cache line
volatile atomic_int num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done = 0;
atomic_int	solvers_done = 0;

// Put general global variables on their own cache line
static int32_t	min_search_depth __attribute__ ((aligned(64))) = 0;
static int	write_metrics = 0;
static int	nthreads = 0;
static int	nkeys = 0;
static uint32_t hash_collisions = 0;

// We build the solutions directly as a character array to write out when done
static char	solutions[MAX_SOLUTIONS * 30] __attribute__ ((aligned(64)));

#include "utilities.h"

// The role of this function is to re-arrange the key set according to all
// words containing the least frequently used letter, and then scanning
// the remainder and so on until all keys have been assigned to sets
// It achieves this by swapping keys in the key set, and inserting values
// to ensure each set is properly aligned for vectorized scanning
void
setup_frequency_sets()
{
	register struct frequency *f = frq;
	register uint32_t *kp = keys;

	qsort(f, 26, sizeof(*f), by_frequency_lo);

	// Now set up our scan sets by lowest frequency to highest
	for (int i = 0; i < 26; i++, f++) {
		register uint32_t mask, *ks, key;

		if (i == 7)
			rescan_frequencies(i, kp);

		mask = f->m;
		f->s = kp;
		for (ks = kp; (key = *ks); ks++) {
			if (key & mask) {
				*ks = *kp;
				*kp++ = key;
			}
		}
		register uint32_t nfk = kp - f->s;
		f->l = nfk;

		// Update the min_search_depth if needed
		if (nfk > 0)
			min_search_depth = i - 3;

		// Ensure 32-byte alignments for the set
		// We "poison" the alignment values with all bits set
		while ((nfk++ % 8)) {
			*ks++ = *kp;
			*kp++ = (uint32_t)(~0);
		}

		// Ensure key set is 0 terminated for next loop
		*ks = 0;
	}
	for (int i = 0; i < 26; i++)
		frq[i].e = frq + min_search_depth;
} // setup_frequency_sets


// ********************* SOLVER ALGORITHM ********************

static void
add_solution(uint32_t *solution)
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
	__m256i vkeys = _mm256_load_si256((__m256i *)set);
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

	for (register struct frequency *e = f->e + depth; f < e; f++) {
		if (mask & f->m)
			continue;
		for (register uint32_t *set = f->s, *end = f->s + f->l; set < end; set += 8) {
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
void
solve_work()
{
	uint32_t solution[6] __attribute__((aligned(64)));
	register struct frequency *f = frq;
	register int32_t pos;

	// Solve starting with least frequent set
	while ((pos = atomic_fetch_add(&f->pos, 1)) < f->l)
		find_solutions(1, f + 1, 0, 0, solution, f->s[pos]);

	// Solve after skipping least frequent set
	f++;
	while ((pos = atomic_fetch_add(&f->pos, 1)) < f->l)
		find_solutions(1, f + 1, 0, 1, solution, f->s[pos]);

	atomic_fetch_add(&solvers_done, 1);
} // solve_work

void
solve()
{
	// Instruct any waiting worker-threads to start solving
	solve_go = 1;

	// The main thread also participates in finding solutions
	solve_work();

	// Wait for all solver threads to finish up
	while(solvers_done < nthreads)
		usleep(1);
} // solve

void *
work_pool(void *arg)
{
	struct worker *work = (struct worker *)arg;
	int wn = work - workers;

	if (pthread_detach(pthread_self()))
		perror("pthread_detach");

	while (read_go == 0)
		usleep(1);

	if (wn < max_readers)
		file_reader(work);
	
	while (solve_go == 0)
		usleep(1);

	solve_work();

	return NULL;
} // work_pool


// ********************* MAIN SETUP AND OUTPUT ********************

int
main(int argc, char *argv[])
{
	struct timespec t0[1], t1[1], t2[1], t3[1], t4[1], t5[1];
	char file[256];
	pthread_t tid[1];

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

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t0);

	for (int i = 1; i < nthreads; i++)
		pthread_create(tid, NULL, work_pool, workers + i);

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t1);

	read_words(file);

	if (write_metrics) clock_gettime(CLOCK_MONOTONIC, t2);

	setup_frequency_sets();

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
		printf("%c set_length = %d\n", c, frq[i].l);
	}
	printf("\n\n");

	printf("Num Unique Words  = %8d\n", nkeys);
	printf("Hash Collisions   = %8u\n", hash_collisions);
	printf("Number of threads = %8d\n", nthreads);

	printf("\nNUM SOLUTIONS = %d\n", num_sol);

	printf("\nTIMES TAKEN :\n");
	print_time_taken("Total", t0, t5);
	printf("\n");
	print_time_taken("Spawn Workers", t0, t1);
	print_time_taken("File Load", t1, t2);
	print_time_taken("Frequency Set Build", t2, t3);
	print_time_taken("Main Algorithm", t3, t4);
	print_time_taken("Emit Results", t4, t5);

	exit(0);
} // main
