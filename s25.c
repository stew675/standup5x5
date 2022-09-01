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

#define	MAX_SOLUTIONS	8192
#define	MAX_WORDS	8192
#define	MAX_THREADS	  64

static const char	*solution_filename = "solutions.txt";

// Worker thread state
static struct worker {
	char     *start;
	char     *end;
} workers[MAX_THREADS] __attribute__ ((aligned(64)));

// Character frequency recording
static struct frequency {
	uint32_t  *s;		// Pointer to set
	uint32_t   m;		// Mask (1 << (c - 'a'))
	uint32_t  tm1;		// Tier 1 Mask
	uint32_t  tm2;		// Tier 2 Mask
	uint32_t  toff1;	// Tier offset 1
	uint32_t  toff2;	// Tier offset 2
	uint32_t  toff3;	// Tier offset 3
	int32_t    f;		// Frequency
	int32_t    l;		// Length of set
	atomic_int pos;		// Position within a set
	uint32_t  pad[5];	// Pad to 64 bytes
} frq[26] __attribute__ ((aligned(64)));

// Keep frequently modified atomic variables on their own CPU cache line
atomic_int 	num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done = 0;
atomic_int	solvers_done = 0;

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
	register uint32_t *kp = keys, tm1, tm2;

	qsort(f, 26, sizeof(*f), by_frequency_lo);
	tm1 = frq[25].m;
	tm2 = frq[24].m;

	// Split full keyset into sets organised by the
	// least frequently occurring letter to the most
	for (int i = 0; i < 26; i++, f++) {
		register uint32_t mask, *ks, key;

		if (i == 7)
			rescan_frequencies(i, kp);

		f->tm1 = tm1;
		f->tm2 = tm2;

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

		// 0-terminate this frequency key set
		*ks++ = *kp;
		*kp++ = 0;

		// Ensure key set is 0 terminated for next loop
		*ks = 0;
	}

	// Now organise each set into 2 subsets, that which
	// has tm1 followed by that which does not
	f = frq;
	for (int i = 0; i < 26; i++, f++) {
		register uint32_t mask = f->tm1, *ks, key, *end;

		kp = f->s;
		end = f->s + f->l;
		for (ks = kp; ks < end; ks++) {
			key = *ks;
			if (key & mask) {
				*ks = *kp;
				*kp++ = key;
			}
		}
		f->toff2 = kp - f->s;

		// Now organise each first subset into that which
		// has tm2 followed by that which does not, and
		// then each second subset into that which does
		// not have tm2 followed by that which does

		mask = f->tm2;

		// First Subset has tm2 then not
		kp = f->s;
		end = f->s + f->toff2;
		for (ks = kp; ks < end; ks++) {
			key = *ks;
			if (key & mask) {
				*ks = *kp;
				*kp++ = key;
			}
		}
		f->toff1 = kp - f->s;

		// Second Subset does not have tm2 then has
		kp = f->s + f->toff2;
		end = f->s + f->l;
		for (ks = kp; ks < end; ks++) {
			key = *ks;
			if (!(key & mask)) {
				*ks = *kp;
				*kp++ = key;
			}
		}
		f->toff3 = kp - f->s;
	}
} // setup_frequency_sets


// ********************* SOLUTION FUNCTIONS ********************

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

// Since find_solutions() is the busiest function we keep the loops
// within it as small and tight as possible for the most speed
void
find_solutions(register int depth, register struct frequency *f, register uint32_t *solution,
		register uint32_t mask, register uint32_t key, register int skipped)
{
	solution[depth] = key;
	if (depth == 5)
		return add_solution(solution);
	mask |= key;

	register struct frequency *e = frq + (min_search_depth + depth);
	for (; f < e; f++) {
		if (mask & f->m)
			continue;

		// Determine the values for set and end
		// The !! means we end up with only 0 or 1
		register int mt1 = !!(mask & f->tm1);
		register int mt2 = !!(mask & f->tm2);

		// A branchless calculation of end
		register uint32_t *end = f->s + (mt2 * f->toff3) + (!mt2 * f->l);

		// A branchless calculation of set
		mt2 &= !mt1;
		register uint32_t *set = f->s + ((mt1 & !mt2) * f->toff2) + (mt2 * f->toff1);

		while (set < end)
			if (!((key = *set++) & mask))
				find_solutions(depth + 1, f + 1, solution, mask, key, skipped);

		if (skipped)
			return;

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
		find_solutions(1, f + 1, solution, 0, f->s[pos], 0);

	// Solve after skipping least frequent set
	f++;
	while ((pos = atomic_fetch_add(&f->pos, 1)) < f->l)
		find_solutions(1, f + 1, solution, 0, f->s[pos], 1);

	atomic_fetch_add(&solvers_done, 1);
} // solve_work

void
solve()
{
	// Instruct worker pool to start solving
	start_solvers();

	// Do work ourselves too!
	solve_work();

	// Wait for any other threads to finish up
	while(solvers_done < nthreads)
		usleep(1);
} // solve

int
main(int argc, char *argv[])
{
	struct timespec t1[1], t2[1], t3[1], t4[1], t5[1];
	char file[256];

	// Copy in a default file-name
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
		printf("%c set_length=%4d  toff1=%4d, toff2=%4d, toff[3]=%4d\n",
			c, frq[i].l, frq[i].toff1, frq[i].toff2, frq[i].toff3);
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
