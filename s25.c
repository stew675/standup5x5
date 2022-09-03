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

#include "utilities.h"

// ********************* SOLUTION FUNCTIONS ********************

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

// Since find_solutions() is the busiest function we keep the loops
// within it as small and tight as possible for the most speed
void
find_solutions(int depth, struct frequency *f, uint32_t *solution,
		uint32_t mask, uint32_t key, int skipped)
{
	solution[depth] = key;
	if (depth == 5)
		return add_solution(solution);
	mask |= key;

	struct frequency *e = frq + (min_search_depth + depth);
	for (; f < e; f++) {
		if (mask & f->m)
			continue;

		struct tier *t = f->sets + !!(mask & f->tm1) + (2 * !!(mask & f->tm2)) + (4 * !!(mask & f->tm3));

		// Determine the values for set and end
		// The !! means we end up with only 0 or 1
		int mf = !!(mask & f->tm4);	// Mask First
		int ms = !!(mask & f->tm5);	// Mask Second

		// A branchless calculation of end
		uint32_t *end = t->s + (ms * t->toff3) + (!ms * t->l);

		// A branchless calculation of set
		ms &= !mf;
		uint32_t *set = t->s + ((mf & !ms) * t->toff2) + (ms * t->toff1);

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
	struct frequency *f = frq;
	struct tier *t = f->sets;
	int32_t pos;

	// Solve starting with least frequent set
	while ((pos = atomic_fetch_add(&f->pos, 1)) < t->l)
		find_solutions(1, f + 1, solution, 0, t->s[pos], 0);

	// Solve after skipping least frequent set
	f++;
	t = f->sets;
	while ((pos = atomic_fetch_add(&f->pos, 1)) < t->l)
		find_solutions(1, f + 1, solution, 0, t->s[pos], 1);

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

	setup_frequency_sets(1);

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
