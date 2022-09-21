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

// NUM_POISON must be defined before include utilities.h
#define	NUM_POISON 0

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
//		assert(wd != NULL);

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
	for (uint32_t *set, *end; f < e; f++) {
		if (mask & f->m)
			continue;

		CALCULATE_SET_AND_END;

		uint32_t ks[1024], *kp = ks;

		while (set < end)
			kp += !((*kp = *set++) & mask);

		for (*kp = 0, kp = ks; (key = *kp++); )
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
	while ((pos = atomic_fetch_add(&set0pos, 1)) < t->l)
		find_solutions(1, f + 1, solution, 0, t->s[pos], 0);

	// Solve after skipping least frequent set
	f++;
	t = f->sets;
	while ((pos = atomic_fetch_add(&set1pos, 1)) < t->l)
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
