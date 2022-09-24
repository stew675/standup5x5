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
add_solution(uint32_t *sp)
{
	char *so = solutions + (atomic_fetch_add(&num_sol, 1) << 5);

	*(uint64_t *)so = *(uint64_t *)hash_lookup(*sp++);
	so[5] = '\t'; so += 6;

	*(uint64_t *)so = *(uint64_t *)hash_lookup(*sp++);
	so[5] = '\t'; so += 6;

	*(uint64_t *)so = *(uint64_t *)hash_lookup(*sp++);
	so[5] = '\t'; so += 6;

	*(uint64_t *)so = *(uint64_t *)hash_lookup(*sp++);
	so[5] = '\t'; so += 6;

	*(uint64_t *)so = *(uint64_t *)hash_lookup(*sp);
	so[5] = ' '; so[6] = ' '; so[7] = '\n';
} // add_solution

void
find_skipped(uint32_t mask, uint32_t *sp)
{
	if (__builtin_popcount(mask) == 26)
		return add_solution(sp - 4);

	struct frequency *f = frq + __builtin_ctz(~mask);
	uint32_t ks[1024] __attribute__((aligned(64)));
	uint32_t key, *set, *end, *kp = ks;

	CALCULATE_SET_AND_END;

	while (set < end)
		kp += !((*kp = *set++) & mask);

	for (sp++, *kp = 0, kp = ks; (*sp = key = *kp++); )
		find_skipped(mask | key, sp);
} // find_solutions

// Since find_solutions() is the busiest function we keep the loops
// within it as small and tight as possible for the most speed
void
find_solutions(uint32_t mask, uint32_t *sp)
{
	if (__builtin_popcount(mask) == 25)
		return add_solution(sp - 4);

	struct frequency *f = frq + __builtin_ctz(~mask);
	uint32_t ks[1024] __attribute__((aligned(64)));
	uint32_t key, *set, *end, *kp = ks;

	CALCULATE_SET_AND_END;

	while (set < end)
		kp += !((*kp = *set++) & mask);

	for (sp++, *kp = 0, kp = ks; (*sp = key = *kp++); )
		find_solutions(mask | key, sp);

	find_skipped(mask | f->m, sp - 1);
} // find_solutions

// Thread driver
void
solve_work()
{
	uint32_t solution[6] __attribute__((aligned(64)));
	struct tier *t;
	int32_t pos;

	// Solve starting with least frequent set
	t = frq[0].sets;
	while ((pos = atomic_fetch_add(&set0pos, 1)) < t->l)
		find_solutions((*solution = t->s[pos]), solution);

	// Solve after skipping least frequent set
	t = frq[1].sets;
	while ((pos = atomic_fetch_add(&set1pos, 1)) < t->l)
		find_skipped((*solution = t->s[pos]) | frq[0].m, solution);

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
