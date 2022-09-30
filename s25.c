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
find_skipped(struct frequency *f, uint32_t mask, uint32_t *sp)
{
	if (__builtin_popcount(mask) == 25)
		return add_solution(sp - 5);

	uint32_t ks[1024] __attribute__((aligned(64)));
	uint32_t key, *set, *end, *kp = ks;

	while (mask & f->m) f++;

	CALCULATE_SET_AND_END;

	while (set < end)
		kp += !((*kp = *set++) & mask);

	for (*kp = 0, kp = ks; (*sp = key = *kp++); )
		find_skipped(f + 1,  mask | key, sp + 1);
} // find_solutions

// Since find_solutions() is the busiest function we keep the loops
// within it as small and tight as possible for the most speed
void
find_solutions(struct frequency *f, uint32_t mask, uint32_t *sp)
{
	if (__builtin_popcount(mask) == 25)
		return add_solution(sp - 5);

	uint32_t ks[1024] __attribute__((aligned(64)));
	uint32_t key, *set, *end, *kp = ks;

	while (mask & f->m) f++;

	CALCULATE_SET_AND_END;

	while (set < end)
		kp += !((*kp = *set++) & mask);

	for (*kp = 0, kp = ks; (*sp = key = *kp++); )
		find_solutions(f + 1, mask | key, sp + 1);

	find_skipped(f + 1, mask, sp);
} // find_solutions

// Thread driver
void
solve_work()
{
	uint32_t solution[6] __attribute__((aligned(64)));
	struct tier *t;
	int32_t pos, len;

	// Solve starting with least frequent set
	t = frq[0].tiers;
	len = t->end - t->set;
	while ((pos = atomic_fetch_add(&set0pos, 1)) < len)
		find_solutions(frq + 1, (*solution = t->set[pos]), solution);

	// Solve after skipping least frequent set
	t = frq[1].tiers;
	len = t->end - t->set;
	while ((pos = atomic_fetch_add(&set1pos, 1)) < len)
		find_skipped(frq + 2, (*solution = t->set[pos]), solution);

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
