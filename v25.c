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

// NUM_POISON must be defined before include utilities.h
#define NUM_POISON 8
#define USE_AVX2_SCAN

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

#define vzero _mm256_setzero_si256()

static inline uint32_t *
vscan(uint32_t mask, uint32_t **set, uint32_t *to)
{
	// Find all valid keys
	__m256i vmask = _mm256_set1_epi32(mask);
	__m256i vkeys = _mm256_loadu_si256((__m256i *)*set);
	__m256i vres = _mm256_cmpeq_epi32(_mm256_and_si256(vmask, vkeys), vzero);

	// Blend results into destination.  Store everything valid to destination, or store a zero
	_mm256_storeu_si256((__m256i *)to, _mm256_blendv_epi8(vzero, vkeys, vres));
	*set += 8;

	// Pack the results (remove the zeros)
	// XXX Is there a better way to do this?
	for (uint32_t *ts = to, i = 8; i--; )
		to += !!(*to = *ts++);

	return to;
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
	for (uint32_t *set, *end; f < e; f++) {
		if (mask & f->m)
			continue;

		CALCULATE_SET_AND_END;

		uint32_t ks[1024], *kp = ks;

		// Find all matching keys
		while (set < end)
			kp = vscan(mask, &set, kp);

		// Now make the recursion calls using the packed results
		for (*kp = 0, kp = ks; (key = *kp++); )
			find_solutions(depth + 1, f + 1, mask, skipped, solution, key);

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
	struct tier *t;
	int32_t pos;

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
