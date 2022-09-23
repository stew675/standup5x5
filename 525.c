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

// NUM_POISON must be defined before include utilities.h
#define	NUM_POISON	16
#define USE_AVX2_SCAN

#include "utilities.h"

// ********************* SOLVER ALGORITHM ********************

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
		int skipped, uint32_t *sp, uint32_t key)
{
	*sp++ = key;
	if (depth == 5)
		return add_solution(sp - 5);
	mask |= key;

	for (uint32_t *set, *end; ; f++) {
		if (mask & f->m)
			continue;

		CALCULATE_SET_AND_END;

		for (; set < end; set += 16)
			for (uint16_t vresmask = vscan(mask, set); vresmask; vresmask &= vresmask - 1)
				find_solutions(depth + 1, f + 1, mask, skipped, sp, set[__builtin_ctz(vresmask)]);

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
	struct tier *t;
	int32_t pos;

	// Solve starting with least frequent set
	t = frq[0].sets;
	while ((pos = atomic_fetch_add(&set0pos, 1)) < t->l)
		find_solutions(1, frq + 1, 0, 0, solution, t->s[pos]);

	// Solve after skipping least frequent set
	t = frq[1].sets;
	while ((pos = atomic_fetch_add(&set1pos, 1)) < t->l)
		find_solutions(1, frq + 2, 0, 1, solution, t->s[pos]);

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

#else

int
main()
{
	fprintf(stderr, "AVX-512 not supported on this system\n");
	exit(1);
} // main
#endif
