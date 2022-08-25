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

#define	MAX_SOLUTIONS	8192
#define	MAX_WORDS	8192
#define	MAX_THREADS	  64
#define	MAX_READERS	  11	// Virtual systems don't like too many readers
#define	READ_CHUNK	8192
//#define	HASHSZ		69001	// Also a good value
#define	HASHSZ		39009

static const char	*solution_filename = "solutions.txt";

// Worker thread state
static struct worker {
	pthread_t tid;
	char *start;
	char *end;
} workers[MAX_THREADS] __attribute__ ((aligned(64)));

// Word Hash Entries
static uint32_t	hash_collisions	__attribute__ ((aligned(64))) = 0;
static struct wordhash {
	uint32_t	key;
	uint32_t	pos;
} hashmap[HASHSZ] __attribute__ ((aligned(64)));

// Character frequency recording
static struct frequency {
	uint32_t  *s;		// Pointer to set
	uint32_t  *e;		// Pointer to end of set
	uint32_t   m;		// Mask (1 << (c - 'a'))
	int32_t    f;		// Frequency
	int32_t    l;		// Length of set
	atomic_int pos;		// Position within a set
} frq[26] __attribute__ ((aligned(64)));

static int32_t	min_search_depth __attribute__ ((aligned(64))) = 0;
static int	write_metrics = 0;
static int	nthreads = 0;
static int	nkeys = 0;

// Keep frequently modified atomic variables on their own CPU cache line
atomic_int	num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done = 0;
atomic_int	solvers_done = 0;

// Allow for up to 3.2x the number of unique non-anagram words
static char	words[MAX_WORDS * 16] __attribute__ ((aligned(64)));

// We build the solutions directly as a character array to write out when done
static char	solutions[MAX_SOLUTIONS * 30] __attribute__ ((aligned(64)));

// We add 1024 here to MAX_WORDS to give us extra space to perform vector
// alignments for the AVX functions.  At the very least the keys array must
// be 32-byte aligned, but we align it to a typical system page boundary
static uint32_t	keys[MAX_WORDS + 1024] __attribute__ ((aligned(4096)));

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
		register uint32_t mask = f->m, *ks, key;

		f->s = kp;
		for (ks = kp; (key = *ks); ks++) {
			if (key & mask) {
				*ks = *kp;
				*kp++ = key;
			}
		}
		f->e = kp;

		register uint32_t nfk = kp - f->s;
		f->l = nfk;

		// Update the min_search_depth if needed
		if (nfk > 0)
			min_search_depth = i - 3;

		// Ensure 32-byte alignments for the set
		// We "poison" the alignment values with all bits set
		while ((nfk++ % 16)) {
			*ks++ = *kp;
			*kp++ = (uint32_t)(~0);
		}

		// Ensure key set is 0 terminated for next loop
		*ks = 0;
	}
} // setup_frequency_sets


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

#define ctz(x)		__builtin_ctz(x)
#define vzero		_mm512_setzero_si512()
#define vmask(x)	_mm512_set1_epi32(x)
#define vkeys(y) 	_mm512_load_si512((__m512i *)y)
#define vand(x, y)	_mm512_and_si512(vmask(x), vkeys(y))
#define vscan(x, y)	_mm512_cmpeq_epi32_mask(vand(x, y), vzero)

// find_solutions() which is the busiest loop is kept
// as small and tight as possible for the most speed
void
find_solutions(register int depth, register int setnum, register uint32_t mask,
		register int skipped, register uint32_t *solution, register uint32_t key)
{
	solution[depth] = key;
	mask |= key;
	if (depth == 5)
		return add_solution(solution);

	for (register int e = min_search_depth + depth; setnum < e; setnum++)
		if (!(mask & frq[setnum].m)) {
			for (register uint32_t *set = frq[setnum].s, *end = frq[setnum].e; set < end; set += 16)
				for (register uint16_t vresmask = vscan(mask, set); vresmask; vresmask &= (vresmask - 1))
					find_solutions(depth + 1, setnum + 1, mask, skipped, solution, set[ctz(vresmask)]);
			if (skipped) return;
			skipped = 1;
		}
} // find_solutions

// Top level solver
void *
solve_work(void *arg)
{
	struct worker *work = (struct worker *)arg;

	if (work->tid)
		if (pthread_detach(work->tid))
			perror("pthread_detach");

	// Initial work solver.  Here we overload the use
	// of value of setnum as our skipped argument too
	for (int setnum = 0; setnum < 2; setnum++) {
		register uint32_t pos, key;
		uint32_t solution[6];
		atomic_int *apos = &frq[setnum].pos;

		do {
			pos = atomic_fetch_add(apos, 1);

			if (pos >= frq[setnum].l)
				break;

			key = frq[setnum].s[pos];
			solution[0] = key;
			find_solutions(1, setnum + 1, key, setnum, solution, key);
		} while (1);
	}

	atomic_fetch_add(&solvers_done, 1);
	return NULL;
} // solve_work

void
solve()
{
	for (int i = 1; i < nthreads; i++)
		pthread_create(&workers[i].tid, NULL, solve_work, workers + i);

	// The main thread also participates in finding solutions
	workers[0].tid = 0;
	solve_work(workers);

	// Wait for all solver threads to finish up
	while(solvers_done < nthreads)
		usleep(1);
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
