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

static int	write_metrics = 0;
static int	nthreads = 0;
static int	nkeys = 0;

// Keep frequently modified atomic variables on their own CPU cache line
atomic_int	num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done = 0;
atomic_int	solvers_synced = 0;
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


#if 0
uint32_t num_combos = 0;
int8_t	combos[14950][4];
int8_t cidx[4];
void
combo(int n, int r, int d)
{
	if (r == 0) {
		int8_t *cmb = combos[num_combos++];

		cmb[0] = cidx[3];
		cmb[1] = cidx[2];
		cmb[2] = cidx[1];
		cmb[3] = cidx[0];
		return;
	}

	for (int end = n - r; d <= end; d++) {
		cidx[r-1] = d;
		combo(n, r - 1, d + 1);
	}
} //combo
#endif


uint32_t *dset = NULL, *zset = NULL;

void
setup_frequency_sets()
{
	register struct frequency *f = &frq[0];
	register uint32_t *kp = keys, mask, *ks, key;

	qsort(f, 26, sizeof(*f), by_frequency_hi);

	mask = f->m;
	f->s = kp;
	for (ks = kp; (key = *ks); ks++) {
		if (key & mask) {
			*ks = *kp;
			*kp++ = key;
		}
	}
	f->e = kp;
	f->l = kp - f->s;

	// 0-terminate this frequency key set
	*ks++ = *kp;
	*kp++ = 0;

	// Ensure key set is 0 terminated for next loop
	*ks = 0;

	frq[1].s = kp;
	frq[1].e = ks;
	frq[1].l = ks - kp;

	zset = keys;
	dset = kp;
} // setup_frequency_sets


// ********************* SOLUTION FUNCTIONS ********************
static inline void
add_solution_real(uint32_t *solution)
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
} // add_solution_real

static inline void
add_solution(uint32_t key0, uint32_t key1, uint32_t key2, uint32_t key3, uint32_t key4)
{
	uint32_t solution[6];

	solution[1] = key0;
	solution[2] = key1;
	solution[3] = key2;
	solution[4] = key3;
	solution[5] = key4;
	add_solution_real(solution);
} // add_solution

 
#define FOUR_READY  ((uint32_t)0xEAD1EAD1)

atomic_int num_four = 0;
uint32_t   fourset[6291456] = {0};

static inline void
add_fourset(uint32_t key0, uint32_t key1, uint32_t key2, uint32_t key3)
{
	register uint32_t mask = (key0 | key1 | key2 | key3);
	register uint32_t *f = fourset + (atomic_fetch_add(&num_four, 1) * 6);

	*f++ = mask; *f++ = key0; *f++ = key1; *f++ = key2; *f++ = key3;
	*f = FOUR_READY;
} // add_solution


// s0-s3 = starts of the scanbuf
// s0, s1 remain fixed
// s2, s3 can change depending upon what is found
static inline void
gen_four_set(register uint32_t *s0, register uint32_t *s1, register uint32_t key0)
{
	register uint32_t *k1 = s0, scan, key1;

	while ((key1 = *k1++)) {
		register uint32_t *k2 = k1, *s2 = s1, key2;

		while ((scan = *k2++))
			if (!(key1 & scan))
				*s2++ = scan;

		*s2++ = 0, k2 = s1;
		while ((key2 = *k2++)) {
			register uint32_t *k3 = k2, *s3 = s2, key3;

			while ((scan = *k3++))
				if (!(key2 & scan))
					*s3++ = scan;

			*s3++ = 0, k3 = s2;
			while ((key3 = *k3++)) {
				register uint32_t *k4 = k3;

				while ((scan = *k4++))
					if (!(key3 & scan))
						add_solution(key0, key1, key2, key3, scan);

				add_fourset(key0, key1, key2, key3);
			}
		}
	}
} // gen_four_set

// Top level driver
volatile atomic_int	four_pos = 0;
volatile atomic_size_t	spins = 0;
volatile atomic_int	driver_pos = 0;

// A gormless attempt at a lockless way to process fourset
uint32_t
apply_four()
{
	register uint32_t *zp = zset, key, *fp = fourset, pos, *f, *z;

	do {
		while (four_pos >= num_four) {
			if (solvers_synced >= (nthreads - 1))
				return four_pos;
			atomic_fetch_add(&spins, 1);
		}

		pos = atomic_fetch_add(&four_pos, 1);
		f = fp + (pos * 6);
		z = zp;

		// We over-shot.  Check if all threads are done
		while (pos >= num_four) {
			if (solvers_synced >= (nthreads - 1))
				return pos;
			atomic_fetch_add(&spins, 1);
		}

		// Spin-wait until our data is ready
		while (*(f + 5) != FOUR_READY)
			atomic_fetch_add(&spins, 1);

		while ((key = *z++))
			if (!(key & *f))
				add_solution(key, f[1], f[2], f[3], f[4]);
	} while (1);
} // apply_four


// driver is written to allow for either threaded or non-threaded use
void *
driver(void *arg)
{
	struct worker *work = (struct worker *)arg;
	uint32_t scanbuf[4096];
	register uint32_t *dp = dset, *sp = scanbuf;

	if (work->tid)
		if (pthread_detach(work->tid))
			perror("pthread_detach");

	for (;;) {
		register uint32_t *d = dp + atomic_fetch_add(&driver_pos, 1);
		register uint32_t *s = sp;
		register uint32_t key, scan;

		// Build the scanbuf at this level
		if (!(key = *d++))
			break;
		while ((scan = *d++))
			if (!(key & scan))
				*s++ = scan;
		*s++ = 0;
		if (*scanbuf) {
			madvise(scanbuf, sizeof(scanbuf), MADV_WILLNEED);
			gen_four_set(sp, s, key);
		}
	}

	if (arg) {
		atomic_fetch_add(&solvers_synced, 1);
	}

	apply_four();

	if (arg) {
		atomic_fetch_add(&solvers_done, 1);
	}
	return NULL;
} // driver


void
solve()
{
	for (int i = 1; i < nthreads; i++)
		pthread_create(&workers[i].tid, NULL, driver, workers + i);

	// We have to solve ourselves if we're the only thread
	if (nthreads < 2) {
		workers[0].tid = 0;
		driver(workers);
	}

	// Process fourset while waiting
	apply_four();

	while(solvers_done < (nthreads - 1))
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

	printf("Num Unique Words    = %8d\n", nkeys);
	printf("Hash Collisions     = %8u\n", hash_collisions);
	printf("Number of threads   = %8d\n", nthreads);
	printf("Number of Four Sets = %8d\n", num_four);
	printf("Zero Set Size       = %8d\n", frq[0].l);
	printf("Driver Set Size     = %8d\n", frq[1].l);

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
