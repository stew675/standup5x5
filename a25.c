// A solution to the Stand Up Maths Unique 5x5 word problem
//
// Author: Stew Forster (stew675@gmail.com)	Date: Aug 2022

#include <sys/mman.h>
#include <sys/stat.h>
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

#define	HASHSZ	75011		// Prime number of 75K hash entries

#define MAX_THREADS	64

static int nthreads = 0;

#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)

// A very simple 5 letter word copy
static inline void
wcp(register char *a, register const char *b)
{
	*a++ = *b++; *a++ = *b++; *a++ = *b++; *a++ = *b++; *a = *b;
} // wcp


// A very simple for-purpose hash map implementation
struct wh {
	uint32_t	key;
	char		wd[6];
};

struct wh hashmap[HASHSZ] = {0};
uint32_t collisions = 0;

#define get_hash(x)	(((x) * 3) % HASHSZ)
int
hash_insert(uint32_t key, const char *wd)
{
	key &= 0x3ffffff;

	if (key == 0)
		return 0;

	register struct wh *h = hashmap + get_hash(key), *e = hashmap + HASHSZ;
	register int32_t col = 0;
	do {
		// Check if duplicate key
		if (h->key == key)
			return 0;

		// Check if we can insert at this position
		if (h->key == 0)
			break;

		col++;
		h++;

		if (col == HASHSZ)
			assert(col < HASHSZ);

		if (h == e)
			h -= HASHSZ;
	} while (1);

	collisions += col;

	// Now insert at hash location
	h->key = key;
	wcp(h->wd, wd);

	return 1;
} // hash_insert


char *
hash_lookup(uint32_t key)
{
	key &= 0x3ffffff;

	if (key == 0)
		return NULL;

	register struct wh *h = hashmap + get_hash(key), *e = hashmap + HASHSZ;
	register int32_t col = 0;
	do {
		// Check the not-in-hash scenario
		if (h->key == 0)
			return NULL;

		// Check if entry matches the key
		if (h->key == key)
			break;

		col++;
		h++;

		if (col == HASHSZ)
			assert(col < HASHSZ);

		if (h == e)
			h -= HASHSZ;
	} while (1);

	collisions += col;
	return h->wd;
} // hash_lookup


static inline uint32_t
calc_key(register const char *wd)
{
	register uint32_t key = 0;
	key  = (1 << (*wd++ - 'a'));
	key |= (1 << (*wd++ - 'a'));
	key |= (1 << (*wd++ - 'a'));
	key |= (1 << (*wd++ - 'a'));
	key |= (1 << (*wd - 'a'));

	return key;
} // calc_key


static inline uint32_t
add_word_to_set(register const char *wd)
{
	// Reject word if it has duplicate characters
	for (int p = 0; p < 4; p++)
		for (int r = p + 1; r < 5; r++)
			if (wd[p] == wd[r])
				return 0;

	uint32_t key = calc_key(wd);

	// Returns 0 (false) if key already exists. Since
	// the key is a bitmap and order independent this
	// has the effect of quickly eliminating anagrams
	if (hash_insert(key, wd))
		return key;

	return 0;
} // add_word_to_set


uint32_t zset[4096] = {0};
uint32_t dset[4096] = {0};
uint32_t nkeys = 0;

void
read_words(char *path)
{
	register uint32_t *d = dset, *z = zset, nk = 0;

	struct stat statbuf;

	int fd = open(path, O_RDONLY);
	if (fd < 0)
		handle_error("open");
	if (fstat(fd, &statbuf) < 0)
		handle_error("fstat");

	int len = statbuf.st_size;

	char *addr = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
	if (addr == MAP_FAILED)
		handle_error("mmap");

	// Safe to close file now.  mapping remains until munmap() is called
	close(fd);

	for (register char *s = addr, *e = addr + len; s < e; s++) {
		register char *w = s, c;

		c = *s;   if (c < 'a' || c > 'z') continue;
		c = *++s; if (c < 'a' || c > 'z') continue;
		c = *++s; if (c < 'a' || c > 'z') continue;
		c = *++s; if (c < 'a' || c > 'z') continue;
		c = *++s; if (c < 'a' || c > 'z') continue;
		c = *++s;

		if (c < 'a' || c > 'z') {
			// We have a 5 letter word
			uint32_t key = add_word_to_set(w);

			if (key) {
				if (key & 0x1) {	// a
					*z++ = key;
				} else {
					*d++ = key;
				}
				nk++;
			}
		} else {
			// Advance to next non [a-z] character
			for (c = *++s; c >= 'a' && c <= 'z'; c = *++s);
		}
	}

	// Release our file mapping
	munmap(addr, len);

	nkeys = nk;
} // read_words


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


#define FOUR_READY  ((uint32_t)0xEAD1EAD1)
atomic_int num_sol = 0;
uint32_t   solutions[2048][5];

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

static inline void
add_solution(uint32_t key0, uint32_t key1, uint32_t key2, uint32_t key3, uint32_t key4)
{
	register uint32_t *s = solutions[atomic_fetch_add(&num_sol, 1)];

	*s++ = key0; *s++ = key1; *s++ = key2; *s++ = key3; *s = key4;
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
pthread_t	tids[MAX_THREADS];
volatile atomic_int	four_pos = 0;
volatile atomic_size_t	spins = 0;
atomic_int	driver_pos = 0;
atomic_int	threads_synced = 0;
atomic_int	threads_done = 0;

// A gormless attenpt at a lockless way to process fourset
uint32_t
apply_four()
{
	register uint32_t *zp = zset, key, *fp = fourset, pos, *f, *z;

	do {
		while (four_pos >= num_four) {
			if (threads_synced >= nthreads)
				return four_pos;
			atomic_fetch_add(&spins, 1);
		}

		pos = atomic_fetch_add(&four_pos, 1);
		f = fp + (pos * 6);
		z = zp;

		// We over-shot.  Check if all threads are done
		while (pos >= num_four) {
			if (threads_synced >= nthreads)
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
	uint32_t scanbuf[4096];
	register uint32_t *dp = dset, *sp = scanbuf;

	if (arg) {
		pthread_t tid = *(pthread_t *)arg;
		if (pthread_detach(tid))
			perror("pthread_detach");
	}

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
		atomic_fetch_add(&threads_synced, 1);
	}

	apply_four();

	if (arg) {
		atomic_fetch_add(&threads_done, 1);
	}
	return NULL;
} // driver


void
solve(int threads)
{
	if (threads > MAX_THREADS)
		threads = MAX_THREADS;
	if (threads < 0)
		threads = 0;

	if (threads == 0) {
		driver(NULL);
		return;
	}

	nthreads = threads;
	for (uintptr_t i = 0; i < nthreads; i++) {
		void *arg = (void *)(tids + i);
		pthread_create(tids + i, NULL, driver, arg);
	}

	// Process fourset while waiting
	apply_four();

	while(threads_done < nthreads)
		usleep(1);
} // solve


void
emit_solutions()
{
	FILE *fp = fopen("solutions.txt", "w");
	if (fp == NULL)
		return;

	fprintf(fp, "SOLUTIONS:\n\n");
	for (int i = 0; i < num_sol; i++) {
		fprintf(fp, "%.5s ", hash_lookup(solutions[i][0]));
		fprintf(fp, "%.5s ", hash_lookup(solutions[i][1]));
		fprintf(fp, "%.5s ", hash_lookup(solutions[i][2]));
		fprintf(fp, "%.5s ", hash_lookup(solutions[i][3]));
		fprintf(fp, "%.5s\n", hash_lookup(solutions[i][4]));
	}
	fprintf(fp, "\n");

	fclose(fp);
} // emit_solutions


int
main(int argc, char *argv[])
{
	struct timespec ts[1], te[1];
	int64_t time_taken = 1000000000LL;

	read_words("words_alpha.txt");
//	read_words("nyt_wordle.txt");

	clock_gettime(CLOCK_MONOTONIC, ts);
	solve(14);
	clock_gettime(CLOCK_MONOTONIC, te);

	time_taken *= (te->tv_sec - ts->tv_sec);
	time_taken += (te->tv_nsec - ts->tv_nsec);

	emit_solutions();

	printf("Num Solutions = %d\n", num_sol);
	printf("Num Four Sets = %d\n", num_four);
	printf("Num Unique Word = %d\n", nkeys);
	printf("Hash Collisions = %u\n", collisions);
	printf("Number of threads = %d\n", nthreads);

	printf("Time Taken = %ld.%06lus\n", time_taken / 1000000000, (time_taken % 1000000000) / 1000);
	if (nthreads)
		printf("Number of Thread Contention Busy Spins = %lu\n", spins);

	exit(0);
} // main
