#include <immintrin.h>

#define HASHBITS              15
#define MAX_SOLUTIONS       8192
#define MAX_WORDS           8192
#define MAX_THREADS           16
#define MAX_READERS            8	// No more than 8 ever needed

static const char	*solution_filename = "solutions.txt";

// Worker thread state
static struct worker {
	char     *start;
	char     *end;
} workers[MAX_THREADS] __attribute__ ((aligned(64)));

// Set Pointers (8 bytes in size)
static struct tier {
	uint32_t	*set;	// Offset to start of set
	uint32_t	*end;	// Offset to end of set
} tiers[26][64] __attribute__ ((aligned(64)));

// Character frequency recording
static struct frequency {
	// Mask (1 << (c - 'a'))
	uint32_t	m	__attribute__ ((aligned(64)));
	uint32_t	f;		// Frequency

	uint32_t	tm1;		// Tiered Mask 1
	uint32_t	tm2;		// Tiered Mask 2
	uint32_t	tm3;		// Tiered Mask 3
	uint32_t	tm4;		// Tiered Mask 4
	uint32_t	tm5;		// Tiered Mask 5
	uint32_t	tm6;		// Tiered Mask 6

	struct tier	*tiers;		// The set tiers

	uint32_t	tmm;		// Logical OR of tm1..tm4

	uint16_t	b;		// char - 'a'
	uint16_t	ready;		// Ready to set up

} frq[26] __attribute__ ((aligned(64)));

// Keep atomic variables on their own CPU cache line
atomic_int 	num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	setup_set	__attribute__ ((aligned(64))) = 0;
atomic_int	setups_done	__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done	__attribute__ ((aligned(64))) = 0;
atomic_int	solvers_done	__attribute__ ((aligned(64))) = 0;
atomic_int	first_rdr_done	__attribute__ ((aligned(64))) = 0;
atomic_int	set0pos		__attribute__ ((aligned(64))) = 0;
atomic_int	set1pos		__attribute__ ((aligned(64))) = 0;

// Put volatile thread sync variables on their own CPU cache line
static volatile int	workers_start	__attribute__ ((aligned(64))) = 0;
static volatile int	go_solve	__attribute__ ((aligned(64))) = 0;
static volatile int	num_readers	__attribute__ ((aligned(64))) = 0;

// Put all general global variables together on their own CPU cache line
static uint32_t hash_collisions __attribute__ ((aligned(64))) = 0;
static int	write_metrics = 0;
static int	nthreads = 0;
static int	nkeys = 0;

// We build the solutions directly as a character array to write out when done
static char     solutions[MAX_SOLUTIONS * 32] __attribute__ ((aligned(64)));

// Allow for up to 3x the number of unique non-anagram words
static char     words[MAX_WORDS * 24] __attribute__ ((aligned(64)));
static uint32_t wordkeys[MAX_WORDS * 3] __attribute__ ((aligned(64)));

// We add 1024 here to MAX_WORDS to give us extra space to perform vector
// alignments for the AVX functions.  At the very least the keys array must
// be 32-byte aligned, but we align it to 64 bytes anyway
static	uint32_t	keys[MAX_WORDS + 1024] __attribute__ ((aligned(64)));
static	uint32_t	tkeys[26][MAX_WORDS * 2] __attribute__ ((aligned(64)));
static	uint32_t	unmap[32] __attribute__((aligned(64)));

// Per-reader frequency collation stats.  We set to 32, instead of just 26, to
// ensure readers aren't sharing CPU cache lines (which are 64 bytes wide)
static	uint32_t	cfs[MAX_READERS][32] __attribute__((aligned(64))) = {0};

static void solve();
static void solve_work();
static void set_tier_offsets(struct frequency *f);

void
print_time_taken(char *label, struct timespec *ts, struct timespec *te)
{
	int64_t time_taken = 1000000000LL;	// Number of ns in 1s
	time_taken *= (te->tv_sec - ts->tv_sec);
	time_taken += (te->tv_nsec - ts->tv_nsec);

	printf("%-20s = %ld.%06lus\n", label, time_taken / 1000000000,
				     (time_taken % 1000000000) / 1000);
} // print_time_taken
 
//********************* INIT FUNCTIONS **********************

static void
frq_init()
{
	memset(frq, 0, sizeof(frq));

	for (int b = 0; b < 26; b++) {
		frq[b].tiers = tiers[b];
		frq[b].m = (1UL << b);	// The bit mask
	}
} // frq_init

//********************* UTILITY FUNCTIONS **********************

// Determine number of threads to use
int
get_nthreads()
{
	int ncpus = sysconf(_SC_NPROCESSORS_ONLN);

	if (ncpus < 2)
		return 1;

	if (ncpus < 5)
		return ncpus;

	if (ncpus < 9)
		return ncpus - 1;

	// Generally speaking, not much to be gained beyond 20 threads
	if ((ncpus - 2) > 20)
		return 20;

	return ncpus - 2;
} // get_nthreads

// Given a 5 letter word, calculate the bit-map representation of that word
static inline uint32_t
calc_key(const char *wd)
{
	uint32_t one = 1, mask = 0x1F;
	uint32_t key = (one << (wd[0] & mask)) |
                       (one << (wd[1] & mask)) |
                       (one << (wd[2] & mask)) |
                       (one << (wd[3] & mask)) |
                       (one << (wd[4] & mask));
	return key >> 1;
} // calc_key

//********************* HASH TABLE FUNCTIONS **********************

// A very simple for-purpose hash map implementation.  Used to
// lookup words given the key representation of that word
// I've included 3 decent key_hash() functions here. All should
// work decently for most English 5-letter words @ HASHBITS = 15

#define	HASHSZ          (1 << HASHBITS)
#define HASHMASK        (HASHSZ - 1)
#define key_hash(x)	(((x * 5287) ^ (x >> 11)) & HASHMASK)
//#define key_hash(x)	(((x * 13334) ^ x ^ (x >> 12)) & HASHMASK)
//#define key_hash(x)	(x ^ (x >> 6) ^ (x >> 10) ^ (~x >> 1)) & HASHMASK

// Key Hash Entries
// We keep keys and positions in separate array because faster to initialise
uint32_t keymap[HASHSZ] __attribute__ ((aligned(64)));
uint32_t posmap[HASHSZ] __attribute__ ((aligned(64)));

static void
hash_init()
{
	memset(keymap, 0, sizeof(keymap));
} // hash_init

uint32_t
hash_insert(uint32_t key, uint32_t pos)
{
	uint32_t col = 0, hashpos = key_hash(key);

	do {
		// Check if we can insert at this position
		if (keymap[hashpos] == 0)
			break;

		// Check if duplicate key
		if (keymap[hashpos] == key)
			return 0;

		if (++hashpos == HASHSZ)
			hashpos = 0;

		col++;
	} while (1);

	// Now insert at hash location
	keymap[hashpos] = key;
	posmap[hashpos] = pos << 3;

	hash_collisions += col;

	return 1;
} // hash_insert

const char *
hash_lookup(uint32_t key)
{
	uint32_t hashpos = key_hash(key);

	do {
		// Check for a match
		if (keymap[hashpos] == key)
			break;

		// Check the not-in-hash scenario
		if (keymap[hashpos] == 0)
			return NULL;

		if (++hashpos == HASHSZ)
			hashpos = 0;
	} while (1);

	return words + posmap[hashpos];
} // hash_lookup
#undef key_hash

// Just a handy debugging function which was used when developing the
// 5 letter word extraction bit masking algorithm within find_words()
void
print_bits32(char *label, uint32_t v)
{
	printf("%s ", label);
	for (int i = 0; i < 32; i++) {
		if (v >> 31)
			printf("1");
		else
			printf("0");
		v <<= 1;
	}
	printf("\n");
} // print_bits32

void
print_bits(char *label, uint64_t v)
{
	printf("%s ", label);
	for (int i = 0; i < 64; i++) {
		if (v >> 63)
			printf("1");
		else
			printf("0");
		v <<= 1;
	}
	printf("\n");
} // print_bits

// ********************* FILE READER ********************

#define READ_CHUNK        65536		// Appears to be optimum

void
find_words(char *s, char *e, uint32_t rn)
{
	char *fives[(READ_CHUNK / 6) + 1] __attribute__((aligned(64)));
	char **fivep = fives;
	char a = 'a', z = 'z';
	int64_t msbset = 0x8000000000000000;
	uint32_t *cf = cfs[rn];

#ifdef __AVX2__
	// AVX512 is about 10% faster than AVX2 for processing the words
	// Use AVX512 if the current platform supports it

	// Prepare 2 constant vectors with all a's and all z's
#ifdef __AVX512F__
	__m512i avec = _mm512_set1_epi8(a);
	__m512i zvec = _mm512_set1_epi8(z);
#else
	__m256i avec = _mm256_set1_epi8(a);
	__m256i zvec = _mm256_set1_epi8(z);
#endif

	e -= 64;
	while (s < e) {
#ifdef __AVX512F__
		// Unaligned load of a vector with the next 64 characters
		__m512i wvec = _mm512_loadu_si512((const __m512i_u *)s);

		// Find the lower-case letters in the word vector
		// wmask will have a 0-bit for every lower-case letter in the vector
		uint64_t wmask = _mm512_cmp_epi8_mask(wvec, avec, _MM_CMPINT_LT) |
					  _mm512_cmp_epi8_mask(zvec, wvec, _MM_CMPINT_LT);
#else
		// Here we're in AVX2 mode.  Emulate AVX512 mode by doing 2 loads

		// Unaligned load of 2 vectors with the next 64 characters
		__m256i wvec1 = _mm256_loadu_si256((const __m256i_u *)s);
		__m256i wvec2 = _mm256_loadu_si256((const __m256i_u *)(s + 32));

		// Find the lower-case letters in the word vector
		// wmask will have a 0-bit for every lower-case letter in the vector
		__m256i wres = _mm256_or_si256(_mm256_cmpgt_epi8(avec, wvec1),
						_mm256_cmpgt_epi8(wvec1, zvec));
		uint32_t wmask1 = _mm256_movemask_epi8(wres);

		// Load and process another 32 characters
		wres = _mm256_or_si256(_mm256_cmpgt_epi8(avec, wvec2),
					_mm256_cmpgt_epi8(wvec2, zvec));
		uint64_t wmask = _mm256_movemask_epi8(wres);

		// Merge the results of the two loads
		wmask = (wmask << 32) | wmask1;
#endif
		// Handle lines over 64 characters in length.  Jump ahead just
		// far enough such that we won't accidentally feed the last 5
		// characters from an overly long line into the next pass
		// !wmask is never true for words_alpha.txt, so the CPU branch
		// predictor should never get this wrong
		if (!wmask) {
			s += 58;
			continue;
		}

		// Calculate where to start the next loop pass and invalidate
		// everything after the last non-lower case letter
		char *ns = s + 64;
		uint32_t nlz = __builtin_clzll(wmask);
		ns -= nlz;
		wmask |= (msbset >> nlz);

		// Get the 1's complement of wmask. ocwm will have a 1-bit set
		// for every valid lower-case letter than was in the vector.
		uint64_t ocwm = ~wmask;

		// Isolate all words of <=5 characters
		wmask = (wmask >> 5) & ((wmask << 1) | 1);

		// Reset bit 0 in all words with less than 5 characters
		ocwm = ((ocwm >> 1) & (ocwm >> 2));

		// Intersect the two
		wmask &= (ocwm & (ocwm >> 2));

		// wmask will now contain a 1 bit located at the
		// start of every word with exactly 5 letters

		// Process all 5 letter words in the vector
		while (wmask) {
			// Get a pointer to the start of the 5 letter word
			char *w = s + __builtin_ctzll(wmask);

			// Add word to our list
			*fivep = w;

			// Advance list if word has no duplicate characters
			fivep += (__builtin_popcount(calc_key(w)) == 5);

			// Unset the lowest bit
			wmask &= (wmask - 1);
		}
		s = ns;
	}
	e += 64;
#endif

	// Scalar code to find 5 words. This also
	// handles residuals from the vector loop
	for (char c, *w = s; s < e; w = s) {
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;

		// We've now found 5 [a..z] characters in a row
		c = *s++;
		if ((c < a) || (c > z))  {
			*fivep = w;
			fivep += (__builtin_popcount(calc_key(w)) == 5);
		}

		// Just quickly find the next line
		while (c != '\n')
			c = *s++;
	}

	// Bulk process all found unique 5 words
	// If no words to process, return now
	int num = fivep - fives;
	if (num == 0)
		return;

	// Bulk reserve where to place the words
	int pos = atomic_fetch_add(&num_words, num);
	fivep = fives;
	while (num--) {
		char *w = *fivep++;

		// Copy word to word table as a single 64-bit copy
		*(uint64_t *)(words + (pos << 3)) = *(uint64_t *)w;

		// Copy key to wordkeys array
		uint32_t key = calc_key(w);
		wordkeys[pos++] = key;

		// Get character frequencies
		cf[__builtin_ctz(key)]++; key &= key - 1;
		cf[__builtin_ctz(key)]++; key &= key - 1;
		cf[__builtin_ctz(key)]++; key &= key - 1;
		cf[__builtin_ctz(key)]++; key &= key - 1;
		cf[__builtin_ctz(key)]++;
	}
} // find_words

//#define FILE_READER_TIMES

void
file_reader(struct worker *work)
{
	uint32_t rn = work - workers;
#ifdef FILE_READER_TIMES
	struct timespec t1[1], t2[1];
	clock_gettime(CLOCK_MONOTONIC, t1);
#endif

	// The e = s + (READ_CHUNK + 1) below is done because each reader
	// (except the first) only starts at a newline.  If the reader
	// starts at the very start of a 5 letter word, that means that it
	// would skip that word.  By adding the extra 1 here, the reader
	// processing the chunk before it can catch the word that may
	// have been skipped by the reader ahead of it
	do {
		char *s = work->start;
		s += atomic_fetch_add(&file_pos, READ_CHUNK);
		char *e = s + (READ_CHUNK + 1);

		if (s > work->end)
			break;
		if (e > work->end)
			e = work->end;

		// Make sure to only start after a newline
		// if we are not at the start of the file
		if (s > work->start)
			while ((s < e) && (*s++ != '\n'));

		find_words(s, e, rn);
	} while (1);

#ifdef FILE_READER_TIMES
	clock_gettime(CLOCK_MONOTONIC, t2);
	print_time_taken("Find Words", t1, t2);
#endif
	atomic_fetch_add(&readers_done, 1);
} // file_reader

//#define HASH_TABLE_TIMES

uint64_t
process_words()
{
	uint64_t spins = 0;

#ifdef HASH_TABLE_TIMES
	struct timespec t1[1], t2[1];
	clock_gettime(CLOCK_MONOTONIC, t1);
#endif

	// We do hash_init() and frq_init() here after the reader threads
	// start. This speeds up application load time as the OS needs to
	// clear less memory on startup.  Also, by initialising here, we
	// avoid blocking other work while initialisation occurs.

	// Build hash table and final key set
	hash_init();
	for (uint32_t *k = keys, key, pos = 0; ;) {
		if (pos >= num_words) {
			if (readers_done < num_readers) {
				spins++;
				asm("nop");
				continue;
			}
			if (pos >= num_words) {
				nkeys = k - keys;
				*k = 0;
				break;
			}
		}

		while ((key = wordkeys[pos]) == 0) {
			spins++;
			asm("nop");
		}

		*k = key;
		k += hash_insert(key, pos++);
	}

#ifdef HASH_TABLE_TIMES
	clock_gettime(CLOCK_MONOTONIC, t2);
	print_time_taken("Hash Insert", t1, t2);
#endif

	// All readers are done.  Collate character frequency stats
	frq_init();
	for (int rn = 0; rn < num_readers; rn++)
		for (int c = 0; c < 26; c++)
			frq[c].f += cfs[rn][c];

	return spins;
} // process_words

void
start_solvers()
{
	go_solve = 1;
} // start_solvers


// We create a worker pool like this because on virtual systems, especially
// on WSL, thread-creation is very expensive, so we only want to do it once
void *
work_pool(void *arg)
{
	struct worker *work = (struct worker *)arg;
	int worker_num = work - workers;

	if (pthread_detach(pthread_self()))
		perror("pthread_detach");

	// Wait until told to start
	while (!workers_start)
		asm("nop");

	if (worker_num < num_readers)
		file_reader(work);

#ifndef NO_FREQ_SETUP
	while (1) {
		int set_num = atomic_fetch_add(&setup_set, 1);

		if (set_num >= 26)
			break;

		set_tier_offsets(frq + set_num);
	}
#endif

	// Not gonna lie.  This is ugly.  We're busy-waiting until we get
	// told to start solving.  It shouldn't be for too long though...
	// I tried many different methods but this was always the fastest
	while (!go_solve)
		asm("nop");

	solve_work();
	return NULL;
} // work_pool

void
spawn_readers(char *start, size_t len)
{
	char *end = start + len;

	num_readers = (len / READ_CHUNK) + 1;

	if (num_readers > MAX_READERS)
		num_readers = MAX_READERS;
	if (num_readers > nthreads)
		num_readers = nthreads;
	if (num_readers < 1)
		num_readers = 1;

	for (int i = 0; i < num_readers; i++) {
		workers[i].start = start;
		workers[i].end = end;
	}

	// Need to zero out the word table so that the main thread can
	// detect when a word key has been written by a reader thread
	if (num_readers > 1)
		memset(wordkeys, 0, sizeof(wordkeys));

	// Start any waiting workers
	workers_start = 1;

	// Check if main thread must do reading
	if (num_readers < 2)
		file_reader(workers);
	else
		atomic_fetch_add(&readers_done, 1);

	// The main thread processes the words as the reader threads find them
	process_words();
} // spawn_readers

// File Reader.  We use mmap() for efficiency for both reading and processing
void
read_words(char *path)
{
	int fd;


	if ((fd = open(path, O_RDONLY)) < 0) {
		perror("open");
		exit(EXIT_FAILURE);
	}

	struct stat statbuf[1];
	if (fstat(fd, statbuf) < 0) {
		perror("fstat");
		exit(EXIT_FAILURE);
	}

	size_t len = statbuf->st_size;
	char *addr = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0);
	if (addr == MAP_FAILED) {
		perror("mmap");
		exit(EXIT_FAILURE);
	}

	// Safe to close file now.  mapping remains until munmap() is called
	close(fd);

	// Start file reader threads
	spawn_readers(addr, len);

	// We don't explicitly call munmap() as this can be slowish on some systems
	// Instead we'll just let the process terminate and it'll get unmapped then
} // read_words


// ********************* RESULTS WRITER ********************

// Solutions exists as a single character array assembled
// by the solver threads We just need to write it out.
void
emit_solutions()
{
	ssize_t len = num_sol << 5, written = 0;

	int solution_fd;
	if ((solution_fd = open(solution_filename, O_WRONLY | O_CREAT, 0644)) < 0) {
		fprintf(stderr, "Unable to open %s for writing\n", solution_filename);
		return;
	}

	// Truncate the file size if needed
	struct stat statbuf[1];
	if (fstat(solution_fd, statbuf) < 0) {
		perror("fstat");
		exit(EXIT_FAILURE);
	}

	if ((statbuf->st_size > len) && (ftruncate(solution_fd, len) < 0)) {
		perror("ftruncate");
		fprintf(stderr, "WARNING: Unable to truncate %s to %ld bytes\n",
			solution_filename, len);
	}

	// We loop here to handle any short writes that might occur
	while (written < len) {
		ssize_t ret = write(solution_fd, solutions + written, len - written);
		if (ret < 0) {
			fprintf(stderr, "Error writing to %s\n", solution_filename);
			perror("write");
			break;
		}
		written += ret;
	}

	close(solution_fd);
} // emit_solutions


// ********************* FREQUENCY HANDLER ********************

int
by_frequency_hi(const void *a, const void *b)
{
	return ((struct frequency *)b)->f - ((struct frequency *)a)->f;
} // by_frequency_hi



#ifdef __BMI2__

#define CALCULATE_SET_AND_END				\
	{						\
		struct tier *t = f->tiers +		\
			_pext_u32(mask, f->tmm);	\
		set = t->set;				\
		end = t->end;				\
	}

#else

#define CALCULATE_SET_AND_END				\
	{						\
		uint32_t tnum = !!(mask & f->tm1);	\
		tnum |= (!!(mask & f->tm2) << 1);	\
		tnum |= (!!(mask & f->tm3) << 2);	\
		tnum |= (!!(mask & f->tm4) << 3);	\
		tnum |= (!!(mask & f->tm5) << 4);	\
		tnum |= (!!(mask & f->tm6) << 5);	\
		struct tier *t = f->tiers + tnum;	\
		set = t->set;				\
		end = t->end;				\
	}

#endif


void
setup_tkeys(struct frequency *f, uint32_t t0_toff1,
		uint32_t t0_toff2, uint32_t t0_toff3)
{
	struct tier	*t0 = f->tiers;
	uint32_t	t0_len = t0->end - t0->set;
	uint32_t	*kp = t0->end + NUM_POISON;
	uint32_t	tm1 = f->tm1, tm2 = f->tm2;
	uint32_t	tm3 = f->tm3, tm4 = f->tm4;
	uint32_t	masks[16];

	// Define the mask bitmaps for splitting the sets
	masks[0]  = 0;
	masks[1]  = tm1;
	masks[2]  = tm2;
	masks[3]  = tm2 | tm1;
	masks[4]  = tm3;
	masks[5]  = tm3 | tm1;
	masks[6]  = tm3 | tm2;
	masks[7]  = tm3 | tm2 | tm1;

	masks[8]  = tm4;
	masks[9]  = tm4 | tm1;
	masks[10] = tm4 | tm2;
	masks[11] = tm4 | tm2 | tm1;
	masks[12] = tm4 | tm3;
	masks[13] = tm4 | tm3 | tm1;
	masks[14] = tm4 | tm3 | tm2;
	masks[15] = tm4 | tm3 | tm2 | tm1;

	// Create key arrays for each tier set mask
	for (uint32_t i = 0; i < 16; i++) {
		uint32_t toff1, toff2, toff3;
		uint32_t *ts, *te;
		uint32_t mask = masks[i];

		if (i == 0) {
			ts = t0->set;
			toff1 = t0_toff1;
			toff2 = t0_toff2;
			toff3 = t0_toff3;
			te = t0->end;
		} else {
			uint32_t *ks = t0->set;

			ts = kp;

			for (uint32_t len = t0_toff1; len--; )
				kp += !((*kp = *ks++) & mask);
			toff1 = kp - ts;

			for (uint32_t len = t0_toff2 - t0_toff1; len--; )
				kp += !((*kp = *ks++) & mask);
			toff2 = kp - ts;

			for (uint32_t len = t0_toff3 - t0_toff2; len--; )
				kp += !((*kp = *ks++) & mask);
			toff3 = kp - ts;

			for (uint32_t len = t0_len - t0_toff3; len--; )
				kp += !((*kp = *ks++) & mask);

			te = kp;

			for (uint32_t p = NUM_POISON; p--; )
				*kp++ = (uint32_t)(~0);
		}

		for (uint32_t tm5 = 0; tm5 < 2; tm5++) {
			for (uint32_t tm6 = 0; tm6 < 2; tm6++) {
				uint32_t tnum = i + ((tm5 << 4) | (tm6 << 5));
				struct tier *tts = f->tiers + tnum;

				if (tm5) {
					if (tm6) {
						tts->set = ts + toff2;
						tts->end = ts + toff3;
					} else {
						tts->set = ts + toff2;
						tts->end = te;
					}
				} else {
					if (tm6) {
						tts->set = ts + toff1;
						tts->end = ts + toff3;
					} else {
						tts->set = ts;
						tts->end = te;
					}
				}
			}
		}
	}
} // setup_tkeys

// This function looks like it's doing a lot, but because of good spatial
// and temportal localities each call typically takes ~1us on words_alpha
static void
set_tier_offsets(struct frequency *f)
{
	uint32_t key, mask, len;
	uint32_t *ks, *kp;

	// Wait here until all data is ready
	while (!f->ready)
		asm("nop");

	// "poison" NUM_POISON ending values with all bits set
	struct tier *t = f->tiers;
	ks = t->end;
	for (int p = NUM_POISON; p--; )
		*ks++ = (uint32_t)(~0);

	// Skip first set.  Nothing uses its subsets
	if (f == frq)
		goto set_tier_offsets_done;

	// "aeious" are the best static defaults
	f->tm1 = 1 << ('a' - 'a');
	f->tm2 = 1 << ('e' - 'a');
	f->tm3 = 1 << ('i' - 'a');
	f->tm4 = 1 << ('o' - 'a');
	f->tm5 = 1 << ('s' - 'a');
	f->tm6 = 1 << ('u' - 'a');

	// f->tmm is the aggregate of all individual tier masks
	f->tmm = f->tm1 | f->tm2 | f->tm3 | f->tm4 | f->tm5 | f->tm6;

	// Normalise the f->tm* masks according to bit order
	mask = f->tmm;

	f->tm1 = 1 << __builtin_ctz(mask); mask &= mask - 1;
	f->tm2 = 1 << __builtin_ctz(mask); mask &= mask - 1;
	f->tm3 = 1 << __builtin_ctz(mask); mask &= mask - 1;
	f->tm4 = 1 << __builtin_ctz(mask); mask &= mask - 1;
	f->tm5 = 1 << __builtin_ctz(mask); mask &= mask - 1;
	f->tm6 = 1 << __builtin_ctz(mask);

	// Organise full set into 2 subsets, that which
	// has tm5 followed by that which does not

	uint32_t toff1, toff2, toff3;
	mask = f->tm5;

	// First subset has tm5, and then not
	ks = kp = t->set;
	len = t->end - kp;
	for (; len--; ++ks)
		if ((key = *ks) & mask) {
			*ks = *kp;
			*kp++ = key;
		}
	toff2 = kp - t->set;

	// Now organise the first tm5 subset into that which
	// has tm6 followed by that which does not, and then
	// the second tm5 subset into that which does not
	// have tm6 followed by that which does

	mask = f->tm6;

	// First tm5 subset has tm6 then not
	ks = kp = t->set;
	len = toff2;
	for (; len--; ++ks)
		if ((key = *ks) & mask) {
			*ks = *kp;
			*kp++ = key;
		}
	toff1 = kp - t->set;

	// Second tm5 subset does not have tm6 then has
	ks = kp = t->set + toff2;
	len = t->end - kp;
	for (; len--; ++ks)
		if (!((key = *ks) & mask)) {
			*ks = *kp;
			*kp++ = key;
		}
	toff3 = kp - t->set;

	setup_tkeys(f, toff1, toff2, toff3);

set_tier_offsets_done:
	// Mark as done
	atomic_fetch_add(&setups_done, 1);
} // set_tier_offsets

// Specialised frequency sort, since we only need to swap the first 8 bytes
// of frequency sets at this point in time and each frequency set structure
// can be many hundreds of bytes, which wastes time if qsort is used
void
fsort()
{
	for (int i = 1; i < 26; ++i)
		for (int j = i; j; --j) {
			if (frq[j].f == 0)
				break;
			if (frq[j - 1].f && (frq[j].f > frq[j - 1].f))
				break;
			// Swap first 8 bytes
			uint64_t tmp = *(uint64_t *)(frq + j);
			*(uint64_t *)(frq + j) = *(uint64_t *)(frq + (j - 1));
			*(uint64_t *)(frq + (j - 1)) = tmp;
		}

	// Set the bit indices and the unmap table
	for (int i = 0, one = 1; i < 26; i++) {
		frq[i].b = __builtin_ctz(frq[i].m);
		unmap[frq[i].b] = (one << i);
	}
} // fsort

// The role of this function is to re-arrange the key set according to all
// words containing the least frequently used letter, and then scanning the
// remainder and so on until all keys have been assigned to sets. It achieves
// this by swapping keys in the key set and padding for any AVX operations
//
// Despite looking CPU and memory intensive, this function utilises strong
// spatial and temporal locality principles, and so runs in ~42us in practise
// without factoring in the calls to set_tier_offsets()
void
setup_frequency_sets()
{
	fsort();

	// Setup for key spray
	uint32_t *bp[32] __attribute__((aligned(64)));
	for (uint32_t i = 0; i < 26; i++)
		bp[i] = tkeys[i];

	// Spray keys to buckets
	for (uint32_t *kp = keys, key; (key = *kp++); ) {
		uint32_t mk = unmap[__builtin_ctz(key)];
		uint32_t k = key & (key - 1);

		mk |= unmap[__builtin_ctz(k)]; k &= k - 1;
		mk |= unmap[__builtin_ctz(k)]; k &= k - 1;
		mk |= unmap[__builtin_ctz(k)]; k &= k - 1;
		mk |= unmap[__builtin_ctz(k)];

		*bp[__builtin_ctz(mk)]++ = key;
	}

	// Start worker threads
	for (int i = 0; i < 26; i++) {
		struct frequency *f = frq + i;

		f->tiers[0].set = tkeys[i];
		f->tiers[0].end = bp[i];

		// Instruct any waiting worker thread to start setup
		// but we have to do it ourselves if single threaded
		f->ready = 1;
		if (nthreads == 1)
			set_tier_offsets(f);
	}

	// Wait for all setups to complete
	while(setups_done < 26)
		asm("nop");
} // setup_frequency_sets

#ifndef DONT_INCLUDE_MAIN

// ********************* MAIN SETUP AND OUTPUT ********************

int
main(int argc, char *argv[])
{
	struct timespec t1[1], t2[1], t3[1], t4[1], t5[1];
	char file[256];
	pthread_t tid[1];

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

	for (int i = 1; i < nthreads; i++)
		pthread_create(tid, NULL, work_pool, workers + i);

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
		struct tier *t = frq[i].tiers;
		char c = 'a' + __builtin_ctz(frq[i].m);
		printf("%c set_length=%4ld\n", c, t->end - t->set);
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

#endif
