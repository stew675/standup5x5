#define	HASHSZ          30383		// Emperically derived optimum
#define MAX_READERS        15    	// Virtual systems don't like too many readers
#define	MAX_SOLUTIONS    8192
#define	MAX_WORDS        8192
#define	MAX_THREADS        64

static const char	*solution_filename = "solutions.txt";

// Worker thread state
static struct worker {
	char     *start;
	char     *end;
} workers[MAX_THREADS] __attribute__ ((aligned(64)));

// Set Pointers
struct tier {
	uint32_t	*s;	// Pointer to set
	uint32_t	l;	// Length of set
	uint32_t	toff1;	// Tiered Offset 1
	uint32_t	toff2;	// Tiered Offset 2
	uint32_t	toff3;	// Tiered Offset 3
};

// Character frequency recording
static struct frequency {
	struct tier	sets[16];	// Tier Sets
	uint32_t	m;		// Mask (1 << (c - 'a'))
	uint32_t	tm1;		// Tiered Mask 1
	uint32_t	tm2;		// Tiered Mask 2
	uint32_t	tm3;		// Tiered Mask 3
	uint32_t	tm4;		// Tiered Mask 4
	uint32_t	tm5;		// Tiered Mask 5
	uint32_t	tm6;		// Tiered Mask 6

	int32_t		f;		// Frequency
	atomic_int	pos;		// Position within a set
	uint32_t	pad[7];		// Pad to 448 bytes (7 * 64)
} frq[26] __attribute__ ((aligned(64)));

// Keep frequently modified atomic variables on their own CPU cache line
atomic_int 	num_words	__attribute__ ((aligned(64))) = 0;
atomic_int	file_pos	__attribute__ ((aligned(64))) = 0;
atomic_int	num_sol		__attribute__ ((aligned(64))) = 0;
atomic_int	readers_done = 0;
atomic_int	solvers_done = 0;

static int32_t	min_search_depth __attribute__ ((aligned(64))) = 0;
static int	write_metrics = 0;
static int	nthreads = 0;
static int	nkeys = 0;
static uint32_t hash_collisions = 0;

// We build the solutions directly as a character array to write out when done
static char     solutions[MAX_SOLUTIONS * 30] __attribute__ ((aligned(64)));

// Key Hash Entries
// We keep keys and positions in separate array because faster to initialise
uint32_t keymap[HASHSZ] __attribute__ ((aligned(64)));
uint32_t posmap[HASHSZ] __attribute__ ((aligned(64)));

// Allow for up to 3x the number of unique non-anagram words
static char     words[MAX_WORDS * 15] __attribute__ ((aligned(64)));
static uint32_t wordkeys[MAX_WORDS * 3] __attribute__ ((aligned(64)));

// We add 1024 here to MAX_WORDS to give us extra space to perform vector
// alignments for the AVX functions.  At the very least the keys array must
// be 32-byte aligned, but we align it to 64 bytes anyway
static uint32_t	keys[MAX_WORDS + 1024] __attribute__ ((aligned(64)));
static uint32_t	tkeys[MAX_WORDS * 10] __attribute__ ((aligned(64)));

// Here we pad the frequency counters to 32, instead of 26.  With the 64-byte
// alignment, this ensures all counters for each reader exist fully in 2 cache
// lines independent to each reader, thereby minimising CPU cache contention
static uint32_t cf[MAX_READERS][32] __attribute__ ((aligned(64))) = {0};

void
print_time_taken(char *label, struct timespec *ts, struct timespec *te)
{
	int64_t time_taken = 1000000000LL;	// Number of ns in 1s
	time_taken *= (te->tv_sec - ts->tv_sec);
	time_taken += (te->tv_nsec - ts->tv_nsec);

	printf("%-20s = %ld.%06lus\n", label, time_taken / 1000000000, (time_taken % 1000000000) / 1000);
} // print_time_taken
 
//********************* INIT FUNCTIONS **********************

static void
frq_init()
{
	memset(frq, 0, sizeof(frq));

	for (int b = 0; b < 26; b++)
		frq[b].m = (1UL << b);	// The bit mask
} // frq_init

static void
hash_init()
{
	memset(keymap, 0, sizeof(keymap));
} // hash_init


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
	uint32_t one = 1, a = 'a';
	uint32_t key;

	key  = (one << (*wd++ - a));
	key |= (one << (*wd++ - a));
	key |= (one << (*wd++ - a));
	key |= (one << (*wd++ - a));
	key |= (one << (*wd   - a));

	return key;
} // calc_key

//********************* HASH TABLE FUNCTIONS **********************

// A very simple for-purpose hash map implementation.  Used to
// lookup words given the key representation of that word
#define key_hash(x)	(x % HASHSZ)
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

		// Handle full hash table condition
		if (++col == HASHSZ)
			return 0;

		if (++hashpos == HASHSZ)
			hashpos -= HASHSZ;
	} while (1);

	// Now insert at hash location
	keymap[hashpos] = key;
	posmap[hashpos] = pos * 5;

	hash_collisions += col;

	return key;
} // hash_insert

const char *
hash_lookup(uint32_t key, const char *wp)
{
	uint32_t col = 0, hashpos = key_hash(key);

	do {
		// Check for a match
		if (keymap[hashpos] == key)
			break;

		// Check the not-in-hash scenario
		if (keymap[hashpos] == 0)
			return NULL;

		// Handle full hash table condition
		if (++col == HASHSZ)
			return NULL;

		if (++hashpos == HASHSZ)
			hashpos -= HASHSZ;
	} while (1);

	hash_collisions += col;

	return wp + posmap[hashpos];
} // hash_lookup
#undef key_hash


// ********************* FILE READER ********************

void
process_five_word(char *w, uint32_t *ft)
{
	uint32_t key = calc_key(w);
	if (__builtin_popcount(key) == 5) {
		int pos = atomic_fetch_add(&num_words, 1);

		// Get the key into the list as soon as possible
		// to prevent holding up the hash table builder
		wordkeys[pos] = key;

		// Update frequency table and
		// copy word at the same time
		char a = 'a';
		char *to = words + (5 * pos);
		ft[(*to++ = *w++) - a]++;
		ft[(*to++ = *w++) - a]++;
		ft[(*to++ = *w++) - a]++;
		ft[(*to++ = *w++) - a]++;
		ft[(*to   = *w  ) - a]++;
	}
} // process_five_word

// Because of setup overheads, AVX2 Scanning
// benefits from larger read chunk sizes

#ifdef USE_AVX2_SCAN
#define READ_CHUNK        32768		// Appears to be optimum
#else
#define READ_CHUNK        10240		// Appears to be optimum
#endif

void
find_words(char *s, char *e, uint32_t rn)
{
	char a = 'a', z = 'z', nl = '\n';
	uint32_t *ft = cf[rn];

#ifdef USE_AVX2_SCAN
	// The following code makes the assumption that all
	// words will immediately follow a newline character

#ifdef __AVX512F__
	// AVX512 is about 15-20% faster than AVX2 for processing the words

	// Prepare 3 constant vectors with newlines, a's and z's
	__m512i nvec = _mm512_set1_epi8(nl);
	__m512i avec = _mm512_set1_epi8(a);
	__m512i zvec = _mm512_set1_epi8(z);

	e -= 64;
	while (s < e) {
		int pos = 0;

		// Unaligned load of a vector with the next 32 characters
		__m512i wvec = _mm512_loadu_si512((const __m512i_u *)s);

		// Find the newlines in the word vector
		// nmask will have a 1-bit for every newline in the vector
		uint64_t nmask = _mm512_cmpeq_epi8_mask(nvec, wvec);

		// Find the lower-case letters in the word vector
		// wmask will have a 0-bit for every lower-case letter in the vector
		uint64_t wmask = _mm512_cmp_epi8_mask(wvec, avec, _MM_CMPINT_LT) |
					  _mm512_cmp_epi8_mask(zvec, wvec, _MM_CMPINT_LT);

		// Get number of newlines in the vector  NB: __builtin_popcount()
		// against a 64-bit value is SLOW, so we do it this way for speed
		int nls = __builtin_popcount((uint32_t)(nmask >> 32)) +
				   __builtin_popcount((uint32_t)(nmask & 0xffffffff));

		// Handle words over 32 characters in length
		if (nls == 0) {
			for (s += 32; s < e && *s++ != nl; );
			continue;
		}

		// Process all complete words in the vector
		while (nls--) {
			// Process word if it has exactly 5 letters
			if (__builtin_ctz((uint32_t)(wmask >> pos)) == 5)
				process_five_word(s + pos, ft);

			// Get position of next word
			pos += (__builtin_ctz((uint32_t)(nmask >> pos)) + 1);
		}
		s += pos;
	}
	e += 64;
#else
	// Prepare 3 constant vectors with newlines, a's and z's
	__m256i nvec = _mm256_set1_epi8(nl);
	__m256i avec = _mm256_set1_epi8(a);
	__m256i zvec = _mm256_set1_epi8(z);

	e -= 32;
	while (s < e) {
		int pos = 0;

		// Unaligned load of a vector with the next 32 characters
		__m256i wvec = _mm256_loadu_si256((const __m256i_u *)s);

		// Find the newlines in the word vector
		// nmask will have a 1-bit for every newline in the vector
		__m256i nres = _mm256_cmpeq_epi8(nvec, wvec);
		uint32_t nmask = _mm256_movemask_epi8(nres);

		// Find the lower-case letters in the word vector
		// wmask will have a 0-bit for every lower-case letter in the vector
		__m256i wres = _mm256_or_si256(_mm256_cmpgt_epi8(avec, wvec),
					       _mm256_cmpgt_epi8(wvec, zvec));
		uint32_t wmask = _mm256_movemask_epi8(wres);

		// Get number of newlines in the vector
		int nls = __builtin_popcount(nmask);

		// Handle words over 32 characters in length
		if (nls == 0) {
			for (s += 32; s < e && *s++ != nl; );
			continue;
		}

		// Process all complete words in the vector
		while (nls--) {
			// Process word if it has exactly 5 letters
			if (__builtin_ctz(wmask >> pos) == 5)
				process_five_word(s + pos, ft);

			// Get position of next word
			pos += (__builtin_ctz(nmask >> pos) + 1);
		}
		s += pos;
	}
	e += 32;
#endif
#endif

	char c, *w;
	for (w = s; s < e; w = s) {
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;
		c = *s++; if ((c < a) || (c > z)) continue;

		// We've now found 5 [a..z] characters in a row
		c = *s++; if ((c < a) || (c > z)) 
			process_five_word(w, ft);

		// Just quickly find the next line
		while (c != nl)
			c = *s++;
	}
} // find_words

void
file_reader(struct worker *work)
{
	uint32_t rn = work - workers;

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

	atomic_fetch_add(&readers_done, 1);
} // file_reader

static int num_readers = 0;

uint64_t
process_words()
{
	uint64_t spins = 0;

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

		if (hash_insert(key, pos))
			*k++ = key;

		pos++;
	}

	// All readers are done.  Collate character frequency stats
	frq_init();
	for (int c = 0; c < 26; c++)
		for (int r = 0; r < num_readers; r++)
			frq[c].f += cf[r][c];

	return spins;
} // process_words

volatile int go_solve = 0;
static void solve_work();

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

	if (worker_num < num_readers)
		file_reader(work);

	if (worker_num == num_readers)
		for (int i = worker_num + 1; i < nthreads; i++) {
			pthread_t tid[1];
			pthread_create(tid, NULL, work_pool, workers + i);
		}
	
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
	int main_must_read = 1;
	char *end = start + len;
	pthread_t tid[1];

	num_readers = len / (READ_CHUNK << 3);

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

	if (num_readers > 1) {
		// Need to zero out the word table so that the main thread can
		// detect when a word key has been written by a reader thread
		// The main thread doesn't do reading if num_readers > 1
		memset(wordkeys, 0, sizeof(wordkeys));

		// Spawn reader threads
		for (int i = 1; i < num_readers; i++)
			pthread_create(tid, NULL, work_pool, workers + i);

		if (num_readers > 3)
			main_must_read = 0;
	}

	// Spawn a thread that will create the rest of the worker threads
	if (num_readers < nthreads)
		pthread_create(tid, NULL, work_pool, workers + num_readers);

	// Check if main thread must do reading
	if (main_must_read)
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
	ssize_t len = num_sol * 30, written = 0;

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

// Sort lowest to highest, but treat 0 values as "infinite"
int
by_frequency_lo(const void *a, const void *b)
{
	if (((struct frequency *)a)->f == ((struct frequency *)b)->f)
		return 0;
	if (((struct frequency *)a)->f == 0)
		return 1;
	if (((struct frequency *)b)->f == 0)
		return -1;
	return ((struct frequency *)a)->f - ((struct frequency *)b)->f;
} // by_frequency_lo

int
by_frequency_hi(const void *a, const void *b)
{
	return ((struct frequency *)b)->f - ((struct frequency *)a)->f;
} // by_frequency_hi

// This function looks like it's doing a lot, but because of good spatial
// and temportal localities each call typically takes ~1us on words_alpha
void
set_tier_offsets(struct frequency *f)
{
	struct tier *t = f->sets;

	f->tm1 = frq[25].m;
	f->tm2 = frq[24].m;
	f->tm3 = frq[23].m;
	f->tm4 = frq[22].m;
	f->tm5 = frq[21].m;
	f->tm6 = frq[20].m;

	uint32_t key, mask, len;
	uint32_t *ks, *kp;

	// Organise full set into 2 subsets, that which
	// has tm1 followed by that which does not

	mask = f->tm5;

	// First subset has tm1, and then now
	ks = kp = t->s;
	len = t->l;
	for (; len--; ks++) {
		key = *ks;
		if (key & mask) {
			*ks = *kp;
			*kp++ = key;
		}
	}
	t->toff2 = kp - t->s;

	// Now organise the first tm1 subset into that which
	// has tm2 followed by that which does not, and then
	// the second tm1 subset into that which does not
	// have tm2 followed by that which does

	mask = f->tm6;

	// First tm1 subset has tm2 then not
	ks = kp = t->s;
	len = t->toff2;
	for (; len--; ks++) {
		key = *ks;
		if (key & mask) {
			*ks = *kp;
			*kp++ = key;
		}
	}
	t->toff1 = kp - t->s;

	// Second tm1 subset does not have tm2 then has
	ks = kp = t->s + t->toff2;
	len = t->l - t->toff2;
	for (; len--; ks++) {
		key = *ks;
		if (!(key & mask)) {
			*ks = *kp;
			*kp++ = key;
		}
	}
	t->toff3 = kp - t->s;
} // set_tier_offsets

void
setup_tkeys(struct frequency *f, int num_poison)
{
	static uint32_t	ntkeys = 0;
	struct tier	*t0 = f->sets;
	uint32_t	tm1 = f->tm1, tm2 = f->tm2;
	uint32_t	tm3 = f->tm3, tm4 = f->tm4;
	uint32_t	*kp = tkeys + ntkeys, *ks, len, key;
	uint32_t	masks[16];

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

	for (uint32_t mask, i = 1; i < 16; i++) {
		struct tier *ts = f->sets + i;
		mask = masks[i];

		ks = t0->s;
		ts->s = kp;

		len = t0->toff1;
		while (len--)
			if (!((key = *ks++) & mask))
				*kp++ = key;
		ts->toff1 = kp - ts->s;

		len = t0->toff2 - t0->toff1;
		while (len--)
			if (!((key = *ks++) & mask))
				*kp++ = key;
		ts->toff2 = kp - ts->s;

		len = t0->toff3 - t0->toff2;
		while (len--)
			if (!((key = *ks++) & mask))
				*kp++ = key;
		ts->toff3 = kp - ts->s;

		len = t0->l - t0->toff3;
		while (len--)
			if (!((key = *ks++) & mask))
				*kp++ = key;
		ts->l = kp - ts->s;

		for (uint32_t p = num_poison; p--; )
			*kp++ = (uint32_t)(~0);
	}

	ntkeys = kp - tkeys;
//	printf("ntkeys = %u\n", ntkeys);
} // setup_tkeys


// The role of this function is to re-arrange the key set according to all
// words containing the least frequently used letter, and then scanning the
// remainder and so on until all keys have been assigned to sets. It achieves
// this by swapping keys in the key set and padding for any AVX operations
//
// Despite looking CPU and memory intensive, this function utilises strong
// spatial and temporal locality principles, and so runs in ~42us in practise
// without factoring in the calls to set_tier_offsets()
void
setup_frequency_sets(int num_poison)
{
	struct frequency *f = frq;
	struct tier *t;
	uint32_t *kp = keys;

	qsort(f, 26, sizeof(*f), by_frequency_lo);

	// Now set up our scan sets by lowest frequency to highest
	for (int i = 0; i < 26; i++, f++) {
		uint32_t mask = f->m, *ks = kp, key;
		t = f->sets;

		for (t->s = kp; (key = *ks); ks++)
			if (key & mask) {
				*ks = *kp;
				*kp++ = key;
			}

		// Calculate the set length
		t->l = kp - t->s;

		// Update the min_search_depth if needed
		if (t->l > 0)
			min_search_depth = i - 3;

		// We can't do aligned AVX loads with the tiered approach.
		// "poison" num_poison ending values with all bits set
		// and ensure key set is 0 terminated for next loop
		for (uint32_t p = num_poison; p--; ) {
			*ks++ = *kp;
			*kp++ = (uint32_t)(~0);
		}
		*ks = 0;

		// Calculate our tier offsets
		set_tier_offsets(f);
		setup_tkeys(f, num_poison);
	}
} // setup_frequency_sets
