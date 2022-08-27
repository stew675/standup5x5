#define	HASHSZ            39009
#define READ_CHUNK        10240
#define MAX_READERS          14    // Virtual systems don't like too many readers

// Word Hash Entries
static struct wordhash {
	uint32_t	key;
	uint32_t	pos;
} hashmap[HASHSZ] __attribute__ ((aligned(64)));

// Allow for up to 3x the number of unique non-anagram words
static char	words[MAX_WORDS * 15] __attribute__ ((aligned(64)));
static uint32_t	wordkeys[MAX_WORDS * 3] __attribute__ ((aligned(64)));

// We add 1024 here to MAX_WORDS to give us extra space to perform vector
// alignments for the AVX functions.  At the very least the keys array must
// be 32-byte aligned, but we align it to a typical system page boundary
static uint32_t	keys[MAX_WORDS + 1024] __attribute__ ((aligned(64)));

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
	memset(hashmap, 0, sizeof(hashmap));
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
calc_key(register const char *wd)
{
	register uint32_t one = 1, a = 'a';
	register uint32_t key;

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
#define key_hash(x)	((((size_t)x) << 26) % HASHSZ)
uint32_t
hash_insert(register uint32_t key, register uint32_t pos)
{
	register struct wordhash *h = hashmap + key_hash(key);
	register struct wordhash *e = hashmap + HASHSZ;
	register uint32_t col = 0;

	do {
		// Check if duplicate key
		if (h->key == key)
			return 0;

		// Check if we can insert at this position
		if (h->key == 0)
			break;

		col++;
		h++;

		// Handle full hash table condition
		if (col == HASHSZ)
			assert(col < HASHSZ);

		if (h == e)
			h -= HASHSZ;
	} while (1);

	hash_collisions += col;

	// Now insert at hash location
	h->key = key;
	h->pos = pos;

	return key;
} // hash_insert

const char *
hash_lookup(register uint32_t key, register const char *wp)
{
	register struct wordhash *h = hashmap + key_hash(key);
	register struct wordhash *e = hashmap + HASHSZ;
	register uint32_t col = 0;
	do {
		// Check the not-in-hash scenario
		if (h->key == 0)
			return NULL;

		// Check if entry matches the key
		if (h->key == key)
			break;

		col++;
		h++;

		// Handle full hash table condition
		if (col == HASHSZ)
			assert(col < HASHSZ);

		if (h == e)
			h -= HASHSZ;
	} while (1);

	hash_collisions += col;
	return wp + (h->pos * 5);
} // hash_lookup
#undef key_hash


// ********************* FILE READER ********************
void
find_words(register char *s, register char *e)
{
	register char *w, *wp = words, a = 'a', z = 'z';

	while (s < e) {
		w = s;
		if (*s++ < a) continue;
		if (*s++ < a) continue;
		if (*s++ < a) continue;
		if (*s++ < a) continue;
		if (*s++ < a) continue;

		// Rescanning twice is fast because the < a above
		// eliminates almost everything not of interest,
		// and the word will be in the L1 cache already
		if (*s < a || *s > z) {
			s = w;
			if (*s++ > z) continue;
			if (*s++ > z) continue;
			if (*s++ > z) continue;
			if (*s++ > z) continue;
			if (*s++ > z) continue;

			// Reject words without 5 unique [a..z] characters
			register uint32_t key = calc_key(w);
			if (__builtin_popcount(key) == 5) {
				register int pos = atomic_fetch_add(&num_words, 1);
				register char *to = wp + (5 * pos);

				memcpy(to, w, 5);
				wordkeys[pos] = key;
			}
		}

		// Just quickly find the next line
		while (*s++ != '\n');
	}
} // find_words

void *
file_reader(void *arg)
{
	struct worker *work = (struct worker *)arg;

	if (work->tid)
		if (pthread_detach(work->tid))
			perror("pthread_detach");

	// The e = s + (READ_CHUNK + 1) below is done because each reader
	// (except the first) only starts at a newline.  If the reader
	// starts at the very start of a 5 letter word, that means that it
	// would skip that word.  By adding the extra 1 here, the reader
	// processing the chunk before it can catch the word that may
	// have been skipped by the reader ahead of it
	do {
		register char *s = work->start;
		s += atomic_fetch_add(&file_pos, READ_CHUNK);
		register char *e = s + (READ_CHUNK + 1);

		if (s > work->end)
			break;
		if (e > work->end)
			e = work->end;

		// Make sure to only start after a newline
		// if we are not at the start of the file
		if (s > work->start)
			while ((s < e) && (*s++ != '\n'));

		find_words(s, e);
	} while (1);

	atomic_fetch_add(&readers_done, 1);
	return NULL;
} // file_reader

void
process_words(int nr)
{
	register uint32_t spins = 0, *k = keys;
	register char *wp = words, *w;
	uint32_t freqs['z' + 1] = {0};

	// We do these initialisions here after the reader threads start
	// The speeds up application load time as the OS needs to clear
	// less memory on startup.  Also, by initialising here, we avoid
	// blocking other work while initialisation occurs.
	frq_init();
	hash_init();

	for (register uint32_t key, pos = 0; ;) {
		if (pos >= num_words) {
			if (readers_done < nr) {
				spins++;
				usleep(1);
				continue;
			}
			if (pos >= num_words)
				break;
		}

		if ((key = wordkeys[pos]) == 0)
			continue;

		if (hash_insert(key, pos) == 0) {
			pos++;
			continue;
		}

		// We have a non-anagram word. Update letter frequencies
		w = wp + (5 * pos++);
		freqs[(int)*w++]++; freqs[(int)*w++]++;
		freqs[(int)*w++]++; freqs[(int)*w++]++;
		freqs[(int)*w]++;
		*k++ = key;
	}
	nkeys = k - keys;
	*k = 0;
	for (int i = 0; i < 26; i++)
		frq[i].f = freqs[i + 'a'];
} // process_words

void
spawn_readers(char *start, size_t len)
{
	int nr = nthreads;
	char *end = start + len;

	if (nr > MAX_READERS)
		nr = MAX_READERS;

	for (uintptr_t i = 0; i < nr; i++) {
		workers[i].start = start;
		workers[i].end = end;
	}

	if (nr > 1) {
		// Need to zero out the word table so that the main thread can
		// detect when a full word has been written by a reader thread
		memset(wordkeys, 0, sizeof(wordkeys));

		for (int i = 1; i < nr; i++)
			pthread_create(&workers[i].tid, NULL, file_reader, workers + i);

		// The main thread doesn't do reading if nthreads > 1
		nr--;
	} else {
		// The main thread is the only thread so it has to read the file
		workers[0].tid = 0;
		file_reader(workers);
	}

	// The main thread processes the words as the reader threads find them
	process_words(nr);
} // spawn_readers

// File Reader.  Here we use mmap() for efficiency for both reading and processing
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

void
rescan_frequencies(int start, register uint32_t *set)
{
	register uint32_t key;
	struct frequency *map[26];

	for (int i = start; i < 26; i++)
		map[__builtin_ctz(frq[i].m)] = frq + i;

	// Reset the frequencies we haven't processed yet
	for (int i = start; i < 26; i++)
		frq[i].f = 0;

	while ((key = *set++)) {
		while (key) {
			int i = __builtin_ctz(key);
			map[i]->f++;
			key ^= ((uint32_t)1 << i);
		}
	}

	qsort(frq + start, 26 - start, sizeof(*frq), by_frequency_hi);
} // rescan_frequencies

// ********************* MAIN SETUP AND OUTPUT ********************

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
