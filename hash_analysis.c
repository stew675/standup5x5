#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "keys.h"

#define	HASHSZ            204800

uint32_t keymap[HASHSZ] __attribute__ ((aligned(64)));
uint32_t posmap[HASHSZ] __attribute__ ((aligned(64)));
uint32_t hash_collisions;


void
print_time_taken(char *label, struct timespec *ts, struct timespec *te)
{
	int64_t time_taken = 1000000000LL;	// Number of ns in 1s
	time_taken *= (te->tv_sec - ts->tv_sec);
	time_taken += (te->tv_nsec - ts->tv_nsec);

	printf("%-20s = %ld.%06lus\n", label, time_taken / 1000000000, (time_taken % 1000000000) / 1000);
} // print_time_taken
 

static void
hash_init(size_t len)
{
	memset(keymap, 0, len * sizeof(*keymap));
	hash_collisions = 0;
} // hash_init


//********************* HASH TABLE FUNCTIONS **********************

// A very simple for-purpose hash map implementation.  Used to
// lookup words given the key representation of that word
uint32_t
hash_insert(register uint32_t key, register uint32_t pos, register uint32_t hashsz)
{
	register uint32_t col = 0, hashpos = key % hashsz;

	do {
		if (keymap[hashpos] == 0)
			break;

		// Check if duplicate key
		if (keymap[hashpos] == key)
			return 0;

		// Handle full hash table condition
		if (++col == hashsz)
			return 0;

		if (++hashpos == hashsz)
			hashpos -= hashsz;
	} while (1);

	// Now insert at hash location
	keymap[hashpos] = key;
	posmap[hashpos] = pos * 5;

	hash_collisions += col;

	return key;
} // hash_insert

uint32_t
hash_lookup(register uint32_t key, register uint32_t hashsz)
{
	register uint32_t col = 0, hashpos = key % hashsz;

	do {
		if (keymap[hashpos] == key)
			break;

		if (keymap[hashpos] == 0)
			return UINT32_MAX;

		// Handle full hash table condition
		if (++col == hashsz)
			return UINT32_MAX;

		if (++hashpos == hashsz)
			hashpos -= hashsz;
	} while (1);

	hash_collisions += col;

	return posmap[hashpos];
} // hash_insert


int
main(int argc, char *argv[])
{
	struct timespec ts[1], te[1];
	int nkeys = sizeof(keys) / sizeof(*keys);
	int64_t min_time_taken = 1000000000LL;
	uint32_t min_col = 10000;

	register uint32_t hashsz = atoi(argv[1]);
	if (hashsz < (nkeys + 100))
		exit(1);
//	for (register uint32_t hashsz = nkeys + 100; hashsz < HASHSZ; hashsz++) {
		int64_t time_taken = 1000000000LL;
		uint32_t key, *k;

		clock_gettime(CLOCK_MONOTONIC, ts);

		hash_init(hashsz);

		k = keys;
		while ((key = *k++))
			hash_insert(key, (uint32_t)(k - keys), hashsz);

		k = keys;
		while ((key = *k++))
			hash_lookup(key, hashsz);

		clock_gettime(CLOCK_MONOTONIC, te);

		time_taken *= (te->tv_sec - ts->tv_sec);
		time_taken += (te->tv_nsec - ts->tv_nsec);

		if (time_taken < min_time_taken) {
			printf("HASHSZ = %-10u  Collisions = %-8u   %ld\n", hashsz, hash_collisions, time_taken);
			min_time_taken = time_taken;
		}

#if 0
		if (hash_collisions < min_col) {
			min_col = hash_collisions;
			printf("HASHSZ = %-10u  Collisions = %-8u   %ld\n", hashsz, hash_collisions, time_taken);
		}
#endif
//	}
	return 0;
} // main
