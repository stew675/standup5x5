#include <sys/types.h>
#include <stdio.h>
#include <stdint.h>

uint32_t sets[26][4096];

uint32_t num_combos = 0;
int8_t	combos[14950][4];
int8_t cidx[26];

// Combination Grinder implementing classic nCr
//
// Notable Arguments and Variables
//
// N = Total number of entries (the n part of nCr)
// R = Subset length (the r part of nCr)
// D = Current Depth (ranges from 0..R)
// S = Start Point of this recursion
// E = End Point of this recursion

void
combo(int N, int R, int D, int S, uint32_t *sets, uint32_t mask, uint32_t *solution)
{
	if (D == R) {
		// commit solution
		for (int i = 0; i < R; i++)
			printf("%5d", cidx[i]);
		printf("\n");
		return;
	}

	uint32_t *set, key;
	int E = N - R + D + 1;
	for (int i = S; i < E; i++) {
		cidx[D] = sets[i];

		set = sets[i];
		while ((key = *set++)) {
			if (key & mask)
				continue;
			solution[D] = key;
			combo(N, R, D + 1, i + 1, sets, mask | key, solution);
		}
	}
} //combo


void
main()
{
	uint32_t solution[5];

	for (int i = 0; i < 15; i++)
		sets[i] = i;

	combo(6, 3, 0, 0, sets, 0, solution);
} // main
