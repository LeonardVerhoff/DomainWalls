#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "ran.h"

// Adapted from Newman & Barkema, Monte Carlo Methods in Statistical Physics.

#define SIZE 1279
#define OFFSET 216

#define A 2416
#define C 374441
#define M 1771875
#define CONV1 2423.9674
#define CONV2 (1.0 / (UINT_MAX + 1.0))
#define CONV3 (1.0 / (INT_MAX + 1.0))

static unsigned int vec[SIZE];
static int p, pp;

/*** Initialize the vector. ***/
void init_ran(int seed)
{
  static int already_initialized = 0;
  if (already_initialized)
    return;

  int i;

  // No seed? Then take it from the clock to make it different each time.
  if (!seed)
    seed = time(NULL);
  else
    printf("Seed for random number generator: %d.\n", seed);

  for (i = 0; i < SIZE; i++)
  {
    seed = (A * seed + C) % M;
    vec[i] = CONV1 * seed;
  }
  p = 0;
  pp = OFFSET;

  already_initialized = 1;
}

/*** This does the real job ***/
unsigned int uran() // Integer 0 ... 2^32 - 1
{
  if (--pp < 0)
    pp = SIZE - 1;
  if (--p < 0)
    p = SIZE - 1;

  vec[p] += vec[pp];
  return vec[p];
}

double dran() // Double [0, 1)
{
  return CONV2 * uran();
}

double dran_sign() // Double [-1, 1)
{
  return CONV3 * (int)uran();
}

int iran()
{
  return uran() & INT_MAX; // Integer 0 ... 2^31 - 1
}

int iran_sign()
{
  return (int)uran(); // Integer -2^31 ... 2^31 - 1
}
