#include <limits.h>

#define IRAN_MAX INT_MAX

void init_ran(int iseed); // Initialize

int iran_sign();    // Integer IRAN_MIN...IRAN_MAX
int iran();         // Integer 0...IRAN_MAX
double dran();      // Double [0, 1)
double dran_sign(); // Double [-1, 1)
