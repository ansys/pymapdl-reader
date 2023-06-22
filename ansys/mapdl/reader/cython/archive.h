#include <stdio.h>
#include <stdlib.h>

#ifdef __linux__
  #include <stdint.h>
#endif


int write_nblock(FILE*, const int, const int, const int*, const double*,
                 const double*, int);
int write_nblock_float(FILE*, const int, const int, const int*, const float*,
                       const float*, int);
int write_eblock(FILE*, const int, const int*, const int*, const int*,
                 const int*, const int*, const uint8_t*, const int*,
                 const int*, const int*, const int*);
