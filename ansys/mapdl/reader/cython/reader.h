#include <stdio.h>

#ifdef __linux__
  #include <stdint.h>
#endif


int read_nblock(char*, int*, double*, int, int*, int, int64_t*);
int read_eblock(char*, int*, int*, int, int, int64_t*);
int read_nblock_from_nwrite(char*, int*, double*, int);
int write_array_ascii(const char*, const double*, int nvalues);
