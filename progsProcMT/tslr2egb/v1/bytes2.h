#ifndef BYTES2_H
#define BYTES2_H

#include <cctype>
#include <iostream>
#include <iomanip>
//#include <byteswap.h>

#ifdef WORDS_BIGENDIAN
    const bool bswap=true;
#else
    const bool bswap=false;
#endif

typedef short int2;
typedef long int4;
typedef long long int8;
typedef double float8;

bool correct_types(void);

//#ifndef _BYTESWAP_H
int2 bswap_16(int2 datum);
int4 bswap_32(int4 datum);
int8 bswap_64(int8 datum);
//#endif

#endif BYTES2_H
