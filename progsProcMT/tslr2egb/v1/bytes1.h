#ifndef BYTES1_H
#define BYTES1_H

#include <cctype>
#include <iostream>
#include <iomanip>
#include <byteswap.h>

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
//long bswap_32(long datum);
//short bswap_16(short datum);
//#endif

#endif BYTES1_H
