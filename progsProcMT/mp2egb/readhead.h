#ifndef READHEAD_H
#define READHEAD_H

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <inttypes.h>
#include <endian.h>
#include "declination.h"
#include "count.h"

# if __BYTE_ORDER == __LITTLE_ENDIAN
const bool bswap=false;
# else
const bool bswap=true;
# endif

using namespace std;

typedef float  float4;

typedef struct {
    int16_t day,
         month,
         year,
         hour,
         min;
} hdtime;

typedef struct {
    hdtime stime;         // start time
    int16_t   snum;          // Serial Number of LRMT unit
    int32_t   srpm;          // Number of samples per minute
    string site;          // Site ID of the measurement
    string srvy;          // Survey ID
    int16_t   awin;          // Current data acquisition window
    int16_t   exaz;          // Ex sensor azimuth
    int16_t   eyaz;          // Ey sensor azimuth
    int16_t   exln;          // Ex dipole length, m
    int16_t   eyln;          // Ey dipole length, m
    int16_t   hxaz;          // Hx sensor azimuth
    int16_t   fs;            // A to D factor = 305 microvolts/bit (full scale)
    int16_t   itgain;        // hardware telluric gain
    double magfac;        // magnetometer factor.
    int16_t   gain[5];       // gains. 0:Hx 1:Hy 2:Hz 3:Ex 4:Ey
    int32_t   bsex;          // MAG sensor base line setting
    int32_t   bsey;          // MAG sensor base line setting
    int32_t   bsez;          // MAG sensor base line setting
    float  decl;          // magnetic declination
    double latg;          // Latitude
    double lngg;          // longitude
    unsigned short nchan; // number of channels
    unsigned long sizehdr;// size of header [bytes]
    unsigned long sizefile;// size of file [bytes]
} table;

table readhead(string instr);

#endif // READHEAD_H
