#ifndef READTBL_H
#define READTBL_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include "bytes.h"
#include "declination.h"

using namespace std;

typedef struct {                   /* AMX Time/Date Structure             */
      unsigned char sec;           /* seconds  (0-59)                     */
      unsigned char min;           /* minutes  (0-59)                     */
      unsigned char hr;            /* hours    (0-23)                     */
      unsigned char dy;            /* day      (1-31)                     */
      unsigned char mn;            /* month    (1-12)                     */
      unsigned char yr;            /* year     (0-99)                     */
      unsigned char ow;            /* day of week (Mon=1 to Sun=7)        */
      unsigned char cen;           /* 0 if time/date is incorrect         */
                                   /* century if time/date is correct     */
} amxtds;

typedef struct {
    int4   snum;          // Serial Number of MTU unit
    int4   srpm;          // Number of samples per minute
    string site;          // Site ID of the measurement
    string srvy;          // Survey ID
    int4   awin;          // Current data acquisition window
    int4   egn;           // Gain for "E" channels. 3,12,48
    int4   hgn;           // Gain for "H" channels. 2, 6, 8
    float8 eazm;          // Ex sensor azimuth
    float8 exln;          // Ex dipole length, m
    float8 eyln;          // Ey dipole length, m
    float8 hazm;          // Hx sensor azimuth
    int4   bsex;          // MAG sensor base line setting
    int4   bsey;          // MAG sensor base line setting
    int4   bsez;          // MAG sensor base line setting
    int4   elev;          // Elevation altitude of site [m]
    float  decl;          // magnetic declination
    double latg;          // Latitude
    double lngg;          // longitude
    float8 hnom;          // Coil nominal gain, mV/nT
    short  fieldtype;     // =0 for box calibration, =1 for system calibration
    unsigned short nchan; // number of channels
} table;

table  readtbl(string tblstr);
bool   positioning(ifstream& tbl, const char* label);
double convlat(const char* LATG);
double convlon(const char* LNGG);

#endif // READTBL_H
