#ifndef COUNT_H
#define COUNT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "readhead.h"
#include "calendar.h"
using namespace std;

typedef struct {
    long jul; // julian day
    long min; // minutes after julian day
    long sec; // seconds after julian min
} jtime;

typedef struct {
    jtime swin;
    jtime ewin;
    jtime sdata;
    jtime edata;
    unsigned long scount;
    unsigned long ecount;
    unsigned long count;
    unsigned long imark;
    long srpm;
    unsigned char nchan;
    unsigned long sizehdr;
} timecount;

#endif // COUNT_H
