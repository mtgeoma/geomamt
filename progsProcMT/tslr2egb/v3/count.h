#ifndef COUNT_H
#define COUNT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "calendar.h"
#include "bytes.h"
#include "readtbl.h"

using namespace std;

typedef struct {
    long jul; // julian day
    long min; // minutes after julian day
} jtime;

typedef struct {
    jtime swin;
    jtime ewin;
    jtime sdata;
    jtime edata;
    long scount;
    long ecount;
    long count;
    long srpm;
    unsigned char nchan;
    unsigned char tagsize;
} timecount;

timecount startcount(string instr, table &tbl);

#endif // COUNT_H
