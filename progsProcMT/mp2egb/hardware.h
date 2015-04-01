#ifndef HARDWARE_H
#define HARDWARE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "readhead.h"
using namespace std;

struct calibration {
    double gain,
           lpass,
           hpass;
};

void readhardware(string calstr, table tbl, vector<calibration> &hardware);

#endif // HARDWARE_H
