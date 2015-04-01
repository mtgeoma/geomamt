#ifndef GDS_H
#define GDS_H

#include <complex>
#include <vector>
#include <string>
using namespace std;

typedef struct {
    double T;
    complex<double> x;
    complex<double> y;
    double xvar;
    double yvar;
} tipper;

vector<tipper> read_egbert (string instr);
vector<tipper> read_edi (string instr);
vector<tipper> read_jones (string instr);

#endif // GDS_H
