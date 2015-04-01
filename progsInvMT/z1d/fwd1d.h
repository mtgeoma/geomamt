#ifndef FWD1D_H
#define FWD1D_H

#include <cstdlib>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;

#ifndef U0_CONST
#define U0_CONST
const double u0=4.e-7*M_PI;
#endif

typedef struct {
  double rho;
  double h;
} camada;

complex<double> z1d(double T, vector<camada> C);

#endif // FWD_H
