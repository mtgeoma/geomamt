#include "fwd1d.h"
#include "T1de.h"

complex<double> z1d(double T, vector<camada> C) {
  // Zp (impedancia calculada pela produtoria)

  vector<double> rho(C.size());
  vector<double> h(C.size());
  for(int i=0; i<C.size(); i++) {
    rho[i]=C[i].rho;
    h[i]=C[i].h;
  }

  complex<double> yfit;
  vector<complex<double> > dZdrho(C.size());
  z1dinv(T, rho, h, yfit, dZdrho);

  return yfit;
}
