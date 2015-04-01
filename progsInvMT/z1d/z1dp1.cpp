#include "fwd1d.h"
#include "T1d.h"

complex<double> z1d(double T, vector<camada> C) {
  // Zp (impedancia calculada pela produtoria)

  T1d Mp(C[0].rho,C[0].h,T);

  for(int i=1; i<C.size()-1; i++) {
    T1d Mi(C[i].rho,C[i].h,T);
    Mp*=Mi;
  }
  complex<double> Zp(1.,1.);
  Zp*=sqrt(M_PI*u0*C[C.size()-1].rho/T);
  return (Zp*Mp.xx()+Mp.xy())/(Zp*Mp.yx()+Mp.yy());
}
