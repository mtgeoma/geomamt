#include "fwd1d.h"

complex<double> z1d(double T, vector<camada> C) {
  // Zr (impedancia calculada recursivamente) é inicializada com a
  // impedancia intrinseca da última camada (semi-espaço infinito)

  complex<double> Zr(1.,1.);
  Zr*=sqrt(M_PI*u0*C[C.size()-1].rho/T);

  for(int i=C.size()-2; i>=0; i--) {
    complex<double> Zi(1.,1.);
    Zi*=sqrt(M_PI*u0*C[i].rho/T);
    complex<double> exp0=exp(-2.*Zi*C[i].h/C[i].rho);
    Zr=Zi*(Zr*(1.+exp0)+Zi*(1.-exp0))/(Zr*(1.-exp0)+Zi*(1.+exp0));
  }
  return Zr;
}
