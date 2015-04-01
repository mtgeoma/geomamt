#ifndef T1D_H
#define T1D_H

#include<cmath>
#include<complex>
using namespace std;

class T1d {
public:
  complex<double> m11, m12, m21, m22;
  T1d(double rho, double h, double T) {
    complex<double> z(1.,1.);
    z*=sqrt(M_PI*(4.e-7*M_PI)*rho/T);
    complex<double> exp0=exp(-2.*z*h/rho);
    m11=1.+exp0;
    m12=(1.-exp0)*z;
    m21=(1.-exp0)/z;
    m22=m11;
  }
  T1d& operator*=(T1d a) {
    complex<double> m[4]={m11,m12,m21,m22};
    m11=m[0]*a.m11+m[1]*a.m21;
    m12=m[0]*a.m12+m[1]*a.m22;
    m21=m[2]*a.m11+m[3]*a.m21;
    m22=m[2]*a.m12+m[3]*a.m22;
    return *this;
  }
};

#endif // T1D_H
