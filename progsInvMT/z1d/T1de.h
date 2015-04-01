#ifndef T1DE_H
#define T1DE_H

#include<cmath>
#include<complex>
#include<vector>

using namespace std;

#ifndef U0_CONST
#define U0_CONST
const double u0=4.e-7*M_PI;
#endif

void z1dinv(double T, vector<double> rho, vector<double> h,
            complex<double> &yfit, vector<complex<double> > &dZdrho);

class T1d {
  double rho, h;
  double z, a;
  complex<double> t11, t12, t21, t22;
public:
  T1d(double rho1=1., double h1=0., double T=1.) {
    reset(rho1,h1,T);
  }
  void reset(double rho1, double h1, double T) {
    rho=rho1;
    h=h1;
    z=sqrt(2.*M_PI*u0*rho/T);
    a=-sqrt(2.)*z*h/rho;
    t11=complex<double>(1.+exp(a)*cos(a),exp(a)*sin(a));
    t12=complex<double>( 1.-exp(a)*(cos(a)-sin(a)),1.-exp(a)*(cos(a)+sin(a)))*(z/sqrt(2.));
    t21=complex<double>(1.-exp(a)*(cos(a)+sin(a)),-1.+exp(a)*(cos(a)-sin(a)))/(z*sqrt(2.));
    t22=t11;
  }
  complex<double> xx(){return t11;}
  complex<double> xy(){return t12;}
  complex<double> yx(){return t21;}
  complex<double> yy(){return t22;}
  T1d& operator*=(T1d a) {
    complex<double> m[4]={t11,t12,t21,t22};
    t11=m[0]*a.xx()+m[1]*a.yx();
    t12=m[0]*a.xy()+m[1]*a.yy();
    t21=m[2]*a.xx()+m[3]*a.yx();
    t22=m[2]*a.xy()+m[3]*a.yy();
    return *this;
  }
  void delrho() { // derivada de T1d
    t11= (complex<double>(z,z)/sqrt(2.))*h*(complex<double>(1.,1.)*exp(a))/pow(rho,2);
    t12=-(complex<double>(z,z)/sqrt(2.))*t11+t12/(2.*rho);
    t21=-t11/(complex<double>(z,z)/sqrt(2.))-t21/(2.*rho);
    t22=t11;
  }
};

#endif // T1DE_H
