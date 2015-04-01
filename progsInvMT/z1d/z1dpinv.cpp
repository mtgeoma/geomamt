#include "T1de.h"
#include <cstdlib>
#include<iostream>

void z1dinv(double T, vector<double> rho, vector<double> h,
            complex<double> &yfit, vector<complex<double> > &dZdrho) {
  // Zp (impedancia calculada pela produtoria)
  
  if(rho.size()!=dZdrho.size()) {
      cerr<<"rho e dZdrho devem ter o mesmo tamanho"<<endl;
      exit(1);
  }
  // inicializa vetor da derivada da produtoria
  vector<T1d> delT1d(rho.size()-1);
  for(int i=0;i<rho.size()-1;i++)
    delT1d[i].reset(rho[i],h[i],T);
  delT1d[0].delrho();
  T1d Mp(rho[0],h[0],T);

  for(int i=1; i<rho.size()-1; i++) {
    T1d Mi(rho[i],h[i],T);
    Mp*=Mi;
    for(int j=0;j<rho.size()-1;j++)
      if(j!=i) delT1d[j]*=Mi;
    Mi.delrho();
    delT1d[i]*=Mi;
  }
  complex<double> Zp(1.,1.);
  Zp*=sqrt(M_PI*u0*rho[rho.size()-1]/T);
  complex<double> Znum=Zp*Mp.xx()+Mp.xy();
  complex<double> Zden=Zp*Mp.yx()+Mp.yy();
  yfit=Znum/Zden;
  for(int i=0;i<rho.size()-1;i++) {
    dZdrho[i]=Zden*(Zp*delT1d[i].xx()+delT1d[i].xy())-Znum*(Zp*delT1d[i].yx()+delT1d[i].yy());
    dZdrho[i]/=pow(Zden,2);
  }
  dZdrho[rho.size()-1]=(Zden*Mp.xx()-Znum*Mp.yx())/pow(Zden,2);
  dZdrho[rho.size()-1]*=Zp/(2.*rho[rho.size()-1]);
}
