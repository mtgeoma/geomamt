/***************************************************************************
                          transferfunction.cpp  -  description
                             -------------------
    begin                : Fri Jan 12 2001
    copyright            : (C) 2001 by Marcelo Banik de Pádua
    email                : banik@dge.inpe.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <iostream>
#include "transferfunction.h"

using namespace std;

double commutIm(complex<double> C1, complex<double> C2) {
	return C1.real()*C2.imag()-C2.real()*C1.imag();
}

double commutRe(complex<double> C1, complex<double> C2) {
	return C1.real()*C2.real()+C2.imag()*C1.imag();
}

TransferFunction::TransferFunction(double To, vector<TFunit> TFo){
	if(TFo.size()!=6) {
		cerr<<"Transfer function at "<<To<<"s don't have 6 elements"<<endl;
		exit(1);
	}

	set_prop_error();
	T=To;
	TF=TFo;

	S1=TF[0].val+TF[3].val;
	S2=TF[1].val+TF[2].val;
	D1=TF[0].val-TF[3].val;
	D2=TF[1].val-TF[2].val;

	S1var=TF[0].var+TF[3].var;
	S2var=TF[1].var+TF[2].var;
	D1var=S1var;
	D2var=S2var;
}

TransferFunction::~TransferFunction(){
}

void TransferFunction::set_prop_error(prop_error p_err) {
	p_error=p_err;
}
double TransferFunction::Period() {
	return T;
}
complex<double> TransferFunction::val(index r, index c) {
	return TF[r*2+c].val;
}
double TransferFunction::err(index r, index c) {
	return sqrt(TF[r*2+c].var);
}
double TransferFunction::coh(index r, index c) {
	return TF[r*2+c].coh;
}
double TransferFunction::weight(index r, index c) {
	return TF[r*2+c].w;
}

double TransferFunction::rho(index r, index c) {
	return 0.2*T*norm(TF[r*2+c].val);
}

double TransferFunction::rhoerr(index r, index c) {
	switch (p_error) {
	case stodt:
		return sqrt(2.*TF[r*2+c].var/(pow(log(10.),2)*norm(TF[r*2+c].val)));
	case gamble: case chave:
		return sqrt(2.*pow(T,2)*norm(TF[r*2+c].val)*TF[r*2+c].var/25.);
	case jones:
		cerr<<"rhoerr() can't be used with Jones' propagation of error"<<endl;
		cerr<<"use rhomax() or rhomin()"<<endl;
		exit(1);
	default:
		cerr<<"unknown propagation of error"<<endl;
		exit(1);
	}
}

double TransferFunction::rhomax(index r, index c) {
	switch (p_error) {
	case stodt:
		return rho(r,c)*pow(10.,sqrt(2.*TF[r*2+c].var/(pow(log(10.),2)*norm(TF[r*2+c].val))));
	case gamble: case chave:
		return rho(r,c)+sqrt(2.*pow(T,2)*norm(TF[r*2+c].val)*TF[r*2+c].var/25.);
	case jones:
		return 0.2*T*pow(abs(TF[r*2+c].val)+sqrt(TF[r*2+c].var),2);
	default:
		cerr<<"unknown propagation of error"<<endl;
		exit(1);
	}
}

double TransferFunction::rhomin(index r, index c) {
	switch (p_error) {
	case stodt:
		return rho(r,c)*pow(10.,-1.*sqrt(2.*TF[r*2+c].var/(pow(log(10.),2)*norm(TF[r*2+c].val))));
	case gamble: case chave:
		return rho(r,c)-sqrt(2.*pow(T,2)*norm(TF[r*2+c].val)*TF[r*2+c].var/25.);
	case jones:
		return 0.2*T*pow(abs(TF[r*2+c].val)-sqrt(TF[r*2+c].var),2);
	default:
		cerr<<"unknown propagation of error"<<endl;
		exit(1);
	}
}

double TransferFunction::phi(index r, index c) {
	return (180./M_PI)*arg(TF[r*2+c].val);
}

double TransferFunction::phierr(index r, index c) {
	switch (p_error) {
	case gamble: case stodt:
		return (180./M_PI)*err(r,c)/(sqrt(2.)*abs(val(r,c)));
	case chave:
		return (180./M_PI)*asin(err(r,c)/abs(val(r,c)));
	case jones:
		return (180./M_PI)*atan(err(r,c)/abs(val(r,c)));
	default:
		cerr<<"unknown propagation of error"<<endl;
		exit(1);
	}
}

double TransferFunction::phimax(index r, index c) {
	switch (p_error) {
	case gamble: case stodt:
		return phi(r,c)+(180./M_PI)*err(r,c)/(sqrt(2.)*abs(val(r,c)));
	case chave:
		return phi(r,c)+(180./M_PI)*asin(err(r,c)/abs(val(r,c)));
	case jones:
		return phi(r,c)+(180./M_PI)*atan(err(r,c)/abs(val(r,c)));
	default:
		cerr<<"unknown propagation of error"<<endl;
		exit(1);
	}
}

double TransferFunction::phimin(index r, index c) {
	switch (p_error) {
	case gamble: case stodt:
		return phi(r,c)-(180./M_PI)*err(r,c)/(sqrt(2.)*abs(val(r,c)));
	case chave:
		return phi(r,c)-(180./M_PI)*asin(err(r,c)/abs(val(r,c)));
	case jones:
		return phi(r,c)-(180./M_PI)*atan(err(r,c)/abs(val(r,c)));
	default:
		cerr<<"unknown propagation of error"<<endl;
		exit(1);
	}
}

double TransferFunction::dD2_D2() {
	return sqrt(D2var)/abs(D2);
}

double TransferFunction::kappa() {
	return abs(S1)/abs(D2);
}

double TransferFunction::kappaerr() {
        return sqrt(S1var/norm(D2)+D2var*norm(S1)/(norm(D2)*norm(D2)));
}

double TransferFunction::mu() {
	return sqrt(abs(commutIm(D1,S2))+abs(commutIm(S1,D2)))/abs(D2);
}

double TransferFunction::muerr() {
        double C=abs(commutIm(D1,S2))+abs(commutIm(S1,D2));
        double Cvar=D1var*norm(S2)+S2var*norm(D1)+S1var*norm(D2)+D2var*norm(S1);
	return sqrt(Cvar/(4.*C*norm(D2))+D2var*C/(norm(D2)*norm(D2)));
}

double TransferFunction::eta() {
	return sqrt(abs(commutIm(D1,S2)-commutIm(S1,D2)))/abs(D2);
}

double TransferFunction::etaerr() {
        double C=abs(commutIm(D1,S2)-commutIm(S1,D2));
        double Cvar=D1var*norm(S2)+S2var*norm(D1)+S1var*norm(D2)+D2var*norm(S1);
	return sqrt(Cvar/(4.*C*norm(D2))+D2var*C/(norm(D2)*norm(D2)));
}
double TransferFunction::Sigma() {
	return (norm(D1)+norm(S2))/norm(D2);
}

double TransferFunction::Sigmaerr() {
  double Sigma=(norm(D1)+norm(S2))/norm(D2);
  return sqrt(4.*norm(D1)*D1var/(norm(D2)*norm(D2))+4.*norm(S2)*S2var/(norm(D2)*norm(D2))+4.*Sigma*Sigma*D2var/norm(D2));
}

double TransferFunction::strike1() {
	double n=2.*real(S2*conj(D1));
	double d=norm(D1)-norm(S2);
        double teta4=atan(n/d);
        // d_d2 2a. derivada dividido por 4 (o que não altera o sinal)
        double d_d2=-(d*cos(teta4)+n*sin(teta4));
        if(d_d2>0)
	    return (180./M_PI)*teta4/4.;
        else {
	    if((teta4/4.)<=0)
	        return (180./M_PI)*teta4/4.+45.;
            else
	        return (180./M_PI)*teta4/4.-45.;
        }
}

double TransferFunction::strike1err() {
	return 1.;
}

double TransferFunction::strike2() {
	double B=commutIm(S1,S2)-commutIm(D1,D2);
	double A=commutIm(S1,D1)+commutIm(S2,D2);
	double strike=(180./M_PI)*atan(B/A)/2.;

	if(strike<-45.)
		strike+=90.;
	else if(strike>45.)
		strike-=90.;
	return strike;
}

double TransferFunction::strike2err() {
	return 1.;
}

// equation 30 in Bahr (1991) corrected by Pracser & Szarka (1999)
double TransferFunction::alfa() {
	double a1=commutIm(S1,D1)+commutIm(S2,D2);
	double a2=commutRe(S1,D1)+commutRe(S2,D2);

	double b1=commutIm(S1,S2)-commutIm(D1,D2);
	double b2=commutRe(S1,S2)-commutRe(D1,D2);

	double c1=commutIm(D1,S2)-commutIm(S1,D2);
	double c2=commutRe(D1,S2)-commutRe(S1,D2);

	double e2=commutRe(S2,S2)-commutRe(D1,D1);
        double f2=2.*commutRe(D1,S2);


        double P1=0.5*(b1*a2+a1*b2+c1*e2)/(a1*a2-c1*c2+c1*f2);

        double P2=0.25*pow((b1*a2+a1*b2+c1*e2)/(a1*a2-c1*c2+c1*f2),2);

        double P3=(b1*b2-c1*c2)/(a1*a2-c1*c2+c1*f2);

        // nao ha raiz negativa
        if((P2-P3)<0)
		return 400.;

        double alfa1=atan(P1+sqrt(P2-P3))/2.;
        double alfa2=atan(P1-sqrt(P2-P3))/2.;

        double delta1=atan( c1/(a2*sin(2.*alfa1)-b2*cos(2.*alfa1)) );
        double delta2=atan( c1/(a2*sin(2.*alfa2)-b2*cos(2.*alfa2)) );

        double alfa;
	if(delta1<delta2)
		alfa=(180./M_PI)*alfa1;
	else
		alfa=(180./M_PI)*alfa2;


	if(alfa<-45.)
		alfa+=90.;
	else if(alfa>45.)
		alfa-=90.;

	return alfa;
}

// equation 13 in Bahr (1991)
double TransferFunction::beta(int signal) {
        vector <TFunit> TFtmp=TF;

	if(signal==1 || signal==2)
          this->rotate(this->strike2());
	else if(signal==3 || signal==4)
          this->rotate(this->alfa());

	if(signal%2==1) {
                double beta=(180./M_PI)*atan(-abs(TF[0].val)/abs(TF[2].val));
                TF=TFtmp;
		return beta;
        }
	else if(signal%2==0) {
                double beta=(180./M_PI)*atan(abs(TF[3].val)/abs(TF[1].val));
                TF=TFtmp;
		return beta;
        }
}

void TransferFunction::rotate(double const teta) {

    // rotate the impedance tensor
	double teta2=(M_PI/180.)*2.*teta;
	TF[0].val=( S1+S2*sin(teta2)+D1*cos(teta2))/2.;
	TF[1].val=( D2-D1*sin(teta2)+S2*cos(teta2))/2.;
	TF[2].val=(-D2-D1*sin(teta2)+S2*cos(teta2))/2.;
	TF[3].val=( S1-S2*sin(teta2)-D1*cos(teta2))/2.;

	TF[0].var=(S1var+S2var*pow(sin(teta2),2)+D1var*pow(cos(teta2),2))/4.;
	TF[1].var=(D2var+D1var*pow(sin(teta2),2)+S2var*pow(cos(teta2),2))/4.;
	TF[2].var=TF[1].var;
	TF[3].var=TF[0].var;

	// and the modified impedance
	complex<double> S2rot=S2*cos(teta2)-D1*sin(teta2);
	complex<double> D1rot=S2*sin(teta2)+D1*cos(teta2);
	double S2rotvar=S2var*pow(cos(teta2),2)+D1var*pow(sin(teta2),2);
	double D1rotvar=S2var*pow(sin(teta2),2)+D1var*pow(cos(teta2),2);

	S2=S2rot;
	D1=D1rot;
	S2var=S2rotvar;
	D1var=D1rotvar;
}

// TF[4]=Tzx    TF[5]=Tzy
double TransferFunction::ind_vector(int par) {
	double Tzx=TF[4].val.real();
	double Tzy=TF[5].val.real();

	if(par==0) { // amplitude
		return sqrt(pow(Tzx,2)+pow(Tzy,2));
	}
	else {  // invert azimut to point to condutor
                double azm=atan2(Tzy,Tzx);
                if(azm<0.)
                    azm+=M_PI;
                else if(azm>0.)
                    azm-=M_PI;

		return (180./M_PI)*azm;
	}
}

void TransferFunction::MPS(double step, double& angle, double& phasesplit) {
// find rotation angle that determine the greatest phase difference between
// the off-diagonal impedance tensor elements
// obs: Zyx shifted to 1st quadrant (+180 degrees)
	if(step<=0) {
            cerr<<"Maximum Phase Split step (defined by -A option)\n"
                <<"must be greater than zero"<<endl;
            exit(1);
        }

	angle=-90.;
        rotate(angle);
        phasesplit=(180./M_PI)*abs(arg(TF[1].val)-(arg(TF[2].val)+M_PI));
        for(double i=step;i<90.;i+=step) {
	    rotate(step);
            double split=(180./M_PI)*abs(arg(TF[1].val)-(arg(TF[2].val)+M_PI));
            if(split>phasesplit) {
                phasesplit=split;
                angle=i;
            }
        }
}
