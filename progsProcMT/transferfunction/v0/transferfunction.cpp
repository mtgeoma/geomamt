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

#include "transferfunction.h"

inline double commutIm(complex<double> C1, complex<double> C2) {
	return C1.real()*C2.imag()-C2.real()*C1.imag();
}

inline double commutRe(complex<double> C1, complex<double> C2) {
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

inline void TransferFunction::set_prop_error(prop_error p_err=stodt) {
	p_error=p_err;
}
double TransferFunction::Period() {
	return T;
}
complex<double> TransferFunction::val(index r, index c) {
	return TF[r*2+c].val;
}
inline double TransferFunction::err(index r, index c) {
	return sqrt(TF[r*2+c].var);
}
inline double TransferFunction::coh(index r, index c) {
	return TF[r*2+c].coh;
}
double TransferFunction::weight(index r, index c) {
	return TF[r*2+c].w;
}
	
inline double TransferFunction::rho(index r, index c) {
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

inline double TransferFunction::phi(index r, index c) {
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
	return 1.;
}

double TransferFunction::mu() {
	return sqrt(abs(commutIm(D1,S2))+abs(commutIm(S1,D2)))/abs(D2);
}

double TransferFunction::muerr() {
	return 1.;
}

double TransferFunction::eta() {
	return sqrt(abs(commutIm(D1,S2)-commutIm(S1,D2)))/abs(D2);
}

double TransferFunction::etaerr() {
	return 1.;
}
double TransferFunction::Sigma() {
	return (norm(D1)+norm(S2))/norm(D2);
}

double TransferFunction::Sigmaerr() {
	return 1.;
}

double TransferFunction::strike1() {
	double n=2.*real(S2*conj(D1));
	double d=norm(D1)-norm(S2);
	return (180./M_PI)*atan(n/d)/4.;
}

double TransferFunction::strike1err() {
	return 1.;
}

double TransferFunction::strike2() {
	double B=commutIm(S1,S2)-commutIm(D1,D2);
	double A=commutIm(S1,D1)+commutIm(S2,D2);
	double strike=(180./M_PI)*atan(B/A)/2.;
	
	if(strike<-45.)
		return strike+90.;
	else if(strike>45.)
		return strike-90.;
	return strike;
}

double TransferFunction::strike2err() {
	return 1.;
}

// equation 30 in Bahr (1991)
// cos(delta)sen(delta) is canceled in the division
double TransferFunction::alfa(int signal) {
	double A1=commutIm(S1,D1)+commutIm(S2,D2);
	double A2=commutRe(S1,D1)+commutRe(S2,D2);
	
	double B1=commutIm(S1,S2)-commutIm(D1,D2);
	double B2=commutRe(S1,S2)-commutRe(D1,D2);
	
	double C1=commutIm(D1,S2)-commutIm(S1,D2);
	double C2=commutRe(D1,S2)-commutRe(S1,D2);
	
	double E2=commutRe(S1,S1)-commutRe(D2,D2);
	
	double sqrt_tmp=pow((B1*A2+A1*B2+C1*E2)/(2.*(A1*A2-C1*C2)),2);
	       sqrt_tmp-=(B1*B2-C1*C2)/(A1*A2-C1*C2);
	       sqrt_tmp=sqrt(sqrt_tmp);
	
	double alfae=(B1*A2+A1*B2+C1*E2)/(2.*(A1*A2-C1*C2));
	if(signal>0)
		alfae=(90./M_PI)*atan(alfae+sqrt_tmp);
	else
		alfae=(90./M_PI)*atan(alfae-sqrt_tmp);
	
	if(isnan(alfae))
		return 400.;
	return alfae;
}

// equation 13 in Bahr (1991)
double TransferFunction::beta(int signal) {
	if(signal==1)
		return (180./M_PI)*atan(-abs(TF[0].val)/abs(TF[2].val));
	else if(signal==2)
		return (180./M_PI)*atan(abs(TF[3].val)/abs(TF[1].val));
	else
		return 0.;
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
	else {
		return (180./M_PI)*atan2(Tzy,Tzx);
	}
}
