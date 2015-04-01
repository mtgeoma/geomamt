/***************************************************************************
                          transferfunction.h  -  description
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

#ifndef TRANSFERFUNCTION_H
#define TRANSFERFUNCTION_H


/**
  *@author Marcelo Banik de Pádua
  */

#include <cmath>
#include <complex>
#include <vector>

using namespace std;

struct TFunit {
	complex<double> val; // value
	double var;          // variance
	double coh;          // coherence
	double w;            // weight
	
	TFunit() {
	val=0.;
	var=0.;
	coh=1.;
	w=1.;
	};
};

inline double commutIm(complex<double> C1, complex<double> C2);
inline double commutRe(complex<double> C1, complex<double> C2);

// in field units [mV/(km*nT)]
// TF[0]=Zxx    TF[1]=Zxy
// TF[2]=Zyx    TF[3]=Zyy
// TF[4]=Tzx    TF[5]=Tzy
// index=row*2+column

// Zmod[0]=S1   Zmod[1]=S2
// Zmod[2]=D1   Zmod[3]=D2

class TransferFunction {
public:
	enum prop_error {stodt, gamble, chave, jones};
	enum index {x, y, z};
  	index r;
  	index c;
private:
	double T;
	vector<TFunit> TF;
	complex<double> S1, S2, D1, D2; // modified impedances
	double S1var, S2var, D1var, D2var;       // and its variance
	prop_error p_error; // type of propagation of errors
public:
	TransferFunction(double To, vector<TFunit> TFo);
	~TransferFunction();
	inline void set_prop_error(prop_error p_err=stodt);
	double TransferFunction::Period();
	complex<double> val(index r, index c);
	inline double err(index r, index c);
	inline double coh(index r, index c);
	double weight(index r, index c);

	inline double rho(index r, index c);
	double rhoerr(index r, index c);
	double rhomax(index r, index c);
	double rhomin(index r, index c);

	inline double phi(index r, index c);
	double phierr(index r, index c);
	double phimax(index r, index c);
	double phimin(index r, index c);

	double dD2_D2();

	double kappa(); // as defined in
	double mu();    // Barh (1990).
	double eta();   //
	double Sigma(); //

	double kappaerr(); // and its
	double muerr();    // errors
	double etaerr();   //
	double Sigmaerr(); //

	double strike1(); // tan(4*alfa)=2*Re(S2*D1)/(norm(D1)-norm(S2))
	double strike1err();

	double strike2(); // tan(2*alfa)=B/A
	double strike2err();

	double alfa(int signal);

	double beta(int signal);

	void rotate(double const teta);

	double ind_vector(int par);
	void MPS(double step, double& angle, double& phasesplit);
};

#endif
