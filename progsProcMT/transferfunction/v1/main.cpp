/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Sex Jan 12 20:35:30 BRST 2001
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream.h>
#include <stdlib.h>
#include <iomanip>
#include<stdio.h>
#include "transferfunction.h"
#include "imputdata.h"
#include "commandline.h"

// argumentos: arquivo tipo extrator
int main(int argc, char *argv[])
{
	CommandLine entrada(argc, argv);
	DataType::ImputData j_data(entrada.file(),entrada.file_format());
	vector<TransferFunction> J;
	j_data.getdata(J);

	if(entrada.rot()!=0. && entrada.extract()!=CommandLine::phasesplit) {
		for(size_t i=0; i<J.size(); i++)
			J[i].rotate(entrada.rot());

		cout<<"Impedance tensor rotated by "<<entrada.rot()<<" degrees"<<endl;
	}

	switch (entrada.extract()) {
	case CommandLine::impedance:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			complex<double> Z=J[i].val(entrada.row(),entrada.col());
			double Zerr=J[i].err(entrada.row(),entrada.col());
			double weight=J[i].weight(entrada.row(),entrada.col());

			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Z.real()
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Z.imag()
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Zerr
				<<setw(entrada.width())<<setprecision(entrada.precision())<<weight
				<<endl;
		}
		break;
	case CommandLine::rhophi:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double rho=J[i].rho(entrada.row(),entrada.col());
			double rhomin=J[i].rhomin(entrada.row(),entrada.col());
			double rhomax=J[i].rhomax(entrada.row(),entrada.col());
			double phi=J[i].phi(entrada.row(),entrada.col());
			double phimin=J[i].phimin(entrada.row(),entrada.col());
			double phimax=J[i].phimax(entrada.row(),entrada.col());
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<rho
				<<setw(entrada.width())<<setprecision(entrada.precision())<<rhomin
				<<setw(entrada.width())<<setprecision(entrada.precision())<<rhomax
				<<setw(entrada.width())<<setprecision(entrada.precision())<<phi
				<<setw(entrada.width())<<setprecision(entrada.precision())<<phimin
				<<setw(entrada.width())<<setprecision(entrada.precision())<<phimax
				<<endl;
		}
		break;
	case CommandLine::mackie:
		bool te;
		if(entrada.row()==TransferFunction::x && entrada.col()==TransferFunction::y)
				te=true;
		else if(entrada.row()==TransferFunction::y && entrada.col()==TransferFunction::x)
				te=false;
		else {
				cerr<<"Mackie output data could be only the xy or yx component"<<endl;
				exit(1);
		}
		printf("%-14d%6.4f\n",J.size(),1.0);
		for(size_t i=J.size()-1; i>0; i--) {
			double T=J[i].Period();
			double rho=J[i].rho(entrada.row(),entrada.col());
			double phi;
         if (te)
					phi=-1.*J[i].phi(entrada.row(),entrada.col());
         else
					phi=-1.*(180.+J[i].phi(entrada.row(),entrada.col()));
			double Z=abs(J[i].val(entrada.row(),entrada.col()));
			double e=J[i].err(entrada.row(),entrada.col());
			double phierr=e/Z;
			printf("%-14.6E(%10.2f,%10.2f)",T,rho,phi);
			printf("%10.4f%10.4f\n",2.*phierr,phierr);
		}
		break;
	case CommandLine::phasesplit:
		for(size_t i=0; i<J.size(); i++) {
			double T=J[i].Period();
                        double angle, phasesplit;
                        J[i].MPS(entrada.rot(),angle, phasesplit);
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<angle
				<<setw(entrada.width())<<setprecision(entrada.precision())<<phasesplit
				<<endl;
                }
		break;
	case CommandLine::skewness:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double skew1=J[i].kappa();
			double strike1=J[i].strike1();
			double skew2=J[i].eta();
			double strike2=J[i].strike2();
			double alfa1=J[i].alfa(1);
			double alfa2=J[i].alfa(-1);
			double beta1=J[i].beta(1);
			double beta2=J[i].beta(2);
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<skew1
				<<setw(entrada.width())<<setprecision(entrada.precision())<<strike1
				<<setw(entrada.width())<<setprecision(entrada.precision())<<skew2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<strike2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<alfa1
				<<setw(entrada.width())<<setprecision(entrada.precision())<<alfa2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<beta1
				<<setw(entrada.width())<<setprecision(entrada.precision())<<beta2
				<<endl;
		}
		break;
	case CommandLine::dist_types:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double Sigma=J[i].Sigma();
			double mu=J[i].mu();
			double dD2_D2=J[i].dD2_D2();
			double kappa=J[i].kappa();
			double eta=J[i].eta();
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Sigma
				<<setw(entrada.width())<<setprecision(entrada.precision())<<mu
				<<setw(entrada.width())<<setprecision(entrada.precision())<<dD2_D2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<kappa
				<<setw(entrada.width())<<setprecision(entrada.precision())<<eta
				<<endl;
		}
		break;
	case CommandLine::ind_vectors:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double amplitude=J[i].ind_vector(0);
			double azimut=J[i].ind_vector(1);
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<azimut
				<<setw(entrada.width())<<setprecision(entrada.precision())<<amplitude
				<<endl;
		}
		break;
	}
/*	cout<<"#        T             RHO         RHOMIN        RHOMAX"
	    <<"        PHI           PHIMIN        PHIMAX"<<endl;
	cout<<"#      T      #    KAPPA    #    SIGMA    #      MU     #    dD2/D2   "
	    <<"#     ETA     #    STRIKE   #"<<endl;
	for(size_t i=0; i<J.size(); i++) {
		// J[i].set_prop_error(TransferFunction::gamble);
		double T=J[i].Period();
		double kappa=J[i].kappa();
		double Sigma=J[i].Sigma();
		double mu=J[i].mu();
		double dD2_D2=J[i].dD2_D2();
		double eta=J[i].eta();
		double strike=J[i].strike();
		cout<<setw(14)<<setprecision(8)<<T<<setw(14)<<setprecision(8)<<kappa \
		    <<setw(14)<<setprecision(8)<<Sigma<<setw(14)<<setprecision(8)<<mu \
		    <<setw(14)<<setprecision(8)<<dD2_D2<<setw(14)<<setprecision(8)<<eta \
		    <<setw(14)<<setprecision(8)<<strike<<endl;
	}
*/
	
	return EXIT_SUCCESS;
}
