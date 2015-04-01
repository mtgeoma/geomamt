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

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include<stdio.h>
#include "transferfunction.h"
#include "inputdata.h"
#include "commandline.h"

// argumentos: arquivo tipo extrator
int main(int argc, char *argv[])
{
	CommandLine entrada(argc, argv);
	DataType::InputData j_data(entrada.file(),entrada.file_format());
	vector<TransferFunction> J;
	j_data.getdata(J);

	if(entrada.rot()!=0. && entrada.extract()!=CommandLine::phasesplit) {
		for(size_t i=0; i<J.size(); i++)
			J[i].rotate(entrada.rot());
	}

	switch (entrada.extract()) {
	case CommandLine::impedance:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			complex<double> Z=J[i].val(entrada.row(),entrada.col());
			double Zerr=J[i].err(entrada.row(),entrada.col());
			double weight=J[i].weight(entrada.row(),entrada.col());

                        cout.setf(ios::scientific);
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
                        cout.setf(ios::scientific);
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
		printf("%-14ld%6.4f\n",J.size(),1.0);
		for(int i=0; i<J.size(); i++) {
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
                        cout.setf(ios::scientific);
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<angle
				<<setw(entrada.width())<<setprecision(entrada.precision())<<phasesplit
				<<endl;
                }
		break;
	case CommandLine::strikes:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double strike1=J[i].strike1();
			double strike2=J[i].strike2();
			double alfa=J[i].alfa();
                        cout.setf(ios::scientific);
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<strike1
				<<setw(entrada.width())<<setprecision(entrada.precision())<<strike2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<alfa
				<<endl;
		}
		break;
	case CommandLine::dist_types:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double kappa=J[i].kappa();
			double Sigma=J[i].Sigma();
			double mu=J[i].mu();
			double dD2_D2=J[i].dD2_D2();
			double eta=J[i].eta();
			double beta1=J[i].beta(1); // com strike2
			double beta2=J[i].beta(2);
			double beta3=J[i].beta(3); // com alfa
			double beta4=J[i].beta(4);
			double Ekappa=J[i].kappaerr();
			double ESigma=J[i].Sigmaerr();
			double Emu=J[i].muerr();
			double Eeta=J[i].etaerr();
                        cout.setf(ios::scientific);
			cout<<setw(entrada.width())<<setprecision(entrada.precision())<<T
				<<setw(entrada.width())<<setprecision(entrada.precision())<<kappa
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Sigma
				<<setw(entrada.width())<<setprecision(entrada.precision())<<mu
				<<setw(entrada.width())<<setprecision(entrada.precision())<<dD2_D2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<eta
				<<setw(entrada.width())<<setprecision(entrada.precision())<<beta1
				<<setw(entrada.width())<<setprecision(entrada.precision())<<beta2
				<<setw(entrada.width())<<setprecision(entrada.precision())<<beta3
				<<setw(entrada.width())<<setprecision(entrada.precision())<<beta4
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Ekappa
				<<setw(entrada.width())<<setprecision(entrada.precision())<<ESigma
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Emu
				<<setw(entrada.width())<<setprecision(entrada.precision())<<Eeta
				<<endl;
		}
		break;
	case CommandLine::ind_vectors:
		for(size_t i=0; i<J.size(); i++) {
			J[i].set_prop_error(entrada.prop_error());
			double T=J[i].Period();
			double amplitude=J[i].ind_vector(0);
			double azimut=J[i].ind_vector(1);
                        cout.setf(ios::scientific);
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
                cout.setf(ios::scientific);
		cout<<setw(14)<<setprecision(8)<<T<<setw(14)<<setprecision(8)<<kappa \
		    <<setw(14)<<setprecision(8)<<Sigma<<setw(14)<<setprecision(8)<<mu \
		    <<setw(14)<<setprecision(8)<<dD2_D2<<setw(14)<<setprecision(8)<<eta \
		    <<setw(14)<<setprecision(8)<<strike<<endl;
	}
*/
	
	return EXIT_SUCCESS;
}
