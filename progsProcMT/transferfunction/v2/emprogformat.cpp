/***************************************************************************
                          emprogformat.cpp  -  description
                             -------------------
    begin                : Sun Jan 14 2001
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

#include "imputdata.h"

using namespace DataType;

struct TFemprog {
	double T;
	vector<TFunit> TF;
	TFemprog() {T=0; TF.resize(6);};
};

struct Ztmp{
	double T;
	double real;
	double imag;
	double var;
	Ztmp() {T=0;};
};

void ImputData::emprog_format()
{
	//const double convsi2field=1./(4.e-4*M_PI);
	const double convsi2field=1.;
	
	enum index {x, y, z};
	
	vector<TFemprog> E;
	string line;
	
	index r,c;
	int nT=0;
	while (getline(ifile,line)) {
		// skip lines until find "IMPEDANCE TENSOR"
		if(line.find("IMPEDANCE TENSOR")==string::npos) continue;
		
		// read a block of data
		getline(ifile,line);
		if(line.find("Zxx")!=string::npos) {
			r=x; c=x;
		}
		else if(line.find("Zxy")!=string::npos) {
			r=x; c=y;
		}
		else if(line.find("Zyx")!=string::npos) {
			r=y; c=x;
		}
		else if(line.find("Zyy")!=string::npos) {
			r=y; c=y;
		}
		else {
			cerr<<"Can't recognize any impedance tensor component in line:"<<endl;
			cerr<<line<<endl;
			exit(1);
		}
		// skip a blank line
		getline(ifile,line);
		
		// impedance tensor block data format:
		// T(s)  dT  Re(Zxx)  dZxx  T(s)  dT  Im(Zxx)  dZxx
		// It's assumed (and not checked) that 2nd T(s) and dZxx are iqual
		// the 1st and that dT==0.
		// Obs: var(Zreal)=var(Zimag)=var(Z)/2
		
		// store one Z componet in a vector
		
		vector<Ztmp> VZcomp;
		
		int nTtmp=0;
		for(;;) {
			double T, dT, var;
			Ztmp Zcomp;
			string tmpline;
			getline(ifile,tmpline);
			istringstream istrline(tmpline.c_str());
			istrline>>Zcomp.T>>dT>>Zcomp.real>>Zcomp.var>>T>>dT>>Zcomp.imag>>var;
			
			Zcomp.real*=convsi2field;
			Zcomp.imag*=convsi2field;
			Zcomp.var=2.*pow(convsi2field*Zcomp.var,2);
			
			// fragile condiction: assume that if T==0, line is a blank line and
			// then is the end of block data.
			if(Zcomp.T==0) break;
			if(nT==0) { // if is the 1st block
				VZcomp.push_back(Zcomp);
			}
			else if(nTtmp<nT) {
				if(E[nTtmp].T!=Zcomp.T) {
					cerr<<nTtmp<<"th period isn't the same in all block data"<<endl;
					exit(1);
				}
				E[nTtmp].TF[r*2+c].val=complex<double>(Zcomp.real,Zcomp.imag);
				E[nTtmp].TF[r*2+c].var=Zcomp.var;
				nTtmp++;
			}
			else {
				cerr<<"block data greater than some previous block"<<endl;
				exit(1);
			}
		}
		
		if(nT==0) { // if is the 1st block
			nT=VZcomp.size();
			E.resize(nT);
			for(size_t i=0; i<E.size(); i++) {
				E[i].T=VZcomp[i].T;
				E[i].TF[r*2+c].val=complex<double>(VZcomp[i].real,VZcomp[i].imag);
				E[i].TF[r*2+c].var=VZcomp[i].var;
			}
		}
		else if(nTtmp!=nT) {
			cerr<<"block with different number of data"<<endl;
			exit(1);
		}
	}
	
	if(!ifile.eof())
		cerr<<"Warning! The file \""<<ifile_name<<"\" was not read until the end"<<endl;
	
	for(size_t i=0;i<E.size();i++) {
		TransferFunction TF(E[i].T, E[i].TF);
		VTF.push_back(TF);
	}
}
