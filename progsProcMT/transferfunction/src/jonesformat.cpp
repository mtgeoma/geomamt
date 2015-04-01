/***************************************************************************
                          jonesformat.cpp  -  description
                             -------------------
    begin                : Sat Jan 13 2001
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

#include "inputdata.h"

using namespace DataType;

struct TFjones {
	double T;
	vector<TFunit> TF;
	TFjones() {T=0; TF.resize(6);};
};

void InputData::jones_format()
{
	// Z=> field or SI units
	// T=> adimensional
	// Z in TF is stored in field units
	const double convsi2field=1./(4.e-4*M_PI);
	
	enum index {x, y, z};
	enum TFunits {field, SI, adimensional};
	
	vector<TFjones> J;
	TFunits units;
	string line;

	// actualy do nothing with COMMENT and INFORMATION BLOCK
	while (getline(ifile,line)) {
		if(line[0]!='#' && line[0]!='>') break;
	}

	// DATA BLOCK:
	// do nothing with 1st record (yet in line): STATION NAME and AZIMUTH(?)

	// actualy only read the impedance and GDS information:
	index r,c;
	int nT=0;
	while (getline(ifile,line)) {
		string data_type;
		istringstream strline(line.c_str());
		strline>>data_type;
		
		if(data_type[0]=='T') {
			units=adimensional;
			if(line[1]=='Z') {
				r=z;
			}
			else {
				cerr<<"unknown label \""<<data_type.substr(0,3)<<'\"'<<endl;
				exit(1);
			}

			if(data_type[2]=='X')
				c=x;
			else if(data_type[2]=='Y')
				c=y;
			else {
				cerr<<"unknown label \""<<data_type.substr(0,3)<<'\"'<<endl;
				exit(1);
			}
		}
		else if(data_type[0]=='Z' || data_type[0]=='R') {
			string data_unit;
			strline>>data_unit;
			if(data_unit.compare("field")==0)
				units=field;
			else if(data_unit.compare("SI")==0)
				units=SI;
			else {  // dcmp2j didn't include unit, I hope is the case
				units=SI;
			}

			if(data_type[1]=='X')
				r=x;
			else if(data_type[1]=='Y')
				r=y;
			else {
				cerr<<"unknown label \""<<data_type.substr(0,3)<<'\"'<<endl;
				exit(1);
			}

			if(data_type[2]=='X')
				c=x;
			else if(data_type[2]=='Y')
				c=y;
			else {
				cerr<<"unknown label \""<<data_type.substr(0,3)<<'\"'<<endl;
				exit(1);
			}
		}
		else {
			cerr<<"unknown label \""<<data_type.substr(0,3)<<'\"'<<endl;
			exit(1);
		}
		
		// actualy all block data must have the same number of periods
		// and they must be the same in all block
		if(!nT) { // nT wasn't yet inicialized
			getline(ifile,line);
			istringstream istrline(line.c_str());
			istrline>>nT;
			if(nT<=0) {
				cerr<<nT<<" is a invalid number of points"<<endl;
				exit(1);
			}
			// resize de vector of transfer functions
			J.resize(nT);
		}
		else { // verify if nT wasn't modified
			int nTtmp;
			getline(ifile,line);
			istringstream istrline(line.c_str());
			istrline>>nTtmp;
			if(nTtmp!=nT) {
				cerr<<"blocks of data with different number of points ";
				cerr<<'('<<nTtmp<<"!="<<nT<<')'<<endl;
				exit(1);
			}
		}
		
		for(int i=0; i<nT; i++) {
			
			getline(ifile,line);
			
			if(data_type[0]=='Z' || data_type[0]=='T') {
                                double T, real, imag, error, weight;
				istringstream istrline(line.c_str());
				istrline>>T>>real>>imag>>error>>weight;
			
				if(T<0) T=-1./T;
				if(units==SI) {
					real*=convsi2field;
					imag*=convsi2field;
					error*=convsi2field;
				}
			
				if(J[i].T==0)
					J[i].T=T;
				else if(T!=J[i].T) {
					cerr<<"The "<<i<<"th period of "<<data_type;
					cerr<<"is different of some previous data type"<<endl;
					cerr<<'('<<T<<"!="<<J[i].T<<')'<<endl;
					exit(1);
				}
			
				J[i].TF[r*2+c].val=complex<double>(real,imag);
				J[i].TF[r*2+c].var=pow(error,2);
				J[i].TF[r*2+c].w=weight;
			}
			if(data_type[0]=='R') {
                                double T, rho, phi, rhomax, rhomin, phimax, phimin, rhow, phiw;
				istringstream istrline(line.c_str());
				istrline>>T>>rho>>phi>>rhomax>>rhomin>>phimax>>phimin>>rhow>>phiw;

                                double error;
                                complex<double> Z;
			
				if(T<0) T=-1./T;
				if(units==SI) {
				  Z=polar(sqrt(rho*4.e-7*M_PI*2.*M_PI/T),phi*(M_PI/180.))*convsi2field;
				  error=(abs(Z)*(M_PI/180.)*(phimax-phimin)/2.);
				}
				else { // field unit
				  Z=polar(sqrt(rho/(0.2*T)),phi*(M_PI/180.));
				  error=(abs(Z)*(M_PI/180.)*(phimax-phimin)/2.);
				}

			
				if(J[i].T==0)
					J[i].T=T;
				else if(T!=J[i].T) {
					cerr<<"The "<<i<<"th period of "<<data_type;
					cerr<<"is different of some previous data type"<<endl;
					cerr<<'('<<T<<"!="<<J[i].T<<')'<<endl;
					exit(1);
				}
			
				if(J[i].TF[r*2+c].val==0.) {
				    J[i].TF[r*2+c].val=Z;
				    J[i].TF[r*2+c].var=pow(error,2);
				    J[i].TF[r*2+c].w=(rhow+phiw)/2.;
                                }
			}
		}
	}
	
	if(!ifile.eof())
		cerr<<"Warning! The file \""<<ifile_name<<"\" was not read until the end"<<endl;
	
	for(size_t i=0;i<J.size();i++) {
		TransferFunction TF(J[i].T, J[i].TF);
		VTF.push_back(TF);
	}
}
