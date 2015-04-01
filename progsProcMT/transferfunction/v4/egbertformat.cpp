/***************************************************************************
                          egbertformat.cpp  -  description
                             -------------------
    begin                : Wed Jan 17 2001
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
	
struct TFegbert {
	double T;
	vector<TFunit> TF;
	TFegbert() {T=0; TF.resize(6);};
};

void InputData::egbert_format()
{
	long nT;
	vector<TFegbert> E;
	
	string line;
	size_t ln;
	// read header
	// skip some lines
	while (getline(ifile,line))
		if((ln=line.find("number of channels "))!=string::npos) break;
	
	// check if no. of channels==5
	string channels(line,ln+19,3);
	if(atoi(channels.c_str())!=5) {
		cerr<<"actualy this program is restricted to 5 channels"<<endl;
		exit(1);
	}
	
	string nfrequencies(line,line.find("number of frequencies ")+22,4);
	if((nT=atoi(nfrequencies.c_str()))<=0) {
		cerr<<"reading number of frequencies"<<endl;
		exit(1);
	}
	
	E.resize(nT);

    // start to read TF, S and N
	size_t n=0;
	while (getline(ifile,line)) {
		if((ln=line.find("period :"))==string::npos) continue;
		
		string strperiod(line,ln+8,12);
		if((E[n].T=atof(strperiod.c_str()))<=0) {
			cerr<<n+1<<" element with T<=0"<<endl;
			exit(1);
		}
		
		while(getline(ifile,line))
			if(line.find("Transfer Functions")!=string::npos) break;
		
		for(size_t i=0; i<3; i++) {
			double TFxre, TFxim, TFyre, TFyim;
			getline(ifile,line);
			istringstream istrline(line.c_str());
			istrline>>TFxre>>TFxim>>TFyre>>TFyim;
			E[n].TF[i*2].val=complex<double>(TFxre,TFxim);
			E[n].TF[i*2+1].val=complex<double>(TFyre,TFyim);
		}
		
		getline(ifile,line);
		if(line.find("Signal Power Matrix")==string::npos) {
			cerr<<"reading data of period "<<E[n].T<<endl;
			exit(1);
		}
		
		// actually only store the var of each element of the impedance tensor
		// Power Matrix
		complex<double> S[3];
		for (size_t i=0;i<2;i++) {
			getline(ifile,line);
			istringstream istrline(line.c_str());
	    	for (size_t j=0;j<=i;j++) {
				double real, imag;
				istrline>>real>>imag;
				S[i+j]=complex<double>(real,imag);
			}
		}
		
		getline(ifile,line);
		if(line.find("Residual Covar")==string::npos) {
			cerr<<"reading data of period "<<E[n].T<<endl;
			exit(1);
		}
		
		complex<double> N[6];
		for (size_t i=0;i<3;i++) {
			getline(ifile,line);
			istringstream istrline(line.c_str());
			for (size_t j=0;j<=i;j++) {
				double real, imag;
				istrline>>real>>imag;
				N[i+j+i/2]=complex<double>(real,imag);
			}
		}
		// calculed the variances
		for(size_t i=0;i<3;i++) {
			for(size_t j=0;j<2;j++)
				E[n].TF[i*2+j].var=N[i*2+i/2].real()*S[j*2].real();
		}
		// reorganize the TF matrix, originally it is:
		// Tzx Tzy
		// Zxx Zxy
		// Zyx Zyy
		TFunit Tzx=E[n].TF[0];
		TFunit Tzy=E[n].TF[1];
		E[n].TF[0]=E[n].TF[2];
		E[n].TF[1]=E[n].TF[3];
		E[n].TF[2]=E[n].TF[4];
		E[n].TF[3]=E[n].TF[5];
		E[n].TF[4]=Tzx;
		E[n].TF[5]=Tzy;
		
		n++;
	}
	if(!ifile.eof())
		cerr<<"Warning! The file \""<<ifile_name<<"\" was not read until the end"<<endl;
	
	for(size_t i=0;i<E.size();i++) {
		TransferFunction TF(E[i].T, E[i].TF);
		VTF.push_back(TF);
	}
}
