/***************************************************************************
                          jones_format.cpp  -  description
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

#include "imputdata.h"

//Apagado pelo KDevelop: void ImputData::jones_format()
//Apagado pelo KDevelop: {
//Apagado pelo KDevelop: 	string line;
//Apagado pelo KDevelop: 		
//Apagado pelo KDevelop: 	// actualy do nothing with COMMENT and INFORMATION BLOCK
//Apagado pelo KDevelop: 	while (getline(ifile,line)) {
//Apagado pelo KDevelop: 		if(line[0]!='#' && line[0]!='>') break;
//Apagado pelo KDevelop: 	}
//Apagado pelo KDevelop: 	
//Apagado pelo KDevelop: 	// DATA BLOCK:
//Apagado pelo KDevelop: 	// do nothing with 1st record: STATION NAME and AZIMUTH(?)
//Apagado pelo KDevelop: 	getline(ifile,line);
//Apagado pelo KDevelop: 	
//Apagado pelo KDevelop: 	// actualy only read the impedance and GDS information:
//Apagado pelo KDevelop: 	index r,c;
//Apagado pelo KDevelop: 	int nT=0;
//Apagado pelo KDevelop: 	while (getline(ifile,line)) {
//Apagado pelo KDevelop: 		string data_type;
//Apagado pelo KDevelop: 		istrstream strline(line.c_str());
//Apagado pelo KDevelop: 		strline>>data_type;
//Apagado pelo KDevelop: 		
//Apagado pelo KDevelop: 		switch (data_type[0]) {
//Apagado pelo KDevelop: 		case 'Z':
//Apagado pelo KDevelop: 			string unit;
//Apagado pelo KDevelop: 			strline>>unit;
//Apagado pelo KDevelop: 			
//Apagado pelo KDevelop: 			switch(data_type[1]) {
//Apagado pelo KDevelop: 			case 'X':
//Apagado pelo KDevelop: 				r=x;
//Apagado pelo KDevelop: 				break;
//Apagado pelo KDevelop: 			case 'Y':
//Apagado pelo KDevelop: 				r=y;
//Apagado pelo KDevelop: 				break;
//Apagado pelo KDevelop: 			default:
//Apagado pelo KDevelop: 				cerr<<"unknown label \""<<substr(data_type,0,3)<<'\"'<<endl;
//Apagado pelo KDevelop: 				exit(1);
//Apagado pelo KDevelop: 			}
//Apagado pelo KDevelop: 			switch(data_type[2]) {
//Apagado pelo KDevelop: 			case 'X':
//Apagado pelo KDevelop: 				c=x;
//Apagado pelo KDevelop: 				break;
//Apagado pelo KDevelop: 			case 'Y':
//Apagado pelo KDevelop: 				c=y;
//Apagado pelo KDevelop: 				break;
//Apagado pelo KDevelop: 			default:
//Apagado pelo KDevelop: 				cerr<<"unknown label \""<<substr(data_type,0,3)<<'\"'<<endl;
//Apagado pelo KDevelop: 				exit(1);
//Apagado pelo KDevelop: 			}
//Apagado pelo KDevelop: 			break;
//Apagado pelo KDevelop: 		case 'T':
//Apagado pelo KDevelop: 			if(line[1]=='Z') {
//Apagado pelo KDevelop: 				r=z;
//Apagado pelo KDevelop: 			}
//Apagado pelo KDevelop: 			else {
//Apagado pelo KDevelop: 				cerr<<"unknown label \""<<substr(line,0,3)<<'\"'<<endl;
//Apagado pelo KDevelop: 				exit(1);
//Apagado pelo KDevelop: 			}
//Apagado pelo KDevelop: 				
//Apagado pelo KDevelop: 			switch(line[2]) {
//Apagado pelo KDevelop: 			case 'X':
//Apagado pelo KDevelop: 				c=x;
//Apagado pelo KDevelop: 				break;
//Apagado pelo KDevelop: 			case 'Y':
//Apagado pelo KDevelop: 				c=y;
//Apagado pelo KDevelop: 				break;
//Apagado pelo KDevelop: 			default:
//Apagado pelo KDevelop: 				cerr<<"unknown label \""<<substr(line,0,3)<<'\"'<<endl;
//Apagado pelo KDevelop: 				exit(1);
//Apagado pelo KDevelop: 			}
//Apagado pelo KDevelop: 		}
//Apagado pelo KDevelop: 		
//Apagado pelo KDevelop: 		// actualy all block data must have the same number of periods
//Apagado pelo KDevelop: 		if(!nT) { // nT wasn't yet inicialized
//Apagado pelo KDevelop: 			getline(ifile,line);
//Apagado pelo KDevelop: 			istrstream istrline(line.c_str());
//Apagado pelo KDevelop: 			istrline>>nT;
//Apagado pelo KDevelop: 			if(nT<=0) {
//Apagado pelo KDevelop: 				cerr<<nT<<" is a invalid number of points"<<endl;
//Apagado pelo KDevelop: 				exit(1);
//Apagado pelo KDevelop: 			}
//Apagado pelo KDevelop: 			// resize de vector of transfer functions
//Apagado pelo KDevelop: 			TF.resize(nT);
//Apagado pelo KDevelop: 		}
//Apagado pelo KDevelop: 		
//Apagado pelo KDevelop: 		{ // to do nTtmp and istrline local
//Apagado pelo KDevelop: 		int nTtmp;
//Apagado pelo KDevelop: 		getline(ifile,line);
//Apagado pelo KDevelop: 		istrstream istrline(line.c_str(),line.size());
//Apagado pelo KDevelop: 		istrline>>nTtmp;
//Apagado pelo KDevelop: 		if(nTtmp!=nT) {
//Apagado pelo KDevelop: 			cerr<<"blocks of data with different number of points ";
//Apagado pelo KDevelop: 			cerr<<'('<<nTtmp<<"!="<<nT<<\')'<<endl;
//Apagado pelo KDevelop: 			exit(1);
//Apagado pelo KDevelop: 		}
//Apagado pelo KDevelop: 		}
//Apagado pelo KDevelop: 		
//Apagado pelo KDevelop: 		for(int i=nT; i<nT; i++) {
//Apagado pelo KDevelop: 			getline(ifile,line);
//Apagado pelo KDevelop: 			istrstream istrline(line.c_str());
//Apagado pelo KDevelop: 		
//Apagado pelo KDevelop: 		}
//Apagado pelo KDevelop: 	}
//Apagado pelo KDevelop: 	
//Apagado pelo KDevelop: }
