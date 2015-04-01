/***************************************************************************
                          ediformat.cpp  -  description
                             -------------------
    begin                : Thu Aug 16 2001
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
#include <iomanip>
#include <string>

using namespace DataType;

void ImputData::edi_format()
{
    int N=-1;  // number of data by block

    vector<double> ReZxx=getedidata(">ZXXR",N);
    vector<double> ImZxx=getedidata(">ZXXI",N);
    vector<double> VarZxx=getedidata(">ZXX.VAR",N);

    vector<double> ReZxy=getedidata(">ZXYR",N);
    vector<double> ImZxy=getedidata(">ZXYI",N);
    vector<double> VarZxy=getedidata(">ZXY.VAR",N);

    vector<double> ReZyx=getedidata(">ZYXR",N);
    vector<double> ImZyx=getedidata(">ZYXI",N);
    vector<double> VarZyx=getedidata(">ZYX.VAR",N);

    vector<double> ReZyy=getedidata(">ZYYR",N);
    vector<double> ImZyy=getedidata(">ZYYI",N);
    vector<double> VarZyy=getedidata(">ZYY.VAR",N);

    vector<double> freq=getedidata(">FREQ",N);

    for(size_t i=0;i<N;i++) {
        double ediT=1./freq[i];
        vector<TFunit> ediTF;
        ediTF.resize(6);

        ediTF[0].val=complex<double>(ReZxx[i],ImZxx[i]);
        ediTF[0].var=VarZxx[i];
        ediTF[1].val=complex<double>(ReZxy[i],ImZxy[i]);
        ediTF[1].var=VarZxy[i];
        ediTF[2].val=complex<double>(ReZyx[i],ImZyx[i]);
        ediTF[2].var=VarZyx[i];
        ediTF[3].val=complex<double>(ReZyy[i],ImZyy[i]);
        ediTF[3].var=VarZyy[i];

        TransferFunction TF(ediT, ediTF);
	VTF.push_back(TF);
    }
}

vector<double> ImputData::getedidata(const char* target, int& ndata) {
    // put at the begin of edi file
    ifile.seekg(0);

    // look for the line with the target string
    string line;
    while (getline(ifile,line))
        if(strstr(line.c_str(),target)!=0) break;

    if(ifile.eof()) {
        cerr<<"can't find target \""<<target<<"\""<<endl;
        exit(1);
    }

    // get the number of data
    string Nstr(line,line.find("//")+2,string::npos);
    int N=atoi(Nstr.c_str());

    if(ndata<0) // if is the first parameter
        ndata=N;
    else if(N!=ndata) { // compare with the previous data
        cerr<<target<<" and previous parameter didn't have the same number of data:"<<endl;
        cerr<<setw(14)<<target<<": "<<setw(4)<<N<<" data"<<endl;
        cerr<<"Previous      : "<<setw(4)<<ndata<<" data"<<endl;
        exit(1);
    }

    // get data
    vector<double> data(N);
    int Nline=N/5;  // number of line-1 to read
    int Nrest=N%5;  // number of data at the last line

    // read line with 5 data
    for(int i=0;i<Nline;i++) {
        getline(ifile,line);
        istringstream strline(line.c_str());

        int j=0;
        while(j<5 && strline>>data[i*5+j]) j++;
        if(strline.fail()) {
            cerr<<"Error reading the "<<i*5+j+1<<" data of target "<<target<<endl;
            exit(1);
        }
    }

    // read line with less than 5 data
    if(Nrest) {
        getline(ifile,line);
        istringstream strline(line.c_str());

        int j=0;
        while(j<Nrest && strline>>data[Nline*5+j]) j++;
        if(strline.fail()) {
            cerr<<"Error reading the "<<Nline*5+j+1<<" data of target "<<target<<endl;
            exit(1);
        }
    }

    return data;
}
