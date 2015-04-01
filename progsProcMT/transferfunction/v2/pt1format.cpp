/***************************************************************************
                          pt1format.cpp  -  description
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
#include <iomanip>

using namespace DataType;

//struct TFemprog {
//	double T;
//	vector<TFunit> TF;
//	TFemprog() {T=0; TF.resize(6);};
//};

//vector<double> getpt1data(const char* target, int& first, int& last);

void ImputData::pt1_format()
{
    int first=-1,  // variables to control the first and last data occurrence,
        last=-1;   // must be the same for all parameters

    vector<double> ReZxx=getpt1data("Zxx Re",first,last);
    vector<double> ImZxx=getpt1data("Zxx Im",first,last);
    vector<double> VarZxx=getpt1data("Zxx VAR",first,last);

    vector<double> ReZxy=getpt1data("Zxy Re",first,last);
    vector<double> ImZxy=getpt1data("Zxy Im",first,last);
    vector<double> VarZxy=getpt1data("Zxy VAR",first,last);

    vector<double> ReZyx=getpt1data("Zyx Re",first,last);
    vector<double> ImZyx=getpt1data("Zyx Im",first,last);
    vector<double> VarZyx=getpt1data("Zyx VAR",first,last);

    vector<double> ReZyy=getpt1data("Zyy Re",first,last);
    vector<double> ImZyy=getpt1data("Zyy Im",first,last);
    vector<double> VarZyy=getpt1data("Zyy VAR",first,last);

    vector<double> log10frq=getpt1data("LOG10(Frq) Hz",first,last);

    for(size_t i=0;i<log10frq.size();i++) {
        double pt1T=pow(10,-log10frq[i]);
        vector<TFunit> pt1TF;
        pt1TF.resize(6);

        pt1TF[0].val=complex<double>(ReZxx[i],ImZxx[i]);
        pt1TF[0].var=VarZxx[i];
        pt1TF[1].val=complex<double>(ReZxy[i],ImZxy[i]);
        pt1TF[1].var=VarZxy[i];
        pt1TF[2].val=complex<double>(ReZyx[i],ImZyx[i]);
        pt1TF[2].var=VarZyx[i];
        pt1TF[3].val=complex<double>(ReZyy[i],ImZyy[i]);
        pt1TF[3].var=VarZyy[i];

        TransferFunction TF(pt1T, pt1TF);
	VTF.push_back(TF);
    }
}

vector<double> ImputData::getpt1data(const char* target, int& pfirst, int& plast) {
    // put at the begin of pt1 file
    ifile.seekg(0);

    int first=-1;
    int last=-1;

    // look for the line with the target string
    string line;
    while (getline(ifile,line))
        if(strstr(line.c_str(),target)!=0) break;

    if(ifile.eof()) {
        cerr<<"can't find target \""<<target<<"\""<<endl;
        exit(1);
    }

    // get the data until find the next parameter
    vector<double> data;
    int n=0;
    while(getline(ifile,line)) {
        if(strstr(line.c_str(),"PARAMETER")!=0) break;

        istringstream strline(line.c_str());
        for(int i=0; i<5; i++) {
            double datum;
            strline>>datum;
            if(datum!=0.) {
                if(first<0) first=n;
                if(last>=0) cerr<<"skipt data no."<<n<<endl;
                data.push_back(datum);
            }
            else if(last<0 && first>=0) last=n;
            n++;
        }
    }

    if(pfirst<0 && plast<0) { // if is the first parameter
        pfirst=first;
        plast=last;
    }// compare first and last occurrence with the previous data
    else if(pfirst!=first || plast!=plast) {
        cerr<<target<<" and previous parameter didn't have the same data range:"<<endl;
        cerr<<setw(14)<<target<<": first: "<<setw(4)<<first<<" last: "<<setw(4)<<last<<endl;
        cerr<<"Previous      : first: "<<setw(4)<<pfirst<<" last: "<<setw(4)<<plast<<endl;
        exit(1);
    }

    return data;
}
