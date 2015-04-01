#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "gds.h"

vector<tipper> read_jones(string instr) {
    // open input file
    ifstream ifile(instr.c_str());
    if(ifile.fail()) {
        cout<<"can't open the "<<instr<<" file"<<endl;
        exit(1);
    }

    string line;
    // read header
    // actualy do nothing with COMMENT and INFORMATION BLOCK
    while (getline(ifile,line)) {
        if(line[0]!='#' && line[0]!='>') break;
    }

    // DATA BLOCK:
    // do nothing with 1st record (yet in line): STATION NAME and AZIMUTH

    int nT=0;
    vector<tipper> VTz;
    while (getline(ifile,line)) {
        string data_type;
        istringstream strline(line.c_str());
        strline>>data_type;

        if(data_type=="TZX" || data_type=="TZY") {
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
                VTz.resize(nT);
                for(int i=0; i<nT; i++) {
                    VTz[i].T=0.;
                }
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

                double T, real, imag, error, weight;
                istringstream istrline(line.c_str());
                istrline>>T>>real>>imag>>error>>weight;

                if(T<0) T=-1./T;

                if(VTz[i].T==0)
                    VTz[i].T=T;
                else if(T!=VTz[i].T) {
                    cerr<<"The "<<i<<"th period of "<<data_type;
                    cerr<<"is different of some previous data type"<<endl;
                    cerr<<'('<<T<<"!="<<VTz[i].T<<')'<<endl;
                    exit(1);
                }
                if(data_type=="TZX") {
                    VTz[i].x=complex<double>(real,imag);
                    VTz[i].xvar=pow(error,2);
                }
                if(data_type=="TZY") {
                    VTz[i].y=complex<double>(real,imag);
                    VTz[i].yvar=pow(error,2);
                }
            }
	}
    }

    if(!nT) {
        cerr<<"Error! Didn't find tipper data in "<<instr<<endl;
        exit(1);
    }

    if(!ifile.eof())
        cerr<<"Warning! The file \""<<instr<<"\" was not read until the end"<<endl;

    ifile.close();
    return VTz;
}
