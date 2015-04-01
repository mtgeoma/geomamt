#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "gds.h"
using namespace std;

vector<double> getedidata(const char* target, int& ndata, ifstream& ifile);

vector<tipper> read_edi (string instr){
    // open input file
    ifstream ifile(instr.c_str());
    if(ifile.fail()) {
        cout<<"can't open the "<<instr<<" file"<<endl;
        exit(1);
    }
    int N=-1;  // number of data by block

    vector<double> freq=getedidata(">FREQ",N,ifile);
    vector<double> TXR=getedidata(">TXR",N,ifile);
    vector<double> TXI=getedidata(">TXI",N,ifile);
    vector<double> TXVAR=getedidata(">TYVAR",N,ifile);
    vector<double> TYR=getedidata(">TYR",N,ifile);
    vector<double> TYI=getedidata(">TYI",N,ifile);
    vector<double> TYVAR=getedidata(">TYVAR",N,ifile);

    vector<tipper> Tz(N);
    for(int i=0;i<N;i++) {
        int index;
        if(freq[0]<freq[1])
	  index=N-(i+1);
        else
	  index=i;
        
        Tz[i].T=1./freq[index];
        Tz[i].x=complex<double>(TXR[index],TXI[index]);
        Tz[i].y=complex<double>(TYR[index],TYI[index]);
        Tz[i].xvar=TXVAR[index];
        Tz[i].yvar=TYVAR[index];
    }
    
    ifile.close();
    return Tz;
}

vector<double> getedidata(const char* target, int& ndata, ifstream& ifile) {
    if(ifile.fail()) {
        cout<<"problem reading data file"<<endl;
        exit(1);
    }

    // put at the begin of edi file
    ifile.seekg(0,ios::beg);

    // look for the line with the target string
    string line;
    while (getline(ifile,line))
        if(strstr(line.c_str(),target)!=0) break;

    if(ifile.eof()) {
        if((target==">TXVAR" || target==">TYVAR") && ndata>0) {
	    ifile.clear();
            vector<double> data(ndata,0.);
            return data;
        }
        else {
            cerr<<"can't find target \""<<target<<"\""<<endl;
            exit(1);
        }
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
