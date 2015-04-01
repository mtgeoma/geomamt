#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "gds.h"

vector<tipper> read_egbert (string instr){
    // open input file
    ifstream ifile(instr.c_str());
    if(ifile.fail()) {
        cout<<"can't open the "<<instr<<" file"<<endl;
        exit(1);
    }

    string line;
    size_t ln;
    // read header
	
    // look for number of channels and frequencies
    while (getline(ifile,line))
	if((ln=line.find("number of channels "))!=string::npos) break;
    string channels(line,ln+19,3);
    int nch=atoi(channels.c_str());

    string nfrequencies(line,line.find("number of frequencies ")+22,4);
    int nT=atoi(nfrequencies.c_str());
    if(nT<=0) {
	cerr<<"reading number of frequencies"<<endl;
	exit(1);
    }
	
	
    // check if the first 2 channels are Hx and Hy
    getline(ifile,line); // skip hedear line: "orientations and tilts of each channel"
    getline(ifile,line);
    if((line.find("Hx"))==string::npos) {
        cerr<<"first channel isn't Hx"<<endl;
        exit(1);
    }
    getline(ifile,line);
    if((line.find("Hy"))==string::npos) {
        cerr<<"second channel isn't Hy"<<endl;
        exit(1);
    }


    // loock for Hz channel
    int hzch;
    for(hzch=1;hzch<=nch-2;hzch++) {
        getline(ifile,line);
        if((line.find("Hz"))!=string::npos)
            break;  
    }
    // check if found Hz channel
    if((line.find("Hz"))==string::npos){
        cerr<<"couldn't find Hz channel in \"orientations and tilts of each channel\""<<endl;
        exit(1);
    }

    // start to read TF, S and N
    vector<tipper> VecTz;

    while (getline(ifile,line)) {
        tipper Tz;
        if((ln=line.find("period :"))==string::npos) continue;

        string strperiod(line,ln+8,13);
        if((Tz.T=atof(strperiod.c_str()))<=0) {
            cerr<<VecTz.size()+1<<" element with T<=0"<<endl;
            exit(1);
        }

        while(getline(ifile,line))
            if(line.find("Transfer Functions")!=string::npos) break;

        // get the tipper parameters
        for(int ch=1;ch<=nch-2;ch++) {
            getline(ifile,line);
            if(ch==hzch) {
                double Tzxre, Tzxim, Tzyre, Tzyim;
                istringstream istrline(line.c_str());
                istrline>>Tzxre>>Tzxim>>Tzyre>>Tzyim;
                Tz.x=complex<double>(Tzxre,Tzxim);
                Tz.y=complex<double>(Tzyre,Tzyim);
            }
        }

        getline(ifile,line);
        if(line.find("Signal Power Matrix")==string::npos) {
            cerr<<"error reading Signal Power Matrix of period "<<Tz.T<<endl;
            exit(1);
        }

        // Power Matrix
        complex<double> S[3];
        for (size_t i=0;i<2;i++) {
            getline(ifile,line);
            istringstream istrline(line.c_str());
            for (size_t j=0;j<=i;j++) {
                double real, imag;
                istrline>>real>>imag;
                S[i*(i+1)/2+j]=complex<double>(real,imag);
            }
        }

        getline(ifile,line);
        if(line.find("Residual Covar")==string::npos) {
            cerr<<"error reading Residual Covariation of period "<<Tz.T<<endl;
            exit(1);
        }

        // Residual Covariation
        vector<complex<double> > N((nch-2)*(nch-1)/2);
        for (size_t i=0;i<nch-2;i++) {
            getline(ifile,line);
            istringstream istrline(line.c_str());
            for (size_t j=0;j<=i;j++) {
                double real, imag;
                istrline>>real>>imag;
                N[i*(i+1)/2+j]=complex<double>(real,imag);
            }
        }
        Tz.xvar=N[hzch*(hzch+1)/2-1].real()*S[0].real();
        Tz.yvar=N[hzch*(hzch+1)/2-1].real()*S[2].real();
        VecTz.push_back(Tz);
    }
    if(!ifile.eof())
        cerr<<"Warning! The file \""<<instr<<"\" was not read until the end"<<endl;
    ifile.close();
    return VecTz;
}
