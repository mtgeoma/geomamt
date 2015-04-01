#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "calendar.h"

using namespace std;

void writets(string instr, string outstr, vector<double> conv, double cda, double cdb) {
    // try to open the input file
    ifstream in(instr.c_str());
    if(in.fail()) {
        cerr<<"can't open "<<instr<<" file"<<endl;
        exit(1);
    }
    
    // try to open the .clk file
    string line;
    double sr, sec;
    int year, month, day, hour, min;
    
    string clkstr=instr.substr(0,instr.find_last_of("."))+".clk";
    ifstream clk(clkstr.c_str());
    if(clk.fail()) {
        cout<<"can't open "<<clkstr<<" file"<<endl
            <<"will get clock parameters in "<<instr<<" file"<<endl;
        if(!getline(in,line)) {
            cerr<<"error reading sample rate from "<<instr<<" file"<<endl;
            exit(1);
        }
        else {
            istringstream iline(line.c_str());
            iline>>sr;
        }
        
        if(!getline(in,line)) {
            cerr<<"error reading sstart time from "<<instr<<" file"<<endl;
            exit(1);
        }
        else {
            istringstream iline(line.c_str());
            iline>>year>>month>>day>>hour>>min>>sec;
        }
        getline(in,line);  // skip reference date
    }
    else {
        if(!getline(clk,line)) {
            cerr<<"error reading sample rate from "<<clkstr<<" file"<<endl;
            exit(1);
        }
        else {
            istringstream iline(line.c_str());
            iline>>sr;
        }
        
        if(!getline(clk,line)) {
            cerr<<"error reading start time from "<<clkstr<<" file"<<endl;
            exit(1);
        }
        else {
            istringstream iline(line.c_str());
            iline>>year>>month>>day>>hour>>min>>sec;
        }
        clk.close();
    }
    
    // use year of 1950 to separate XX from XXI century
    if(year<50)
        year+=2000;
    else
        year+=1900;
    
    // set start time in julian day starting at the begining of current year
    double t0=cal_and_time_to_jul(year,month,day,hour,min,sec);
    t0-=cal_and_time_to_jul(year,1,1,0,0,0.);
    
    // set start time in seconds
    t0*=86400.;
    
    // try to open the output file
    ofstream out(outstr.c_str());
    if(out.fail()) {
        cerr<<"can't open "<<outstr<<" file"<<endl;
        exit(1);
    }
    
    // write data
    int nchan=conv.size();
    double iclk=0.; // internal clock in seconds of input file
    while(getline(in,line)) {
        double tc=iclk+(cda+iclk*cdb)+t0;
        out.setf(ios::fixed);
        out.precision(1);
        out<<setw(14)<<tc;
        
        istringstream iline(line.c_str());
        for(int i=0;i<nchan;i++) {
            double ch;
            iline>>ch;
            ch*=conv[i];
            out.precision(4);
            out<<" "<<setw(12)<<ch;
        }
        out<<endl;
        iclk+=sr;
    }
    out.close();
    in.close();
}
