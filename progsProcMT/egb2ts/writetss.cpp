#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include "boost/date_time/posix_time/posix_time.hpp"

using namespace std;
using namespace boost::gregorian;
using namespace boost::posix_time;

void writets(string instr, string outstr, vector<double> conv, double cda, double cdb) {
    // try to open the input file
    ifstream in(instr.c_str());
    if(in.fail()) {
        cerr<<"can't open "<<instr<<" file"<<endl;
        exit(1);
    }
    
    // try to open the .clk file
    string line;
    double sr;
    int year, month, day, hour, min, sec, micro;
    
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
    
    // set start time
    ptime time_t_epoch(date(1970,1,1));
    date start_date(year,month,day);
    ptime smp_time(start_date,hours(hour)+minutes(min)+seconds(sec)+microseconds(0));

    // decompose sample rate (sr) into seconds and microseconds
    sec=int(sr);
    micro=int((sr-sec)*1.0e6);

    // try to open the output file
    ofstream out(outstr.c_str());
    if(out.fail()) {
        cerr<<"can't open "<<outstr<<" file"<<endl;
        exit(1);
    }
    
    // write data
    int nchan=conv.size();
    while(getline(in,line)) {
        time_duration diff =  smp_time - time_t_epoch;
	out.precision(6);
        out << fixed << diff.total_microseconds()/1.0e6;

        istringstream iline(line.c_str());
        for(int i=0;i<nchan;i++) {
            double ch;
            iline>>ch;
            ch*=conv[i];
            out.precision(4);
            out << "\t" << scientific << ch;
        }
        out << "\n";
        smp_time+=seconds(sec)+microseconds(micro);
    }
    out.close();
    in.close();
}
