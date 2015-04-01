#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <inttypes.h>
#include <byteswap.h>
#include "count.h"

void writeegb(string instr, string outstr,timecount clock) {
    // write the clock file
    string clkstr=outstr.substr(0,outstr.find_last_of("."));
    clkstr+=".clk";
    ofstream clk(clkstr.c_str());
    if(clk.fail()) {
        cout<<"can't open the output clock file "<<clkstr<<endl;
        exit(1);
    }
    else {
        int year, month, day;
        jul_to_cal(clock.sdata.jul,&year,&month,&day);
        // write sample rate in seconds
        clk.setf(ios::fixed);
        clk.precision(5);
        clk<<60./clock.srpm<<endl;
        // write start time
        clk<<" "<<setw(2)<<setfill('0')<<year%100
           <<" "<<setw(2)<<setfill('0')<<month
           <<" "<<setw(2)<<setfill('0')<<day
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.min/60
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.min%60
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.sec<<endl;
        // write reference time (begining of year)
        clk<<" "<<setw(2)<<setfill('0')<<year%100
           <<" "<<setw(2)<<setfill('0')<<1
           <<" "<<setw(2)<<setfill('0')<<1
           <<" "<<setw(2)<<setfill('0')<<0
           <<" "<<setw(2)<<setfill('0')<<0
           <<" "<<setw(2)<<setfill('0')<<0<<endl;
        clk.close();
    }

    // open input MTU-LR file
    ifstream in(instr.c_str(), ios::in | ios::binary);
    if(in.fail()) {
        cout<<"can't open the MTU-LR file "<<instr<<endl;
        exit(1);
    }

    // open output Egbert file
    ofstream out(outstr.c_str());
    if(out.fail()) {
        cout<<"can't open the Egbert file "<<outstr<<endl;
        exit(1);
    }

    // put in the start position
    in.seekg(clock.sizehdr);

    for(long i=0; i<clock.ecount; i++) {
        int16_t data[5];
        in.read(reinterpret_cast <char *> (data),sizeof(data));
        if(i>=clock.scount) {
            if(bswap)
                for(int j=0; j<5; j++)
                    data[j]=bswap_16(data[j]);
            for(int j=0; j<5; j++)
                out<<setw(10)<<data[j];
            out<<endl;
        }
        if(i==clock.imark || i==0 ) {
            in.seekg(12,ios::cur);
            if(i!=0)
                clock.imark+=10*clock.srpm;
        }
    }

    in.close();
    out.close();
}
