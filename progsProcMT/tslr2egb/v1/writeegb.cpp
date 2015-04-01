#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdio.h>
#include "bytes1.h"
#include "count.h"

int4 conv24to32 (int4 datum);

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
        clk.setf(ios::scientific);
        clk<<setprecision(4)<<60./clock.srpm<<endl;
        // write start time
        clk<<" "<<setw(2)<<setfill('0')<<year%100
           <<" "<<setw(2)<<setfill('0')<<month
           <<" "<<setw(2)<<setfill('0')<<day
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.min/60
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.min%60
           <<" 00"<<endl;
        // write reference time (same as start time)
        clk<<" "<<setw(2)<<setfill('0')<<year%100
           <<" "<<setw(2)<<setfill('0')<<month
           <<" "<<setw(2)<<setfill('0')<<day
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.min/60
           <<" "<<setw(2)<<setfill('0')<<clock.sdata.min%60
           <<" 00"<<endl;
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

    // size of a record
    unsigned long recsize=clock.tagsize+clock.nchan*clock.srpm*3u;

    // put in the start position
    in.seekg(clock.scount*recsize);

    for(long i=clock.scount; i<clock.ecount; i++) {
        amxtds time;
        if(!in.read(&time,sizeof(amxtds))) {
            cout<<"error reading time at record number "<<i<<endl;
            exit(1);
        }

        if(time.cen) {  // just check if is the expected time
            long jul=cal_to_jul(time.yr+100*time.cen,time.mn,time.dy);
            long min=time.min+time.hr*60;
            // calculate the actual expected record number
            long rec=(jul-clock.swin.jul)*1440L+(min-clock.swin.min);
            if(rec!=i)
                cout<<"count rec: "<<i<<" expected rec: "<<rec<<endl;
        }
        // positioning to read data
        in.seekg(clock.tagsize-sizeof(amxtds),ios::cur);
        //cout<<clock.scount<<"*"<<recsize<<"+"<<clock.tagsize<<"="
        //    <<clock.scount*recsize+clock.tagsize<<endl;
        //in.seekg(clock.scount*recsize+clock.tagsize);
        for(int j=0;j<clock.srpm;j++) {
            vector<int4> data(clock.nchan);
            for(int k=0;k<clock.nchan;k++) {
                int4 datum;
                in.read(&datum,3);
                //cout<<setw(10)<<hex<<datum;
                data[k]=conv24to32(datum);
            }
            //cout<<endl;
            out<<setw(10)<<data[2]
               <<setw(10)<<data[3]
               <<setw(10)<<data[4]
               <<setw(10)<<data[0]
               <<setw(10)<<data[1]<<endl;
        }
    }
    in.close();
    out.close();
}

int4 conv24to32 (int4 datum) {
    if(bswap)
        datum=bswap_32(datum);

    datum=datum<<8;
    datum/=256;
    return datum;
}
