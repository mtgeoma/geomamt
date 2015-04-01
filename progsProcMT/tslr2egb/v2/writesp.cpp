#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include "readtbl.h"

const double ECV=-1.E6;    // convert V/m to mV/km
const double HCV=1.E9;     // convert T to nT
const long   FS=0x7FFFFF;  // full scale (normalizing the ADN)

void writesp(string spstr, table tbl, vector<string> chstr) {
    // open file
    ofstream sp(spstr.c_str());
    if(sp.fail()) {
        cout<<"can't open the sp file "<<spstr<<endl;
        exit(1);
    }
    cout<<"opened output sp file "<<spstr<<endl;

    if(tbl.fieldtype!=1) {
        cout<<"WARNING: calibration file isn't a system calibration"<<endl
            <<"the following parameters are probably wrong:"<<endl
            <<"magnetic count convertion"<<endl
            <<"electric count convertion"<<endl
            <<"electrode line length"<<endl;
    }

    double hconv=HCV/FS;
    double econv=ECV/FS;

    // write the *.sp file

    // station/run ID
    sp<<tbl.srvy<<tbl.site<<tbl.awin<<endl;
    // station coordinates: lat/lon with decimal deg. frac
    sp<<tbl.latg<<"    "<<tbl.lngg<<endl;
    // geomagnetic declination of site
    sp<<tbl.decl<<endl;
    // number of channels
    sp<<tbl.nchan<<endl;
    // sampling rate on seconds
    sp<<60./tbl.srpm<<endl;
    // clock offset & linear drift coefficients (cda, cdb)
    sp<<"0.  0."<<endl;

    // channel id
    sp<<"Hx"<<endl;
    // sensor orientation; deg. E. of geomag N; vert. tilt
    sp<<tbl.hazm<<"  0."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<hconv<<"  1"<<endl;
    sp<<"RI"<<endl;
    sp<<chstr[2]<<endl;

    // channel id
    sp<<"Hy"<<endl;
    // sensor orientation; deg. E. of geomag N; vert. tilt
    sp<<tbl.hazm+90.<<"  0."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<hconv<<"  1"<<endl;
    sp<<"RI"<<endl;
    sp<<chstr[3]<<endl;

    // channel id
    sp<<"Hz"<<endl;
    // sensor orientation; deg. E. of geomag N; vert. tilt
    sp<<"0.  0."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<hconv<<"  1"<<endl;
    sp<<"RI"<<endl;
    sp<<chstr[4]<<endl;

    // channel id
    sp<<"Ex"<<endl;
    // electrode line length (m), angle, tilt, amp gain
    sp<<1<<"  "<<tbl.eazm<<"  0.  1."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<econv<<"  1"<<endl;
    sp<<"RI"<<endl;
    sp<<chstr[0]<<endl;

    // channel id
    sp<<"Ey"<<endl;
    // electrode line length (m), angle, tilt, amp gain
    sp<<1<<"  "<<tbl.eazm+90.<<"  0.  1."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<econv<<"  1"<<endl;
    sp<<"RI"<<endl;
    sp<<chstr[1]<<endl;

    sp.close();
}
