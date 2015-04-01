#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include "readhead.h"
#include "hardware.h"

void writesp(string spstr, table tbl, vector<calibration> hardware) {
    // open file
    ofstream sp(spstr.c_str());
    if(sp.fail()) {
        cout<<"can't open the sp file "<<spstr<<endl;
        exit(1);
    }
    cout<<"opened output sp file "<<spstr<<endl;

    vector<double> conv(tbl.nchan);
    for(int i=0;i<3;i++)
        conv[i]=(hardware[i].gain*double(tbl.fs))/(tbl.magfac*tbl.gain[i]);
    for(int i=3;i<5;i++)
        conv[i]=-1.*(hardware[i].gain*double(tbl.fs))/(tbl.itgain*tbl.gain[i]);

    // write the *.sp file

    // station/run ID
    sp<<tbl.srvy<<tbl.site<<tbl.awin<<endl;
    // station coordinates: lat/lon with decimal deg. frac
    sp<<tbl.latg<<"    "<<tbl.lngg<<endl;
    // geomagnetic declination of site
    const ios_base::fmtflags old_opt = sp.flags();
    streamsize prec = sp.precision();
    sp<<fixed<<setprecision(0)<<tbl.decl<<setprecision(prec)<<endl;
    sp.flags(old_opt);
    // number of channels
    sp<<tbl.nchan<<endl;
    // sampling rate on seconds
    sp<<60./tbl.srpm<<endl;
    // clock offset & linear drift coefficients (cda, cdb)
    sp<<"0.  0."<<endl;

    // channel id
    sp<<"Hx"<<endl;
    // sensor orientation; deg. E. of geomag N; vert. tilt
    sp<<tbl.hxaz<<"  0."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<conv[0]<<"  1"<<endl;
    sp<<"LR"<<endl;
    sp<<hardware[0].lpass<<endl;

    // channel id
    sp<<"Hy"<<endl;
    // sensor orientation; deg. E. of geomag N; vert. tilt
    sp<<tbl.hxaz+90.<<"  0."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<conv[1]<<"  1"<<endl;
    sp<<"LR"<<endl;
    sp<<hardware[1].lpass<<endl;

    // channel id
    sp<<"Hz"<<endl;
    // sensor orientation; deg. E. of geomag N; vert. tilt
    sp<<"0.  0."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<conv[2]<<"  1"<<endl;
    sp<<"LR"<<endl;
    sp<<hardware[2].lpass<<endl;

    // channel id
    sp<<"Ex"<<endl;
    // electrode line length (m), angle, tilt, amp gain
    sp<<tbl.exln<<"  "<<tbl.exaz<<"  0.  1."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<conv[3]<<"  2"<<endl;
    sp<<"LR"<<endl;
    sp<<hardware[3].lpass<<endl;
    sp<<"HR"<<endl;
    sp<<hardware[3].hpass<<endl;

    // channel id
    sp<<"Ey"<<endl;
    // electrode line length (m), angle, tilt, amp gain
    sp<<tbl.eyln<<"  "<<tbl.eyaz<<"  0.  1."<<endl;
    // count convertion (nT/cont), number of filters
    sp<<conv[4]<<"  2"<<endl;
    sp<<"LR"<<endl;
    sp<<hardware[4].lpass<<endl;
    sp<<"HR"<<endl;
    sp<<hardware[4].hpass<<endl;

    sp.close();
}
