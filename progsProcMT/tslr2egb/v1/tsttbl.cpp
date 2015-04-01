#include <iostream>
#include <fstream>
#include <string>
#include "readtbl.h"

int main(int argc, char* argv[]) {
    string tblstr=argv[1];
    table tbl=readtbl(tblstr);

    cout<<"Serial Number of MTU unit      :"<<tbl.snum<<endl;
    cout<<"Number of samples per minute   :"<<tbl.srpm<<endl;
    cout<<"Site ID of the measurement     :"<<tbl.site<<endl;
    cout<<"Survey ID                      :"<<tbl.srvy<<endl;
    cout<<"Current data acquisition window:"<<tbl.awin<<endl;
    cout<<"Gain for E channels. 3,12,48   :"<<tbl.egn<<endl;
    cout<<"Gain for H channels. 2, 6, 8   :"<<tbl.hgn<<endl;
    cout<<"Ex sensor azimuth              :"<<tbl.eazm<<endl;
    cout<<"Ex dipole length, m            :"<<tbl.exln<<endl;
    cout<<"Ey dipole length, m            :"<<tbl.eyln<<endl;
    cout<<"Hx sensor azimuth              :"<<tbl.hazm<<endl;
    cout<<"MAG sensor base line setting   :"<<tbl.bsex<<endl;
    cout<<"MAG sensor base line setting   :"<<tbl.bsey<<endl;
    cout<<"MAG sensor base line setting   :"<<tbl.bsez<<endl;
    cout<<"Elevation altitude of site [m] :"<<tbl.elev<<endl;
    cout<<"Latitude                       :"<<tbl.latg<<endl;
    cout<<"longitude                      :"<<tbl.lngg<<endl;
    cout<<"declination                    :"<<tbl.decl<<endl;
    cout<<"Coil nominal gain, mV/nT       :"<<tbl.hnom<<endl;
}
