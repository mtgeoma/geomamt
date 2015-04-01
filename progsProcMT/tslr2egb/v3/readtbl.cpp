#include "readtbl.h"
//#include "count.h"

table readtbl(string tblstr) {
    table tbl;

    // open file
    ifstream tblf(tblstr.c_str(), ios::in | ios::binary);
    if(tblf.fail()) {
        cout<<"can't open the tbl file "<<tblstr<<endl;
        exit(1);
    }
    cout<<"opened tbl file "<<tblstr<<endl;

    // look for serial number of MTU unit
    if(positioning(tblf,"SNUM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.snum),sizeof(tbl.snum));
        if(bswap)
            tbl.snum=bswap_32(tbl.snum);
    }
    else  // will try to find in the MTU-LR file
        tbl.snum=-1;

    // look for number of samples per minute
    if(positioning(tblf,"SRPM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.srpm),sizeof(tbl.srpm));
        if(bswap)
            tbl.srpm=bswap_32(tbl.srpm);
    }
    else  // will try to find in the MTU-LR file
        tbl.srpm=-1;

    // read site ID of the measurement
    char  SITE[9];
    if(positioning(tblf, "SITE")) {
        tblf.read(SITE,9);
        tbl.site=SITE;
    }
    else {
        tbl.site="";
        cerr<<"can't find site ID of the measurement in tbl file "<<tblstr<<endl
            <<"will use the output file name as site ID"<<endl<<endl;
    }

    // read survey ID
    char  SRVY[13];
    if(positioning(tblf, "SRVY")) {
        tblf.read(SRVY,13);
        tbl.srvy=SRVY;
    }
    else {
        tbl.srvy="";
        cerr<<"can't find site ID of the measurement in tbl file "<<tblstr<<endl;
    }

    // read current data acquisition window
    if(positioning(tblf,"AWIN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.awin),sizeof(tbl.awin));
        if(bswap)
            tbl.awin=bswap_32(tbl.awin);
    }
    else {
        tbl.awin=0;
        cerr<<"can't find data acquisition window in tbl file "<<tblstr<<endl;
    }

    // read gain for "E" channels. 3,12,48
    if(positioning(tblf, "EGN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.egn),sizeof(tbl.egn));
        if(bswap)
            tbl.egn=bswap_32(tbl.egn);
    }
    else {
        tbl.egn=3;
        cerr<<"can't find gain for E channels in tbl file "<<tblstr<<endl
            <<"will be set to 3"<<endl;
    }

    // read gain for "H" channels. 2, 6, 8
    if(positioning(tblf, "HGN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.hgn),sizeof(tbl.hgn));
        if(bswap)
            tbl.hgn=bswap_32(tbl.hgn);
    }
    else {
        tbl.hgn=2;
        cerr<<"can't find gain for H channels in tbl file "<<tblstr<<endl
            <<"will be set to 2"<<endl;
    }

    // read Ex sensor azimuth
    if(positioning(tblf, "EAZM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.eazm),sizeof(tbl.eazm));
        if(bswap)
            tbl.eazm=(float8)bswap_64((int8)tbl.eazm);
    }
    else {
        tbl.eazm=0;
        cerr<<"can't find Ex sensor azimuth in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }

    // read Ex dipole length, m
    if(positioning(tblf, "EXLN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.exln),sizeof(tbl.exln));
        if(bswap)
            tbl.exln=(float8)bswap_64((int8)tbl.exln);
    }
    else {
        tbl.exln=1.;
        cerr<<"can't find Ex dipole length in tbl file "<<tblstr<<endl
            <<"will be set to 1 "
            <<"(it will work if using system calibration file)"<<endl;
    }


    // read Ey dipole length, m
    if(positioning(tblf, "EYLN")) {
        tblf.read(reinterpret_cast <char *> (&tbl.eyln),sizeof(tbl.eyln));
        if(bswap)
            tbl.eyln=(float8)bswap_64((int8)tbl.eyln);
    }
    else {
        tbl.eyln=1.;
        cerr<<"can't find Ey dipole length in tbl file "<<tblstr<<endl
            <<"will be set to 1 "
            <<"(it will work if using system calibration file)"<<endl;
    }

    // read Hx sensor azimuth
    if(positioning(tblf, "HAZM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.hazm),sizeof(tbl.hazm));
        if(bswap)
            tbl.hazm=(float8)bswap_64((int8)tbl.hazm);
    }
    else {
        tbl.hazm=0;
        cerr<<"can't find Hx sensor azimuth in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }

    // read MAG sensor base line setting
    if(positioning(tblf, "BSEX")) {
        tblf.read(reinterpret_cast <char *> (&tbl.bsex),sizeof(tbl.bsex));
        if(bswap)
            tbl.bsex=bswap_32(tbl.bsex);
    }
    else {
        tbl.bsex=0;
        cerr<<"can't find Hx base line setting in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }

    // read MAG sensor base line setting
    if(positioning(tblf, "BSEY")) {
        tblf.read(reinterpret_cast <char *> (&tbl.bsey),sizeof(tbl.bsey));
        if(bswap)
            tbl.bsey=bswap_32(tbl.bsey);
    }
    else {
        tbl.bsey=0;
        cerr<<"can't find Hy base line setting in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }

    // read MAG sensor base line setting
    if(positioning(tblf, "BSEZ")) {
        tblf.read(reinterpret_cast <char *> (&tbl.bsez),sizeof(tbl.bsez));
        if(bswap)
            tbl.bsez=bswap_32(tbl.bsez);
    }
    else {
        tbl.bsez=0;
        cerr<<"can't find Hz base line setting in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }

    // read elevation altitude of site [m]
    int  ELEV;
    if(positioning(tblf, "ELEV")) {
        tblf.read(reinterpret_cast <char *> (&tbl.elev),sizeof(tbl.elev));
        if(bswap)
            tbl.elev=bswap_32(tbl.elev);
    }
    else {
        tbl.elev=0;
        cerr<<"can't find elevation altitude of site in tbl file "<<tblstr<<endl
            <<"will be set to 0"<<endl;
    }

    // read latitude
    char  LATG[13];
    if(positioning(tblf, "LATG")) {
        tblf.read(LATG,13);
        tbl.latg=convlat(LATG);
    }
    else {
        cerr<<"can't find latitude in tbl file "<<tblstr<<endl
            <<"enter latitude in degrees (north positive):"<<endl;
        cin>>tbl.latg;
    }

    // read longitude
    char  LNGG[13];
    if(positioning(tblf, "LNGG")) {
        tblf.read(LNGG,13);
        tbl.lngg=convlon(LNGG);
    }
    else {
        cerr<<"can't find longitude in tbl file "<<tblstr<<endl
            <<"enter longitude in degrees (east positive):"<<endl;
        cin>>tbl.lngg;
    }

    // calculate declination. Will get STIM as a reference date

    if(positioning(tblf, "STIM")) {
        amxtds  STIM;
        tblf.read(reinterpret_cast <char *> (&STIM),sizeof(amxtds));
        tbl.decl=declination(tbl.latg, tbl.lngg, tbl.elev/1e3, STIM.yr+STIM.cen*100, STIM.mn, STIM.dy);
    }
    else {
        int year, month, day;
        cout<<"enter a reference year month day to calculate declination"<<endl
            <<"(e.g.: 2001 08 12)"<<endl;
        cin>>year>>month>>day;
        tbl.decl=declination(tbl.latg, tbl.lngg, tbl.elev/1e3, year, month, day);
    }

    // read coil nominal gain, mV/nT
    if(positioning(tblf, "HNOM")) {
        tblf.read(reinterpret_cast <char *> (&tbl.hnom),sizeof(tbl.hnom));
        if(bswap)
            tbl.hnom=(float8)bswap_64((int8)tbl.hnom);
    }
    else {
        tbl.hnom=4.;
        cerr<<"can't find coil nominal gain in tbl file "<<tblstr<<endl
            <<"will be set to 4 mV/nT"<<endl;
    }
    tblf.close();

    return tbl;
}

bool positioning(ifstream& tblf, char* label) {
    tblf.seekg(0);

    char labeltmp[5];

    // look for label
    tblf.read(reinterpret_cast <char *> (labeltmp),5);
    while(strcmp(labeltmp,label)!=0 && !tblf.eof()) {
        tblf.seekg(20,ios::cur); // skip to the next label
        tblf.read(labeltmp,5);   // and read it
    }

    // if found label, positioning to get it
    if(strcmp(labeltmp,label)==0) {
        tblf.seekg(7,ios::cur);
        return true;
    }
    else
        return false;
}


double convlat(const char* LATG) {
    string latstr(LATG,13);

    string degstr(latstr,0,2);
    int deg=atoi(degstr.c_str());

    string::size_type coma=latstr.find(",");
    string minstr(latstr,2,coma-2);
    double min=atof(minstr.c_str());

    double signal=1.;
    if(latstr.at(coma+1)=='S')
        signal*=-1.;

    return signal*(deg+min/60.);
}

double convlon(const char* LNGG) {
    string lonstr(LNGG,13);

    string degstr(lonstr,0,3);
    int deg=atoi(degstr.c_str());

    string::size_type coma=lonstr.find(",");
    string minstr(lonstr,3,coma-3);
    double min=atof(minstr.c_str());

    double signal=1.;
    if(lonstr.at(coma+1)=='W')
        signal*=-1.;

    return signal*(deg+min/60.);
}
