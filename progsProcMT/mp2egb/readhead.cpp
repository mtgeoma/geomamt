#include <vector>
#include <byteswap.h>
#include "readhead.h"
#include <iomanip>

table readhead(string instr) {
    table tbl;
    if(sizeof(float4)!=4) {
      cout<<"some type haven't the correct size:\n"
	  <<"size of float4 "<<setw(4)<<sizeof(float4)<<"\n"
	  <<"you will need to edit the readhead.h file."
	  <<"size of float      "<<setw(4)<<sizeof(float)<<"\n"
	  <<"size of double     "<<setw(4)<<sizeof(double)<<"\n"
	  <<"size of long double"<<setw(4)<<sizeof(long double)<<"\n";
      exit(1);
    }
    if(bswap)
        cout<<"WORDS_BIGENDIAN"<<endl;
    else
        cout<<"WORDS_LOWENDIAN"<<endl;
    // open file
    ifstream in(instr.c_str(), ios::in | ios::binary);
    if(in.fail()) {
        cout<<"can't open the LRMT file "<<instr<<endl;
        exit(1);
    }
    cout<<"opened LRMT file "<<instr<<" to read the head"<<endl;

    // project code [Survey ID]
    char prjcod[4];
    in.read(reinterpret_cast <char *> (prjcod),3);
    prjcod[3]='\0';
    tbl.srvy=prjcod;

    // run number [Survey ID+Site ID of the measurement]
    char runnum[8];
    in.read(reinterpret_cast <char *> (runnum),7);
    runnum[7]='\0';
    tbl.site=runnum;
    tbl.site=tbl.site.substr(3,4);  // remove Survey ID

    // window number [Current data acquisition window]
    in.read(reinterpret_cast <char *> (&tbl.awin),sizeof(tbl.awin));
    if(bswap)
        tbl.awin=bswap_16(tbl.awin);

    // skip site description
    in.seekg(20,ios::cur);

    // Ex line length (metres) [Ex dipole length, m]
    in.read(reinterpret_cast <char *> (&tbl.exln),sizeof(tbl.exln));
    if(bswap)
        tbl.exln=bswap_16(tbl.exln);

    // Ey line length (metres) [Ey dipole length, m]
    in.read(reinterpret_cast <char *> (&tbl.eyln),sizeof(tbl.eyln));
    if(bswap)
        tbl.eyln=bswap_16(tbl.eyln);

    // latitude (degrees) (minutes) (seconds)
    int16_t coord[3];
    in.read(reinterpret_cast <char *> (coord),sizeof(coord));
    if(bswap)
        for(short i=0; i<3; i++)
            coord[i]=bswap_16(coord[i]);

    tbl.latg=coord[0]/1.+coord[1]/60.+coord[2]/3600.;

    // latitude direction (N or S)
    char dir;
    in.read(reinterpret_cast <char *> (&dir),sizeof(dir));
    if(dir=='S' || dir=='s' )
        tbl.latg*=-1.;

    // longitude (degrees) (minutes) (seconds)
    in.read(reinterpret_cast <char *> (coord),sizeof(coord));
    if(bswap)
        for(short i=0; i<3; i++)
            coord[i]=bswap_16(coord[i]);

    tbl.lngg=coord[0]/1.+coord[1]/60.+coord[2]/3600.;

    // longitude direction (E or W)
    in.read(reinterpret_cast <char *> (&dir),sizeof(dir));
    if(dir=='W' || dir=='w' )
        tbl.lngg*=-1.;

    // hardware telluric gain [Gain for "E" channels]
    in.read(reinterpret_cast <char *> (&tbl.itgain),sizeof(tbl.itgain));
    if(bswap)
        tbl.itgain=bswap_16(tbl.itgain);

    // magnetometer factor
    tbl.magfac=1.E4;

    // skip telluric filter settling time (min)
    in.seekg(2,ios::cur);

    // Ex azimuth (degrees) [Ex sensor azimuth]
    in.read(reinterpret_cast <char *> (&tbl.exaz),sizeof(tbl.exaz));
    if(bswap)
        tbl.exaz=bswap_16(tbl.exaz);

    // Ey azimuth (degrees) [Ey sensor azimuth]
    in.read(reinterpret_cast <char *> (&tbl.eyaz),sizeof(tbl.eyaz));
    if(bswap)
        tbl.eyaz=bswap_16(tbl.eyaz);

    // Hx azimuth (degrees) [Hx sensor azimuth]
    in.read(reinterpret_cast <char *> (&tbl.hxaz),sizeof(tbl.hxaz));
    if(bswap)
        tbl.hxaz=bswap_16(tbl.hxaz);

    // az direction (from true or magntc north)
    bool magdir;
    char azdir[17];
    in.read(reinterpret_cast <char *> (azdir),17);

    if(strstr(azdir,"MAGNETIC")!=0)
        magdir=true;
    else
        magdir=false;

    // serial number of lunchbox [Serial Number of LRMT unit]
    in.read(reinterpret_cast <char *> (&tbl.snum),sizeof(tbl.snum));
    if(bswap)
        tbl.snum=bswap_16(tbl.snum);

    // skip serial number of magnetometer and serial number of tellurics
    in.seekg(4,ios::cur);

    // start time day
    in.read(reinterpret_cast <char *> (&tbl.stime.day),sizeof(tbl.stime.day));
    if(bswap)
        tbl.stime.day=bswap_16(tbl.stime.day);

    // start time month
    in.read(reinterpret_cast <char *> (&tbl.stime.month),sizeof(tbl.stime.month));
    if(bswap)
        tbl.stime.month=bswap_16(tbl.stime.month);

    // start time year
    in.read(reinterpret_cast <char *> (&tbl.stime.year),sizeof(tbl.stime.year));
    if(bswap)
        tbl.stime.year=bswap_16(tbl.stime.year);
    if(tbl.stime.year<100) // in the form 99 for 1999
        tbl.stime.year+=1900;

    // start time hour
    in.read(reinterpret_cast <char *> (&tbl.stime.hour),sizeof(tbl.stime.hour));
    if(bswap)
        tbl.stime.hour=bswap_16(tbl.stime.hour);

    // start time minute
    in.read(reinterpret_cast <char *> (&tbl.stime.min),sizeof(tbl.stime.min));
    if(bswap)
        tbl.stime.min=bswap_16(tbl.stime.min);


    string model_file="/usr/local/share/IGRF12.COF";
    tbl.decl=declination(tbl.latg, tbl.lngg, 0., tbl.stime.year, tbl.stime.month, tbl.stime.day, model_file);

    if(!magdir) {           // if not in magnetic north (hopefully in true north)
        tbl.exaz+=tbl.decl; // put azimuts in relation to the magnetic north
        tbl.eyaz+=tbl.decl;
        tbl.hxaz+=tbl.decl;
    }

    // skip duration time
    in.seekg(6,ios::cur);

    // magnetometer factor=25600 microvolts/bit (only for baseline)
    int16_t magunit;
    in.read(reinterpret_cast <char *> (&magunit),sizeof(magunit));
    if(bswap)
        magunit=bswap_16(magunit);

    // A to D factor = 305 microvolts/bit
    in.read(reinterpret_cast <char *> (&tbl.fs),sizeof(tbl.fs));
    if(bswap)
        tbl.fs=bswap_16(tbl.fs);

    // number of channels in GROUP1
    int16_t chG1;
    in.read(reinterpret_cast <char *> (&chG1),sizeof(chG1));
    if(bswap)
        chG1=bswap_16(chG1);

    // number of channels in GROUP2
    int16_t chG2;
    in.read(reinterpret_cast <char *> (&chG2),sizeof(chG2));
    if(bswap)
        chG2=bswap_16(chG2);

    if(chG1!=5 || chG2 !=0) {
        cout<<"Sorry, actually this program work only for:"<<endl
            <<"number of channels in GROUP1 = 5; but it is ="<<chG1<<endl
            <<"number of channels in GROUP2 = 0; but it is ="<<chG2<<endl;
        exit(1);
    }
    tbl.nchan=chG1;

    // skip filter type of GROUP1
    in.seekg(2,ios::cur);

    // sample rate
    float4 sr;
    in.read(reinterpret_cast <char *> (&sr),sizeof(sr));
    if(bswap) {
        union {
            int32_t   i;
            float4 f;
        } sr2;

        sr2.f=sr;
        sr2.i=bswap_32(sr2.i);
        sr=sr2.f;
    }

    char srunit[4];
    in.read(reinterpret_cast <char *> (srunit),4);

    if(strstr(srunit,"sec")!=0)  // if sample rate in seconds
        sr=1./sr;                // change to Hz

    tbl.srpm=int32_t(sr*60.);  // sample rate per minute

    // get gains. Channels order must be 1 Hx 2 Hy 3 Hz 4 Ex 5 Ey
    vector<string> chlabels;
    chlabels.push_back("Hx");
    chlabels.push_back("Hy");
    chlabels.push_back("Hz");
    chlabels.push_back("Ex");
    chlabels.push_back("Ey");
    // char* chlabels[5]={"Hx","Hy","Hz","Ex","Ey"};

    for(int i=0;i<5;i++) {
        int16_t chnum;
        in.read(reinterpret_cast <char *> (&chnum),sizeof(chnum));
        if(bswap)
            chnum=bswap_16(chnum);
        
        char chlabel[4];
        in.read(reinterpret_cast <char *> (chlabel),4);
        in.read(reinterpret_cast <char *> (&tbl.gain[i]),sizeof(tbl.gain[i]));
        if(bswap)
            tbl.gain[i]=bswap_16(tbl.gain[i]);
           

        if(chnum!=i+1 && strstr(chlabel,chlabels[i].c_str())==0) {
            cerr<<"error in channel number "<<chnum<<endl;
            exit(1);
        }
    }

    // baseline
    in.read(reinterpret_cast <char *> (&tbl.bsex),sizeof(tbl.bsex));
    in.read(reinterpret_cast <char *> (&tbl.bsey),sizeof(tbl.bsey));
    in.read(reinterpret_cast <char *> (&tbl.bsez),sizeof(tbl.bsez));
    tbl.sizehdr=in.tellg();

    // size of file
    struct stat results;
    if(stat(instr.c_str(),&results)) {
        cerr<<"can't check the size of file "<<instr<<endl;
        exit(1);
    }
    tbl.sizefile=results.st_size;
    in.close();

    return tbl;
}
