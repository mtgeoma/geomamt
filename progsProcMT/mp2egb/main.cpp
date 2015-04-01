#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <dirent.h>
#include "count.h"
#include "hardware.h"

const string dsep="/";

void writesp(string spstr, table tbl, vector<calibration> hardware);
void writeegb(string instr, string outstr, timecount clock);

timecount startcount(table &tbl, bool drop);

bool exist(string target, char type);

int main(int argc, char* argv[])
{
    string instr="", outstr="", calstr="", indir="", outdir="";
    string spstr="";
    string caldir="", datdir="", spdir="";
    bool drop=false;

    if (argc<3 || argc>6) {  // wrong number of parameters
        cerr<<"usage: "<<argv[0]<<" -I<LRMT file> -O<egbert file>"
            <<" [-D<output directory>] [-P<project directory] [-drop]"<<endl
            <<"LRMT file: input binary data from LRMT (.1mp)"<<endl
            <<"egbert file: output Egbert format"<<endl<<endl
            <<"this program try to be smart..."<<endl
            <<"It is assumed that this program is being run from"<<endl
            <<"the database directory of the project"<<endl
            <<"(e.g. /local/mt/lims/nod)"<<endl
            <<"and that all data files reside in subdirectories"<<endl
            <<"with the respective station names"<<endl
            <<"(e.g. /local/mt/lims/nod/nod001 nod002 etc.)"<<endl<<endl
            <<"optionaly, you could define the project directory with the option -P"<<endl<<endl
            <<"if you choose the option -D, it will look for the file paths.cfg"<<endl
            <<"and will try to use the Egbert's directory struct"<<endl
            <<"(if didn't find, will put all the output files in this directory)"<<endl
            <<"-drop will drop points to account for filter settling on data from LiMS"<<endl
            <<"For sampling rate of 30 s, will also add 25 s to start time"<<endl;
        exit (1);
    }
    else {  // get the parameters
        for(int i=1; i< argc; i++) {
            if(argv[i][0]=='-' && argv[i][1]=='I')
                instr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='O')
                outstr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='D')
                outdir=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='P')
                indir=argv[i]+2;
            else if(strstr(argv[i],"-drop")!=0)
                drop=true;
            else {
                cerr<<argv[i]<<" is an unknow parameter"<<endl;
                exit(1);
            }
        }
    }

    // get the directory from name of input file (stn001 in stn001a1.1mp)
    if(indir.length()!=0)  {               // if user enter with the project directory
        if(indir[indir.length()-1]!=dsep[0])  // check if it finish with dsep
            indir+=dsep;
    }
    else
        indir="."+dsep;

    // calibration directory is the same of project directory (indir)
    caldir=indir;

    indir=indir+instr.substr(0,6)+dsep;
    // include extention in the input file name
    if(instr.find(".1mp")==string::npos)
        instr+=".1mp";

    // check if LRMT file exist
    if(!exist(indir+instr,'f')) {
        cerr<<"can't find the input LRMT file "<<indir<<instr<<endl;
        exit(1);
    }

    // look for calibration file (hardware.prj)
    calstr="hardware."+instr.substr(0,3);

    // try to find the cal file
    if(!exist(caldir+calstr,'f')) {
        cout<<"can't find the calibration file "<<caldir<<calstr<<endl
            <<"Enter with calibration file (and its path)"<<endl
            <<"or press ENTER to get the defaults values:"<<endl;
        getline(cin,calstr);
        caldir="";
    }

    // saparate directory from name of egbert file
    if(outstr.find_last_of(dsep)!=string::npos) {
        datdir=outstr.substr(0,outstr.find_last_of(dsep)+1);
        outstr=outstr.substr(outstr.find_last_of(dsep)+1,string::npos);
    }
    else
        datdir="."+dsep;

    // check if option -D was used and if the directory exist
    if(outdir.length()!=0) {
        if (!exist(outdir,'d')) {
            cerr<<"can't open directory "<<outdir<<endl
                <<"permition danied?"<<endl
                <<"remember that '~' don't work as the home directory"<<endl;
                exit(1);
        }
        // check if outdir finish with dsep
        if(outdir[outdir.length()-1]!=dsep[0])
            outdir+=dsep;

        // try to open the path.cfg
        string path=outdir+"paths.cfg";
        ifstream pcfg(path.c_str());
        if(pcfg.fail()) {
            cout<<"can't open the "<<path<<" file"<<endl
                <<"all output files will be placed at "<<outdir<<" directory"<<endl;

            if(datdir.length())
                cout<<"the output directory of option -O will be overlapped "
                    <<"by directory of option -D"<<endl;

            datdir=outdir;
            spdir=outdir;
        }
        else {  // get the data and system parameter directory
            string line;
            getline(pcfg,line);
            istringstream istrline1(line.c_str());
            istrline1>>line;
            datdir=outdir+line+dsep;

            getline(pcfg,line);
            istringstream istrline2(line.c_str());
            istrline2>>line;
            spdir=outdir+line+dsep;
        }
        pcfg.close();
    }

    // write the sp file path
    spstr=spdir+outstr.substr(0,outstr.find_last_of("."));
    spstr+=".sp";

    // read head of LRMT file
    table tbl=readhead(indir+instr);

    // write output
    cout << "baseline " << tbl.bsex*2.56 << " " << tbl.bsey*2.56 << " " << tbl.bsez*2.56 << endl;
    // read calibration file (hardware.prj)
    vector<calibration> hardware;
    readhardware(caldir+calstr, tbl, hardware);

    // start count
    timecount clock=startcount(tbl, drop);

    // write the sp file
    writesp(spstr, tbl, hardware);

    // write the Egbert format
    writeegb(indir+instr,datdir+outstr,clock);
}

// exit look for target of type directory 'd' or file 'f'
// return true if find or false if don't
bool exist(string target, char type) {
    if(type=='d') {
        DIR *dp=opendir(target.c_str());
        if(dp==0)
            return false;
        else {
            (void) closedir (dp);
            return true;
        }
    }
    else if(type=='f') {
        string dir, file;
        // saparate directory from name file
        if(target.find_last_of(dsep)!=string::npos) {
            dir=target.substr(0,target.find_last_of(dsep)+1);
            file=target.substr(target.find_last_of(dsep)+1,string::npos);
        }
        else {
            dir="."+dsep;
            file=target;
        }

        DIR *dp=opendir(dir.c_str());
        if(dp==0)
            return false;
        else {
            struct dirent *ep;
            while (ep = readdir (dp))
                if(file.compare(ep->d_name)==0)
                    return true;

            (void) closedir (dp);
            return false;
        }
    }
    else {
        cerr<<type<<": unknow type"<<endl;
        exit(1);
    }
}
