#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <dirent.h>
#include "bytes.h"
#include "count.h"
#include "readtbl.h"

using namespace std;

const string dsep="/";

void writesp(string spstr, table tbl, vector<string> chstr);
void writeegb(string instr, string outstr,timecount clock);
vector<string> writecal(string calstr, string caldir, string sendir, table &tbl);
bool exist(string target, char type);

int main(int argc, char* argv[])
{
    // verify if the sizes of types are correct
    if(!correct_types())
        exit(1);
    if(bswap)
        cout<<"WORDS_BIGENDIAN"<<endl;
    else
        cout<<"WORDS_LOWENDIAN"<<endl;

    string instr="", outstr="", calstr="", indir="", outdir="";
    string tblstr="", spstr="";
    string caldir="", datdir="", spdir="", sendir;
    string path;

    if (argc<4 || argc>5) {  // wrong number of parameters
        cerr<<"usage: "<<argv[0]<<" -I<MTU-LR file> -C<cal file> -O<egbert file>"
            <<" [-D<output directory>]"<<endl
            <<"MTU-LR file: input binary data from MTU-LR"<<endl
            <<"cal file   : system response calculated by SYSCAL"<<endl
            <<"egbert file: output Egbert format"<<endl<<endl
            <<"this program try to be smart..."<<endl
            <<"it will look for .tbl file in the same directory of the imput file"<<endl
            <<"(if didn't find, will ask for you where it is)"<<endl
            <<"if you choose the option -D, it will look for the file paths.cfg"<<endl
            <<"and will try to use the Egbert's directory struct"<<endl
            <<"(if didn't find, will put all the output files in this directory)"<<endl;
        exit (1);
    }
    else {  // get the parameters
        for(int i=1; i< argc; i++) {
            if(argv[i][0]=='-' && argv[i][1]=='I')
                instr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='C')
                calstr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='O')
                outstr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='D')
                outdir=argv[i]+2;
            else {
                cerr<<argv[i]<<" is an unknow parameter"<<endl;
                exit(1);
            }
        }
    }

    // saparate directory from name of input file
    if(instr.find_last_of(dsep)!=string::npos) {
        indir=instr.substr(0,instr.find_last_of(dsep)+1);
        instr=instr.substr(instr.find_last_of(dsep)+1,string::npos);
    }
    else
        indir="."+dsep;

    // check if MTU-LR file exist
    path=indir+instr;
    if(!exist(path,'f')) {
        cerr<<"can't find the input MTU-LR file "<<path<<endl;
        exit(1);
    }

    // check if tbl file exist
    char awin=instr[instr.length()-1];
    path=path.substr(0,path.find_last_of("."));  // remove the extention
    path=path+".tb"+awin;

    if(!exist(path,'f')) {  // if fail the first attempt
        cout<<"can't find the tbl file "<<path<<endl;

        path[path.length()-1]='l';

        if(!exist(path,'f')) { // prompt the user for the tbl file
            cout<<"can't find the tbl file "<<path<<endl;
            for(;;) { // forever
                cout<<"enter with the tbl file "
                    <<"(just press enter to exit the program)"<<endl;
                getline(cin,path);

                if(path.length()==0)
                    exit(1);

                if(!exist(path,'f'))
                    cout<<"can't find the tbl file "<<path<<endl;
                else
                    break;
            }
        }
    }
    tblstr=path;
    cout<<"found the tbl file "<<tblstr<<endl;

    // saparate directory from name of calibration file
    if(calstr.find_last_of(dsep)!=string::npos) {
        caldir=calstr.substr(0,calstr.find_last_of(dsep)+1);
        calstr=calstr.substr(calstr.find_last_of(dsep)+1,string::npos);
    }
    else
        caldir="."+dsep;

    // try to find the cal file
    path=caldir+calstr;
    if(!exist(path,'f')) {
        cerr<<"can't find the calibration file "<<path<<endl;
        exit(1);
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
        path=outdir+"paths.cfg";
        ifstream pcfg(path.c_str());
        if(pcfg.fail()) {
            cout<<"can't open the "<<path<<" file"<<endl
                <<"all output files will be placed at "<<outdir<<" directory"<<endl;

            if(datdir.length())
                cout<<"the output directory of option -O will be overlapped "
                    <<"by directory of option -D"<<endl;

            datdir=outdir;
            spdir=outdir;
            sendir=outdir;
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

            sendir=outdir+"sensors/";

            if (!exist(sendir,'d'))  // if can't open the sensors directoty
                sendir=caldir;       // will put the cal files at the calibration directory
        }
        pcfg.close();
    }

    // write the sp file path
    spstr=spdir+outstr.substr(0,outstr.find_last_of("."));
    spstr+=".sp";

    // read tbl file
    table tbl=readtbl(tblstr);

    // write the calibration files
    vector<string> chstr=writecal(calstr, caldir, sendir, tbl);

    // start count
    timecount clock=startcount(indir+instr, tbl);

    // write the sp file
    writesp(spstr, tbl, chstr);

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
