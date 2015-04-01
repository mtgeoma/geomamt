#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

const string dsep="/";

vector<double> getconv(string spstr, double error, double &cda, double &cdb);
void writets(string instr, string outstr, vector<double> conv, double cda, double cdb);

int main(int argc, char* argv[])
{
    string instr="", outstr="", spstr="", datadir="."+dsep, spdir="."+dsep;
    double error=0.;

    if (argc<3 || argc>4 ) {  // wrong number of parameters
        cerr<<"usage: "<<argv[0]<<" -I<egbert file> -O<ts file> [-E<error value>]"<<endl
            <<"input : egbert file (ascii)"<<endl
            <<"output: time series with time mark (julian days starting at the begining of year)"<<endl
            <<"error value: only for AP and RI"<<endl
            <<"this program will look for paths.cfg file in the current directory"<<endl
            <<"if didn't find, will assume that all files are at the current directory"<<endl;
        exit (1);
    }
    else {  // get the parameters
        for(int i=1; i< argc; i++) {
            if(argv[i][0]=='-' && argv[i][1]=='I')
                instr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='O')
                outstr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='E')
                error=atof(argv[i]+2);
            else {
                cerr<<argv[i]<<" is an unknow parameter"<<endl;
                exit(1);
            }
        }
    }

    // try to open the path.cfg
    string path="."+dsep+"paths.cfg";
    ifstream pcfg(path.c_str());
    if(pcfg.fail()) {
        cout<<"can't open the "<<path<<" file"<<endl
        <<"will assume that all files are at the current directory"<<endl;
    }
    else {  // get the data and system parameter directory
        string line;
        getline(pcfg,line);
        istringstream istrline1(line.c_str());
        istrline1>>line;
        datadir=datadir+line+dsep;

        getline(pcfg,line);
        istringstream istrline2(line.c_str());
        istrline2>>line;
        spdir=spdir+line+dsep;
        pcfg.close();
    }

    // write the sp file name
    spstr=outstr.substr(0,outstr.find_last_of("."))+".sp";

    // get conversion factors from sp file
    double cda, cdb;
    vector<double> conv;
    conv=getconv(spdir+spstr,error,cda,cdb);

    // write time serie file
    writets(datadir+instr,datadir+outstr,conv,cda,cdb);
}
