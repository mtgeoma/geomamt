#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc !=5) {
        string tab(strlen(argv[0]),' ');
        cerr<<argv[0]<<": read a ascii SN file generated by multmtrn"<<endl
                <<tab<<"  and output the required channel for the required station."<<endl<<endl;
        cerr<<"usege: "<<argv[0]<<" -I<input_file> -D<data> -S<no. station> -C<no. channel>\n"
            <<"-I<input_file> : input SN ascii file\n"
            <<"-D<data>       : choose the data to the output, data could be:"<<endl
            <<"                 egv => Eigenvalues in noise units"<<endl
            <<"                 sgp => SIGNAL POWER"<<endl
            <<"                 inp => INCOHERENT NOISE POWER"<<endl
            <<"                 rnp => RELATIVE INCOHERENT NOISE POWER"<<endl
            <<"-S<no. station>: choose the number station"<<endl
            <<"-C<no. channel>: choose the number channel"<<endl
            <<"obs: no. station is the same order in the array.cfg file (1st station=>1)"<<endl
            <<"     no. channel is the same of the time series file"<<endl
            <<"                 normally Hx=1, Hy=2, Hz=3, Ex=4, Ey=5"<<endl;
        exit(1);
    }

    enum data {egv,sgp,inp,rnp};
    vector<string> datalabel(4);
    // labels as defined in multmtrn.f
    datalabel[egv]="Eigenvalues in noise units";
    datalabel[sgp]="SIGNAL POWER";
    datalabel[inp]="INCOHERENT NOISE POWER";
    datalabel[rnp]="RELATIVE INCOHERENT NOISE POWER";

    ifstream ifile;
    data datum;
    int channel, station;
    // read commands
    for(int i=1; i<argc; i++) {
        if(argv[i][0]=='-' && argv[i][1]=='I') {
            ifile.open(argv[i]+2, ios::in | ios::binary);
            if(!ifile) {
                cerr<<"Couldn'd open the file \""<<argv[i]+2<<'\"'<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='D') {
            if(strcmp(argv[i]+2,"egv")==0)
                datum=egv;
            else if(strcmp(argv[i]+2,"sgp")==0)
                datum=sgp;
            else if(strcmp(argv[i]+2,"inp")==0)
                datum=inp;
            else if(strcmp(argv[i]+2,"rnp")==0)
                datum=rnp;
            else {
                cerr<<'\"'<<argv[i]+2<<"\" is a unknow data output."<<endl
                    <<"run the program whithout option to see the possible data output"<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='S')
            station=atoi(argv[i]+2);
        else if(argv[i][0]=='-' && argv[i][1]=='C')
            channel=atoi(argv[i]+2);
        else {
            cerr<<argv[i]<<" is an unknow parameter"<<endl;
	    exit(1);
        }
    }

    // read no. of channels and no. of periods from first line
    int nC, nT;
    ifile>>nC>>nT;

    // look for the data label
    string line;
    while(getline(ifile,line)) {
        if(line.find(datalabel[datum])!=string::npos)
            break;
    }
    if(line.find(datalabel[datum])==string::npos) {
        cerr<<"couldn't find the data label \""<<datalabel[datum]<<'\"'<<endl;
        exit(1);
    }

    for(int T=0;T<nT;T++) {
        string d;
        // get and write the period
        ifile>>d;
        cout<<" "<<d;
        for(int C=1;C<=nC;C++)  {
            ifile>>d;
            if(C==(station-1)*5+channel)
                cout<<" "<<d<<endl;
        }
    }
    ifile.close();
}
