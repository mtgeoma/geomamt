#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "read_mackie.h"
using namespace std;

vector<double> get_data(string line);

int main(int argc, char* argv[]) {
    string istr;
    string hstr;
    string sstr="";
    
    if (argc<3 || argc>4) {  // wrong number of parameters
        cerr<<"usage: "<<argv[0]<<" -H<header data> -I<input data> [-S<static_shift>]"<<endl
            <<"<header data> : 1st line: Mode Type (tm or te)"<<endl
            <<"                2nd line: station location"<<endl
            <<"<input data>  : Mackie input data file (ascii)"<<endl
            <<"<static_shift>: optionally, could output the Mackie static shift factors to <static_shift> file"<<endl
            <<"OBS: output for standard output"<<endl;
        exit (1);
    }
    else {  // get the parameters
        for(int i=1; i< argc; i++) {
            if(argv[i][0]=='-' && argv[i][1]=='I')
                istr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='H')
                hstr=argv[i]+2;
            else if(argv[i][0]=='-' && argv[i][1]=='S')
                sstr=argv[i]+2;
            else {
                cerr<<argv[i]<<" is an unknow parameter"<<endl;
                exit(1);
            }
        }
    }

    // open header file
    ifstream hfile(hstr.c_str());
    if(hfile.fail()) {
        cout<<"can't open the "<<hstr<<" file"<<endl;
        exit(1);
    }

    // read Mode Type
    string mode_type;
    getline(hfile,mode_type);

    // read station location
    string line;
    getline(hfile,line);
    vector<double> station_loc=get_data(line);
    getline(hfile,line);
    vector<double> T=get_data(line);

    // open input file
    ifstream ifile(istr.c_str());
    if(ifile.fail()) {
        cout<<"can't open the "<<istr<<" file"<<endl;
        exit(1);
    }

    // read number of stations
    int ns;
    getline(ifile,line);
    istringstream iline(line);
    iline>>ns;
    if(ns!=station_loc.size()) {
      cerr<<"number of station ("<<ns<<") not iqual number of station location ("<<station_loc.size()<<") at "<<hstr<<endl;
      exit(1);
    }

    vector< vector<data> > st(ns);
    vector<float> s_shift(ns);

    for(int s=0;s<ns;s++) {
        st[s].resize(T.size());
        // read number of periods and static shift factor
        int nT;
        float static_shift;
        getline(ifile,line);
        istringstream iline(line);
        iline>>nT>>static_shift;
        if(nT>T.size()) {
	  cout<<"number of periods in station "<<s<<" is bigger than number of periods in "<<hstr<<endl;
          exit(1);
        }
        s_shift[s]=static_shift;
        
        // read data
        int i=0;
        for(int t=0;t<nT;t++) {
            getline(ifile,line);
            data datum=read_line_mackie(line);
            while(i<(T.size()-1) && abs(datum.T-T[i])>abs(datum.T-T[i+1])) {
                st[s][i].T=0.;
                i++;
            }
            if(i==(T.size()-1) && t!=(nT-1)) {
		cerr<<"periods in station "<<s<<" couldn't match periods in "<<hstr<<endl;
                exit(1);
            }
            else {
              datum.T=T[i];
              st[s][i]=datum;
              i++;
            }
        }
        for(;i<T.size();i++) {
	  st[s][i].T=0.;
        }
    }
    ifile.close();

    // write static_shift data
    if(sstr.size()!=0){
        // open header file
        ofstream sfile(sstr.c_str());
        if(sfile.fail()) {
            cout<<"can't open the "<<sstr<<" file"<<endl;
            exit(1);
        }
        sfile<<"DISTORTION_INDEX"<<endl;
        for(int s=0;s<ns;s++)
	    sfile<<"1";
        sfile<<endl<<"DISTORTION_PARAMETER"<<endl;
        for(int s=0;s<ns;s++)
	    sfile<<"  "<<log10(s_shift[s]);
        sfile<<endl<<"DISTORTION_INCLUSION index"<<endl;
        for(int t=0;t<T.size();t++) {
            for(int s=0;s<st.size();s++) {
                if(st[s][t].T!=0.)
                    sfile<<"1";
                else
	        sfile<<"0";
            }
            sfile<<endl;
        }
    }

    // write rebooc data
    cout<<"TITLE                   *"<<endl;
    cout<<"MODE_TYPE               "<<mode_type<<endl;
    cout<<"NUMBER_OF_RESPONSE       2"<<endl;
    cout<<"NUMBER_OF_PERIOD"<<setw(10)<<T.size()<<endl;
    cout.setf(ios::scientific);
    cout.precision(4);
    for(int t=0;t<T.size();t++)
        cout<<setw(12)<<T[t];
    cout<<endl;
    cout<<"NUMBER_OF_STATION"<<setw(9)<<st.size()<<endl;
    for(int s=0;s<station_loc.size();s++)
        cout<<setw(12)<<station_loc[s]*1.e3;
    cout<<endl;

    cout<<"DATA_RESPONSE_NO_1* applog"<<endl;
    for(int t=0;t<T.size();t++) {
        for(int s=0;s<st.size();s++) {
	  if(st[s][t].T!=0.)
            cout<<setw(12)<<log10(st[s][t].rho);
          else
            cout<<setw(12)<<0.;
        }
        cout<<endl;
    }

    cout<<"ERROR_RESPONSE_NO_1* 0.1"<<endl;
    for(int t=0;t<T.size();t++) {
        for(int s=0;s<st.size();s++) {
	  if(st[s][t].T!=0.)
	    cout<<setw(12)<<st[s][t].err*2.*log10(exp(1.));
          else
            cout<<setw(12)<<10.*2.*log10(exp(1.));
        }
        cout<<endl;
    }

    cout<<"DATA_RESPONSE_NO_2* phsrad"<<endl;
    for(int t=0;t<T.size();t++) {
        for(int s=0;s<st.size();s++) {
	  if(st[s][t].T!=0.)
	    cout<<setw(12)<<(M_PI/180.)*(st[s][t].phi);
          else
	    cout<<setw(12)<<0.;
        }
        cout<<endl;
    }

    cout<<"ERROR_RESPONSE_NO_2* 0.05"<<endl;
    for(int t=0;t<T.size();t++) {
        for(int s=0;s<st.size();s++) {
	  if(st[s][t].T!=0.)
	    cout<<setw(12)<<st[s][t].err;
          else
	    cout<<setw(12)<<10.;
        }
        cout<<endl;
    }

    cout<<"DATA_INCLUSION_NO_1 index"<<endl;
    for(int t=0;t<T.size();t++) {
        for(int s=0;s<st.size();s++) {
	  if(st[s][t].T!=0.)
	    cout<<"1";
          else
	    cout<<"0";
        }
        cout<<endl;
    }

    cout<<"DATA_INCLUSION_NO_2 index"<<endl;
    for(int t=0;t<T.size();t++) {
        for(int s=0;s<st.size();s++) {
	  if(st[s][t].T!=0.)
	    cout<<"1";
          else
	    cout<<"0";
        }
        cout<<endl;
    }
}


vector<double> get_data(string line) {
    istringstream sline(line);
    vector<double> data;
    double datum;
    while(sline>>datum) {
      data.push_back(datum);
    }
    return data;
}
