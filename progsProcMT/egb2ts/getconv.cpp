#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <complex>
#include <sstream>
#include <vector>

using namespace std;

typedef struct {
  double f;
  complex<double> cmp;
} element;

double filtergain(string filterstr, double error, int type);

vector<double> getconv(string spstr, double error, double &cda, double &cdb) {
    // open sp file
    ifstream sp(spstr.c_str());
    if(sp.fail()) {
        cerr<<"can't open "<<spstr<<" file"<<endl;
        exit(1);
    }

    string line;
    for(int i=0; i<3; i++) // skip first 3 lines
        if(!getline(sp,line)) {
            cerr<<"error reading first 3 lines in "<<spstr<<" file"<<endl;
            exit(1);
        }

    short nchan;
    if(!getline(sp,line)) {
        cerr<<"error reading number of channels in "<<spstr<<" file"<<endl;
        exit(1);
    }
    else {
        istringstream iline(line.c_str());
        iline>>nchan;
    }

    double sr;
    if(!getline(sp,line)) {
        cerr<<"error reading sample rate in "<<spstr<<" file"<<endl;
        exit(1);
    }
    else {
        istringstream iline(line.c_str());
        iline>>sr;
    }

    // clock drift corrections; assume clock time = tc; actual time = t
    // these are related by tc = t + (cda + cdb*t); assumes t in seconds
    if(!getline(sp,line)) {
        cerr<<"error reading clock drift corrections in "<<spstr<<" file"<<endl;
        exit(1);
    }
    else {
        istringstream iline(line.c_str());
        iline>>cda>>cdb;
    }

    // read channels parameters
    vector<double> conv(nchan);
    for(int i=0; i<nchan; i++) {
        if(!getline(sp,line)) {
            cerr<<"error reading channel id "<<i+1
                <<"in "<<spstr<<" file"<<endl;
            exit(1);
        }

        if(line[0]=='H') { // magnetic channel
            getline(sp,line);  // skip sensor orientation
            double scale;
            double filter=1.;
            int nfilter;
            if(!getline(sp,line)) {
                cerr<<"error reading count conversion in "<<spstr<<" file"<<endl;
                exit(1);
            }
            else {
                istringstream iline(line.c_str());
                iline>>scale>>nfilter;

                for(int j=0; j<nfilter; j++) {
                    getline(sp,line);
                    if(error!=0. && (line.find("RI"))) {
                        string filterstr;
                        getline(sp,line);
                        istringstream iline(line.c_str());
                        iline>>filterstr;
                        filterstr="./sensors/"+filterstr;
                        filter*=filtergain(filterstr,error,0);
                    }
                    else if(error!=0. && (line.find("AP"))) {
                        string filterstr;
                        getline(sp,line);
                        istringstream iline(line.c_str());
                        iline>>filterstr;
                        filterstr="./sensors/"+filterstr;
                        filter*=filtergain(filterstr,error,1);
                    }
                    else
                        getline(sp,line);
                }
            }
            conv[i]=scale*filter;
        }
        else if(line[0]=='E') { // electric channel
            double length, angle, tilt, gain;
            if(!getline(sp,line)) {
                cerr<<"error reading electric channel"<<endl;
                exit(1);
            }
            else {
                istringstream iline(line.c_str());
                iline>>length>>angle>>tilt>>gain;
            }

            double scale;
            double filter=1.;
            int nfilter;
            if(!getline(sp,line)) {
                cerr<<"error reading count conversion in "<<spstr<<" file"<<endl;
                exit(1);
            }
            else {
                istringstream iline(line.c_str());
                iline>>scale>>nfilter;

                for(int j=0; j<nfilter; j++) {
                    getline(sp,line);
                    if(error!=0. && (line.find("RI"))) {
                        string filterstr;
                        getline(sp,line);
                        istringstream iline(line.c_str());
                        iline>>filterstr;
                        filterstr="./sensors/"+filterstr;
                        filter*=filtergain(filterstr,error,0);
                    }
                    else if(error!=0. && (line.find("AP"))) {
                        string filterstr;
                        getline(sp,line);
                        istringstream iline(line.c_str());
                        iline>>filterstr;
                        filterstr="./sensors/"+filterstr;
                        filter*=filtergain(filterstr,error,1);
                    }
                    else
                        getline(sp,line);
                }
            }
            conv[i]=scale*filter/(length*gain);
        }
    }

    sp.close();
    return conv;
}

double filtergain(string filterstr, double error, int type) {
    // open filter file
    ifstream filter(filterstr.c_str());
    if(filter.fail()) {
        cerr<<"can't open "<<filterstr<<" file"<<endl;
        exit(1);
    }

    vector<element> function;
    string line, wr[3];
    while(getline(filter,line)) {
        istringstream iline(line.c_str());
        iline>>wr[0]>>wr[1]>>wr[2];
        if(isdigit(wr[0][0])) {
            element el;
            el.f=atof(wr[0].c_str());
            double x=atof(wr[1].c_str());
            double y=atof(wr[2].c_str());

            if(type==0)
                el.cmp=complex<double>(x,y);
            else if(type==1)
                el.cmp=polar<double>(x,y*(M_PI/180.));
            else {
                cerr<<"unknow type"<<endl;
                exit(1);
            }
            function.push_back(el);
        }
    }
    filter.close();
    
    int b=0,
        e=function.size(),
        m=(e+b)/2;

    double gain=1.;
    for(;;) { // forever
        int n;
        double sx=0., sx2=0.;
        double mean1, sdev1;
        for(int i=b;i<m;i++) {
            double x=abs(function[i].cmp);
            sx+=x;
            sx2+=x*x;
        }
        n=m-b+1;
        mean1=sx/n;
        sdev1=sqrt((sx2-(sx*sx)/n)/(n-1));

        sx=0.;
        sx2=0.;
        double mean2, sdev2;
        for(int i=m;i<e;i++) {
            double x=abs(function[i].cmp);
            sx+=x;
            sx2+=x*x;
        }
        n=e-m;
        mean2=sx/n;
        sdev2=sqrt((sx2-(sx*sx)/n)/(n-1));
        if(sdev1<sdev2) {
            b=b;
            e=m;
            m=(e+b)/2;
            if(sdev1/mean1<error || b>=m || m>=e) {
                cout<<b<<":"<<e<<"="<<mean1<<endl;
                return 1./mean1;
            }
        }
        else {
            b=m;
            e=e;
            m=(e+b)/2;
            if(sdev2/mean2<error || b>=m || m>=e) {
                cout<<b<<":"<<e<<"="<<mean2<<endl;
                return 1./mean2;
            }
        }
    }
}
