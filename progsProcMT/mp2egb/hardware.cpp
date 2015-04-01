#include "hardware.h"
#include <cctype>

int idchannel (string label);

void readhardware(string calstr, table tbl, vector<calibration> &hardware) {

    // store default values
    // i=0:Hx 1:Hy 2:Hz 3:Ex 4:Ey
    hardware.resize(tbl.nchan);
    for (int i=0; i<tbl.nchan; i++) {
        hardware[i].gain=1.0;
        hardware[i].lpass=0.2;
        if(i<3)
            hardware[i].hpass=0.0;
        else
            hardware[i].hpass=30000.0;
    }


    // try to open the hardware file
    ifstream in(calstr.c_str());
    if(in.fail()) {
        cout<<"can't open the hardware file "<<calstr<<endl
            <<"will use default values"<<endl;
    }
    else {
        cout<<"opened hardware file "<<calstr<<endl;

        // look for lunchbox
        string line;
        while (getline(in,line)) {
            istringstream istrline(line.c_str());
            int snum;
            istrline>>snum;
            if (snum==tbl.snum)
                break;
        }

        if (in.eof()) {
            cout<<"can't find lunchbox "<<tbl.snum<<endl
                <<"will use default values"<<endl;
        }
        else {
            for(int i=0; i<tbl.nchan; i++) {
                getline(in, line);
                istringstream istrline(line.c_str());
                string label;
                istrline>>label;

                int n=idchannel(label);
                // n<0 is a unknow label identification
                // will be ignored
                if(n>=0) {
                    istrline>>hardware[n].gain
                            >>hardware[n].lpass
                            >>hardware[n].hpass;
                }
            }
        }
        in.close();
    }
}

int idchannel (string label) {
    string::const_iterator p=label.begin();
    string ch;

    while(p!=label.end()) {
        if(isalpha(*p))
            ch+=toupper(*p);
        p++;
    }


    if(ch=="HX" || ch=="H")
        return 0;
    else if(ch=="HY" || ch=="D")
        return 1;
    else if(ch=="HZ" || ch=="Z")
        return 2;
    else if(ch=="EX" || ch=="N")
        return 3;
    else if(ch=="EY" || ch=="E")
        return 4;
    else {
        cout<<"unknow label "<<label<<endl
            <<"it will be ignored"<<endl;
        return -1;
    }
}
