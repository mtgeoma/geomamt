#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "readtbl.h"

vector<string> getwords(string line);

vector<string> writecal(string calstr, string caldir, string sendir, table &tbl) {
    string line;

    line=caldir+calstr;
    ifstream cal(line.c_str());
    if(cal.fail()) {
        cout<<"can't open the calibration file "<<calstr<<endl;
        exit(1);
    }
    cout<<"opened calibration file "<<calstr<<endl;

    // get header
    getline(cal,line);
    vector<string> words=getwords(line);

    // check if MTU serial number is the same from tbl file
    if(atoi(words[1].c_str())!=tbl.snum) {
        cerr<<"FATAL ERROR: diferent MTU serial number:"<<endl
            <<"in calibration file: "<<words[1]<<endl
            <<"in tbl file        : "<<tbl.snum<<endl;
        exit(1);
    }

    // get the field type parameter
    tbl.fieldtype=atoi(words[2].c_str());
    //cout<<"field type: "<<tbl.fieldtype<<endl;
    if(tbl.fieldtype!=0 && tbl.fieldtype!=1) {
        cerr<<tbl.fieldtype<<" is an unknow field type parameter"<<endl;
        exit(1);
    }

    // create the output files
    vector<string> chstr(5);
    chstr[0]=calstr.substr(0,calstr.find_last_of("."))+".ex";
    chstr[1]=calstr.substr(0,calstr.find_last_of("."))+".ey";
    chstr[2]=calstr.substr(0,calstr.find_last_of("."))+".hx";
    chstr[3]=calstr.substr(0,calstr.find_last_of("."))+".hy";
    chstr[4]=calstr.substr(0,calstr.find_last_of("."))+".hz";

    ofstream ch[5];

    for(int i=0; i<5; i++) {
        string chpath=sendir+chstr[i];
        ch[i].open(chpath.c_str());
        if(ch[i].fail()) {
            cout<<"can't open the output calibration file "<<chpath<<endl;
            exit(1);
        }
    cout<<"opened output calibration file "<<chpath<<endl;
    }

    // separate the calibration file
    int n=1;
    while(getline(cal,line)) {
        words=getwords(line);

        if(words.size()!=12) {
            cerr<<"error reading calibration file "<<calstr<<endl
                <<" at line "<<n<<endl;
            exit(1);
        }

        for(int i=0;i<5;i++)
            ch[i]<<setw(14)<<words[0]            // frequency
                 <<setw(14)<<words[i*2+2]        // real part
                 <<setw(14)<<words[i*2+3]<<endl; // imaginary part
        n++;
    }
    return chstr;
}

// get words separated by commas in line
vector<string> getwords(string line) {
    vector<string> words;
    while(line.find(",")!=string::npos) {
        string word=line.substr(0,line.find(","));
        words.push_back(word);
        line=line.substr(line.find(",")+1,string::npos);
    }
    return words;
}
