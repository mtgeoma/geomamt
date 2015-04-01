#include <cstdlib>
#include <iostream>
#include <sstream>
#include "read_mackie.h"

data read_line_mackie(string line, string mode_type) {
    data datatmp;

    string word=line.substr(0,line.find("("));
    datatmp.T=atof(word.c_str());

    word=line.substr(line.find("(")+1,line.find(",")-line.find("(")-1);
    datatmp.rho=atof(word.c_str());

    word=line.substr(line.find(",")+1,line.find(")")-line.find(",")-1);
    datatmp.phi=atof(word.c_str());

    line=line.substr(line.find(")")+1,string::npos);
    istringstream isstream (line);
    double rhoerr, phierr;
    if(mode_type=="tp")
        isstream>>phierr;
    else
        isstream>>rhoerr>>phierr;
    datatmp.err=phierr;

    return datatmp;
}
