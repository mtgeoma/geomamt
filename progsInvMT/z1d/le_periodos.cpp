#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "fwd1d.h"
using namespace std;

vector<double> le_cmdT(string Tstr);

vector<double> le_periodos(string Tstr) {
  vector<double> T;
  ifstream ifile(Tstr.c_str());
  if(ifile.fail()) {
    T=le_cmdT(Tstr);
    if(T.size()==0) {
      cerr<<"nao pode abrir o arquivo "<<Tstr<<endl
          <<"ou formato do comando errado"<<endl;
      exit(1);
    }
    return T;
  }

  string word;
  int w=1;
  while(ifile>>word) {
    double t=atof(word.c_str());
    if(t<=0.) {
      cerr<<"periodo no. "<<w<<" <=0. em "<<Tstr<<endl;
      exit(1);
    }
    T.push_back(t);
    w++;
  }
  return T;
}

vector<double> le_cmdT(string Tstr) {
  vector<double> T;
  int p0=atoi(Tstr.substr(0,Tstr.find("/")).c_str());
  Tstr=Tstr.substr(Tstr.find("/")+1,string::npos);
  int nt=atoi(Tstr.substr(0,Tstr.find("/")).c_str());
  Tstr=Tstr.substr(Tstr.find("/")+1,string::npos);
  int np=atoi(Tstr.substr(0,string::npos).c_str());

  int i=0;
  if(nt>0 && np>0) {
    T.resize(np*nt+1);
    for(int d=p0;d<p0+np;d++)
      for(int t=0;t<nt;t++)
        T[i++]=pow(10.,d+(1.*t)/nt);
  }
  T[i]=pow(10.,p0+np);
  return T;
}
