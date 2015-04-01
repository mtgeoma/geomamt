#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "fwd1d.h"
using namespace std;

vector<camada> le_modelo(string Mstr) {
  ifstream ifile(Mstr.c_str());
  if(ifile.fail()) {
    cerr<<"nao pode abrir o arquivo "<<Mstr<<endl;
    exit(1);
  }

  vector<camada> C;
  string line;
  int l=1;
  while(getline(ifile,line)) {
    double rho=0., h=0.;
    istringstream istr(line);
    string word;
    int w=1;
    while(istr>>word) {
      if(w==1)
        rho=atof(word.c_str());
      else if(w==2)
        h=atof(word.c_str());
      else {
        cerr<<"mais do que 2 colunas na linha "<<l<<" em "<<Mstr<<endl;
        exit(1);
      }
//      cout<<rho<<":"<<h<<endl;
      w++;
    }
    if(rho<=0.) {
      cerr<<"rho<=0 na linha "<<l<<" em "<<Mstr<<endl;
      exit(1);
    }
    if(h<0.) {
      cerr<<"h<0 na linha "<<l<<" em "<<Mstr<<endl;
      exit(1);
    }
    camada c;
    c.rho=rho;
    c.h=h;
    C.push_back(c);

    if(h==0) {
//        cout<<w<<endl;
      if(--w==1)
        break;
      else {
        cerr<<"h==0 na linha "<<l<<" em "<<Mstr<<endl;
        exit(1);
      }
    }
    l++;
  }
  if(getline(ifile,line)) {
    cerr<<"camada sem espessura na linha "<<l<<" antes do final de "<<Mstr<<endl;
    exit(1);
  }
  if(C[C.size()-1].h!=0.) {
      cerr<<"a ultima camada em "<<Mstr<<" deve ter apenas rho definido"<<endl;
    exit(1);
  }
  return C;
}
