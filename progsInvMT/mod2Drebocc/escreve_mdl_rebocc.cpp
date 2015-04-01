#include "modelos.h"

void escreve_mdl_rebocc(string arqE,string arqS, vector< vector<float> > M) {
    ifstream ifile;
    ifile.open(arqE.c_str());
    if(!ifile) {
        cerr<<"Nao pode abrir o arquivo \""<<arqE<<'\"'<<endl;
        exit(1);
    }
    ofstream ofile;
    ofile.open(arqS.c_str());
    if(!ifile) {
        cerr<<"Nao pode abrir o arquivo \""<<arqS<<'\"'<<endl;
        exit(1);
    }

    string line;
    while(getline(ifile,line))
        if(line.find("RESISTIVITY_MODEL")!=string::npos)
            break;
        else
            ofile<<line<<endl;

    ofile<<"RESISTIVITY_MODEL"<<endl;
    ofile.setf(ios::scientific);
    ofile.precision(4);

    for(int i=1;i<M.size();i++) {
        for(int j=1;j<M[0].size();j++) {
            ofile<<"  "<<M[i][j];
        }
        ofile<<endl;
    }

    ifile.close();
    ofile.close();
}
