#include "modelos.h"

vector< vector<float> > le_mdl_rebocc(string arqE, float xo) {
    ifstream ifile;
    ifile.open(arqE.c_str());
    if(!ifile) {
        cerr<<"Nao pode abrir o arquivo \""<<arqE<<'\"'<<endl;
        exit(1);
    }

    string line;
    while(getline(ifile,line))
        if(line.find("NY")!=string::npos)
            break;
    // numero de colunas
    int c=atoi(line.substr(line.find("NY")+2,string::npos).c_str());

    vector<float> L(c+1);  // vetor com a coordenada do centro da celula
    float Lo=0.;
    for(int i=1; i<=c;i++) {
        float l;
        ifile>>l;
        L[i]=Lo+l/2.;
        if(Lo==xo) {
	  L[0]=L[i]; // guarda o centro da celula em que ocorre a origem
        }
        Lo+=l;
    }
    if(L[0]==0.) {
        cerr<<"a coordenada da origem xo="<<xo<<" deve coincidir com alguma grade horizontal"<<endl;
        exit(1);
    }

    while(getline(ifile,line))
        if(line.find("NZB")!=string::npos)
            break;

    // numero de linhas
    int r=atoi(line.substr(line.find("NZB")+3,string::npos).c_str());

    vector< vector<float> > M(r+1);
    // preenche a 1a. linha com as coordenadas horizontais
    M[0].resize(c+1);
    M[0]=L;

    // preenche a 1a. coluna com a coordenada do centro da celula
    Lo=0.;
    for(int i=1; i<=r;i++) {
        M[i].resize(c+1);
        float l;
        ifile>>l;
        M[i][0]=Lo+l/2.;
        Lo+=l;
    }

    while(getline(ifile,line))
        if(line.find("RESISTIVITY_MODEL")!=string::npos)
            break;

    for(int i=1;i<=r;i++) {
        for(int j=1;j<=c;j++) {
            ifile>>M[i][j];
        }
    }
    ifile.close();
    return M;
}
