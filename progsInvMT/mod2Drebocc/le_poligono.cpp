#include "geometria.h"

// decompoe uma linha em um vetor de palavras
vector<string> le_linha(const string line);

vector<poligono> le_poligono(string arqP, float xo) {
    ifstream ifile;
    ifile.open(arqP.c_str());
    if(!ifile) {
        cerr<<"Nao pode abrir o arquivo \""<<arqP<<'\"'<<endl;
        exit(1);
    }

    string line;
    int np=-1;
    vector<poligono> Vp;
    while(getline(ifile,line)) {
        vector<string> L=le_linha(line);
        if(L[0]==">") {
	    np++;
            poligono p;
            Vp.push_back(p);
            float rho=atof(L[1].c_str());
            if(rho<=0.) {
	        cerr<<"resistividade menor ou igual a zero no "<<np+1<<"-esimo segmento."<<endl;
                exit(1);
            }
            Vp[np].rho=rho;
        }
        else {
	    if(np==-1) {
                cerr<<"arquivo "<<arqP<<" nao inicia com >."<<endl;
                exit(1);
            }
            ponto P;
            P.x=atof(L[0].c_str())+xo;
            P.y=atof(L[1].c_str());
            Vp[np].P.push_back(P);
        }
    }

    // verifica se tem mais que 3 pontos e fecha o poligono, se for o caso.
    for(int p=0;p<Vp.size();p++) {
        int n=Vp[p].P.size();

        if(Vp[p].P[0].x==Vp[p].P[n-1].x &&Vp[p].P[0].y==Vp[p].P[n-1].y) {
	    if(Vp[p].P.size()<4) {
	        cerr<<p+1<<"-esimo poligono fechado com menos de 4 pontos."<<endl;
                exit(1);
            }
        }
        else {
	    if(Vp[p].P.size()<3) {
	        cerr<<p+1<<"-esimo poligono aberto com menos de 3 pontos."<<endl;
                exit(1);
            }
            Vp[p].P.push_back(Vp[p].P[0]);
        }
    }
    return Vp;
}

// decompoe uma linha em um vetor de palavras
vector<string> le_linha(const string line) {
    istringstream istr(line.c_str());
    string word;
    vector<string> words;
    while(istr>>word)
        words.push_back(word);

    return words;
}
