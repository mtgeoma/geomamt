#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "modelos.h"
#include "geometria.h"
using namespace std;

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc<4||argc>5) {
        cerr<<argv[0]<<": modifica um modelo 2D do Rebooc"
                     <<" por poligonos pre-definidos."<<endl<<endl;
        cerr<<"uso: "<<argv[0]<<" -I<arquivo_entrada> -O<arquivo_saida> -P<arquivo_poligono> [-X<xo>]"<<endl
            <<"-I<arquivo_entrada> : modelo gerado pelo rebooc"<<endl
            <<"-O<arquivo_saida>   : modelo do rebooc modificado"<<endl
            <<"-P<arquivo_poligono>: arquivo com poligonos definidos no estilo do GMT"<<endl
            <<"-X<xo>              : posicao inicial da grade"<<endl;
        exit(1);
    }

    string ifile, ofile, pfile;
    float xo=0.;

    // le os argumentos da linha de comando
    for(int i=1; i<argc; i++) {
        if(argv[i][0]=='-' && argv[i][1]=='I') {
	    ifile=argv[i]+2;
	}
        else if(argv[i][0]=='-' && argv[i][1]=='O') {
            ofile=argv[i]+2;
	}
        else if(argv[i][0]=='-' && argv[i][1]=='P') {
            pfile=argv[i]+2;
	}
        else if(argv[i][0]=='-' && argv[i][1]=='X') {
	    xo=atof(argv[i]+2);
	}
        else {
            cerr<<"nao reconhece a opcao "<<argv[i]<<endl
	        <<"rode o programa sem opcoes para ver o modo de usar"<<endl;
	    exit(1);
        }
    }

    vector< vector<float> > M;
    M=le_mdl_rebocc(ifile,xo);
    vector<poligono> pol;
    pol=le_poligono(pfile,xo);
    M=modelo(M,pol);
    escreve_mdl_rebocc(ifile,ofile,M);
}
