#include <cstring>
#include <iostream>
#include <fstream>
//#include <strstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "dcmp2gmt.h"

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc<3 || argc>4) {
        cerr<<argv[0]<<": transforma um arquivo *.dcmp"
                     <<" em um formato para o GMT."<<endl<<endl;
        cerr<<"uso: "<<argv[0]<<" -I<arquivo_entrada> -D<opcao> [-Q]"<<endl
            <<"-I<arquivo_entrada>: arquivo *.dcmp"<<endl
            <<"-D<opcao>          : o tipo de dado que se deseja, as seguites opcoes sao validas:"<<endl
            <<"                     [azim]uth, [sh]ear, [ch]annelling, [tw]ist"<<endl
            <<"                     rhoa, rhob, phasea|phia, phaseb|phib, [err]or, skew, anis, phadif|phidif"<<endl
            <<"-Q                 : forca os azimutes para o primeiro quadrante positivo ou negativo."<<endl
            <<"                     so' tem efeito para as opcoes azimuth e channelling"<<endl;
        exit(1);
    }

    ifstream ifile;
    opcoes opcao;

    bool first=false;

    // le os argumentos da linha de comando
    for(int i=1; i<argc; i++) {
        if(argv[i][0]=='-' && argv[i][1]=='I') {
            ifile.open(argv[i]+2);
            if(!ifile) {
                cerr<<"Nao pode abrir o arquivo \""<<argv[i]+2<<'\"'<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='D') {
            if(strstr(argv[i]+2,"azim")!=0)
                opcao=azim;
            else if(strstr(argv[i]+2,"sh")!=0)
                opcao=sh;
            else if(strstr(argv[i]+2,"ch")!=0)
                opcao=ch;
            else if(strstr(argv[i]+2,"tw")!=0)
                opcao=tw;
            else if(strstr(argv[i]+2,"rhoa")!=0)
                opcao=rhoa;
            else if(strstr(argv[i]+2,"rhob")!=0)
                opcao=rhob;
            else if(strstr(argv[i]+2,"phia")!=0 || strstr(argv[i]+2,"phasea")!=0)
                opcao=phia;
            else if(strstr(argv[i]+2,"phib")!=0 || strstr(argv[i]+2,"phaseb")!=0)
                opcao=phib;
            else if(strstr(argv[i]+2,"phadif")!=0 || strstr(argv[i]+2,"phidif")!=0)
                opcao=phidif;
            else if(strstr(argv[i]+2,"err")!=0)
                opcao=err;
            else if(strstr(argv[i]+2,"skew")!=0)
                opcao=skew;
            else if(strstr(argv[i]+2,"anis")!=0)
                opcao=anis;
            else {
                cout<<"nao reconeceu a opcao "<<argv[i]+2<<" do comando -D"<<endl;
                exit(1);
            }
        }
        else if(argv[i][0]=='-' && argv[i][1]=='Q') {
            first=true;
        }
        else {
            cerr<<"nao reconhece a opcao "<<argv[i]<<endl
	        <<"rode o programa sem opcoes para ver o modo de usar"<<endl;
	    exit(1);
        }
    }

    vector<datum> data=le_dados(ifile,opcao);
    ifile.close();

    if( first && (opcao==azim || opcao ==ch))
        data=tofirst(data);

    cout.setf(ios::scientific);
    cout.precision(5);
    for(int i=0; i<data.size(); i++) {
        cout<<setw(14)<<data[i].T
            <<setw(14)<<data[i].val
            <<setw(14)<<data[i].min
            <<setw(14)<<data[i].val
            <<setw(14)<<data[i].val
            <<setw(14)<<data[i].max<<endl;
    }
}
