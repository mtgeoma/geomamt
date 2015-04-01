#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc<3 || argc>5) {
        cerr<<argv[0]<<": transforma um modelo 2D do Rebooc"
                     <<" em um formato para o GMT."<<endl<<endl;
        cerr<<"uso: "<<argv[0]<<" -I<arquivo_entrada> -O<arquivo_saida> [-S<escala> -L]"<<endl
            <<"-I<arquivo_entrada>: modelo gerado pelo rebooc"<<endl
            <<"-O<arquivo_saida>  : arquivo para ser lido pelo GMT"<<endl
            <<"                     este arquivo pode ser plotado"<<endl
            <<"                     usando um comando como:"<<endl
            <<"                     psxy arquivo_saida -R -JX -M -H -L -Ccolors.cpt -K -O -V >> fig.ps"<<endl
            <<"-S<escala>         : as distancias serao divididas pela escala"<<endl
            <<"                     por exemplo, -S1.e3 transforma as distancias de metros para quilometros"<<endl
            <<"-L                 : calcula o log da resistividade"<<endl;
        exit(1);
    }

    ifstream ifile;
    ofstream ofile;
    float escala=1.;
    bool logrho=false;

    // le os argumentos da linha de comando
    for(int i=1; i<argc; i++) {
        if(argv[i][0]=='-' && argv[i][1]=='I') {
            ifile.open(argv[i]+2);
            if(!ifile) {
                cerr<<"Nao pode abrir o arquivo \""<<argv[i]+2<<'\"'<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='O') {
            ofile.open(argv[i]+2);
            if(!ofile) {
                cerr<<"Nao pode abrir o arquivo \""<<argv[i]+2<<'\"'<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='S') {
            escala=atof(argv[3]+2);
        }
        else if(argv[i][0]=='-' && argv[i][1]=='L') {
            logrho=true;
        }
        else {
            cerr<<"nao reconhece a opcao "<<argv[i]<<endl
	        <<"rode o programa sem opcoes para ver o modo de usar"<<endl;
	    exit(1);
        }
    }

    string line;
    while(getline(ifile,line))
        if(line.find("NY")!=string::npos)
            break;
    // numero de colunas
    int c=atoi(line.substr(line.find("NY")+2,string::npos).c_str());

    vector<float> L(c+1);  // vetor com os comprimentos acumulativos
    L[0]=0.; // 1o. nó ocorre na origem
    for(int j=1; j<=c;j++) {
        ifile>>L[j];
        L[j]=L[j]/escala+L[j-1];
    }

    while(getline(ifile,line))
        if(line.find("NZB")!=string::npos)
            break;

    // numero de linhas
    int r=atoi(line.substr(line.find("NZB")+3,string::npos).c_str());

    vector<float> H(r+1);  // vetor com as profundidades acumulativas
    H[0]=0.; // 1o. nó ocorre na origem
    for(int i=1; i<=r;i++) {
        ifile>>H[i];
        H[i]=H[i]/escala+H[i-1];
    }

    while(getline(ifile,line))
        if(line.find("RESISTIVITY_MODEL")!=string::npos)
            break;
    vector< vector<double> > M(r);
    float rhomin=1e14,
          rhomax=0.;

    for(int i=0;i<r;i++) {
        M[i].resize(c);
        for(int j=0;j<c;j++) {
            ifile>>M[i][j];
	    if(M[i][j]<rhomin)
	        rhomin=M[i][j];
	    if(M[i][j]>rhomax)
	        rhomax=M[i][j];
        }
    }
    ifile.close();
    if(logrho)
        ofile<<"# rho min= "<<log10(rhomin)<<"; rho max= "<<log10(rhomax)<<" [log10(rho)]"<<endl;
    else
        ofile<<"# rho min= "<<rhomin<<"; rho max= "<<rhomax<<endl;

    for(int i=0;i<r;i++) {
        for(int j=0;j<c;j++) {
            if(logrho)
                ofile<<"> -Z"<<log10(M[i][j])<<endl;
            else
                ofile<<"> -Z"<<M[i][j]<<endl;

            ofile<<L[j]  <<"  "<<H[i]  <<endl
                 <<L[j+1]<<"  "<<H[i]  <<endl
                 <<L[j+1]<<"  "<<H[i+1]<<endl
                 <<L[j]  <<"  "<<H[i+1]<<endl;
        }
    }
    ofile.close();
}
