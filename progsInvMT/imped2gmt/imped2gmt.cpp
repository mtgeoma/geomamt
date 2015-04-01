#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>

using namespace std;

// decompoe uma linha em um vetor de palavras
vector<string> le_linha(const string line);

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc<3 || argc>4) {
        cerr<<argv[0]<<": transforma o arquivo de impedancia do strike"
                     <<" em um formato para o GMT."<<endl<<endl;
        cerr<<"uso: "<<argv[0]<<" -I<arquivo_entrada> -D<opcao> [-S]"<<endl
            <<"-I<arquivo_entrada>: arquivo *.dcmp"<<endl
            <<"-D<opcao>          : o tipo de dado que se deseja, as seguites opcoes sao validas:"<<endl
            <<"                     dados medidos: xx xy yx yy ou dados estimados: exx exy eyx eyy"<<endl
            <<"-S                 : escalona as impedancias pelo fator 1/sqrt(w*mu)"<<endl
            <<"                     so' tem efeito se as impedancias nao estiverem escalonadas"<<endl;
        exit(1);
    }

    ifstream ifile;
    bool escala=false;
    bool medido=true;
    enum {xx,xy,yx,yy} comp;

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
            int skip=2;
            if(argv[i][2]=='e') {
                medido=false;
                skip=3;
            }
            if(strstr(argv[i]+skip,"xx")!=0)
                comp=xx;
            else if(strstr(argv[i]+skip,"xy")!=0)
                comp=xy;
            else if(strstr(argv[i]+skip,"yx")!=0)
                comp=yx;
            else if(strstr(argv[i]+skip,"yy")!=0)
                comp=yy;
            else {
                cout<<"nao reconeceu a opcao "<<argv[i]+2<<" do comando -D"<<endl;
                exit(1);
            }
        }
        else if(argv[i][0]=='-' && argv[i][1]=='S') {
            escala=true;
        }
        else {
            cerr<<"nao reconhece a opcao "<<argv[i]<<endl
	        <<"rode o programa sem opcoes para ver o modo de usar"<<endl;
	    exit(1);
        }
    }

    string line;

    // pula 1o. linha de comentario
    getline(ifile,line);

    while(getline(ifile,line)) {
        // le o periodo e o fator de escala do tensor
        vector<string> dados=le_linha(line);
        if(dados.size()!=2) {
            cerr<<"erro lendo linha de periodos \""<<line<<'\"'<<endl;
            exit(1);
        }
        double T=atof(dados[0].c_str());
        double fator=atof(dados[1].c_str());

        vector< complex<double> > Zm(4);
        vector< complex<double> > Ze(4);
        vector<double> std(4);

        // le a 1o. linha de dados do tensor de impedancia
        getline(ifile,line);
        dados=le_linha(line);
        if(dados.size()==8) {
            Zm[xx]=complex<double>(atof(dados[0].c_str()),atof(dados[1].c_str()));
            Zm[xy]=complex<double>(atof(dados[2].c_str()),atof(dados[3].c_str()));
            Ze[xx]=complex<double>(atof(dados[4].c_str()),atof(dados[5].c_str()));
            Ze[xy]=complex<double>(atof(dados[6].c_str()),atof(dados[7].c_str()));
        }
        else if(dados.size()==4) {
            Zm[xx]=complex<double>(atof(dados[0].c_str()),atof(dados[1].c_str()));
            Zm[xy]=complex<double>(atof(dados[2].c_str()),atof(dados[3].c_str()));

            getline(ifile,line);
            dados=le_linha(line);
            if(dados.size()!=4) {
                cerr<<"erro lendo 2o. linha do tensor de impedancia \""<<line<<'\"'<<endl;
                exit(1);
            }
            Ze[xx]=complex<double>(atof(dados[0].c_str()),atof(dados[1].c_str()));
            Ze[xy]=complex<double>(atof(dados[2].c_str()),atof(dados[3].c_str()));
        }
        else{
            cerr<<"erro lendo 1o. linha do tensor de impedancia \""<<line<<'\"'<<endl;
            exit(1);
        }

        // le a 2o. linha de dados do tensor de impedancia
        getline(ifile,line);
        dados=le_linha(line);
        if(dados.size()==8) {
            Zm[yx]=complex<double>(atof(dados[0].c_str()),atof(dados[1].c_str()));
            Zm[yy]=complex<double>(atof(dados[2].c_str()),atof(dados[3].c_str()));
            Ze[yx]=complex<double>(atof(dados[4].c_str()),atof(dados[5].c_str()));
            Ze[yy]=complex<double>(atof(dados[6].c_str()),atof(dados[7].c_str()));
        }
        else if(dados.size()==4) {
            Zm[yx]=complex<double>(atof(dados[0].c_str()),atof(dados[1].c_str()));
            Zm[yy]=complex<double>(atof(dados[2].c_str()),atof(dados[3].c_str()));

            getline(ifile,line);
            dados=le_linha(line);
            if(dados.size()!=4) {
                cerr<<"erro lendo 4o. linha do tensor de impedancia \""<<line<<'\"'<<endl;
                exit(1);
            }
            Ze[yx]=complex<double>(atof(dados[0].c_str()),atof(dados[1].c_str()));
            Ze[yy]=complex<double>(atof(dados[2].c_str()),atof(dados[3].c_str()));
        }
        else {
            cerr<<"erro lendo 2o. linha do tensor de impedancia \""<<line<<'\"'<<endl;
            exit(1);
        }

        // le a 1o. linha do desvio padrao
        getline(ifile,line);
        dados=le_linha(line);
        if(dados.size()!=2) {
            cerr<<"erro lendo 1o. linha do desvio padrao \""<<line<<'\"'<<endl;
            exit(1);
        }
        std[xx]=atof(dados[0].c_str());
        std[xy]=atof(dados[1].c_str());

        // le a 2o. linha do desvio padrao
        getline(ifile,line);
        dados=le_linha(line);
        if(dados.size()!=2) {
            cerr<<"erro lendo 2o. linha do desvio padrao \""<<line<<'\"'<<endl;
            exit(1);
        }
        std[yx]=atof(dados[0].c_str());
        std[yy]=atof(dados[1].c_str());

        // acerta o fator de escala
        if(escala && fator==1.0)
            fator=sqrt(T/(2.*M_PI*4.*M_PI*1.E-7));
        else
            fator=1.0;

        // escreve os dados
        cout.setf(ios::scientific);
        cout.precision(5);

        if(medido) {
            cout<<setw(14)<<T
                <<setw(14)<<Zm[comp].real()*fator
                <<setw(14)<<Zm[comp].imag()*fator
                <<setw(14)<<std[comp]*fator<<endl;
        }
        else {
            cout<<setw(14)<<T
                <<setw(14)<<Ze[comp].real()*fator
                <<setw(14)<<Ze[comp].imag()*fator<<endl;
        }
    }
    ifile.close();
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
