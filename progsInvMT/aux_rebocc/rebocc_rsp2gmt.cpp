#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

enum opcoes {rho, phi, rel, img};
vector<bool> read_dif(string strdif,int nT,opcoes opcao,int est);

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc<4 || argc>5) {
        cerr<<argv[0]<<": le a resposta rho ou phi de uma estacao  do modelo 2D do Rebooc"<<endl<<endl;
        cerr<<"uso: "<<argv[0]<<" -I<arquivo_entrada> -D<dado> -S<no. da estacao> [-M<data_inclusion>]"<<endl
            <<"-I<arquivo_entrada>: modelo gerado pelo rebooc"<<endl
            <<"-D<dado>           : dado a ser selecionado, pode ser"<<endl
            <<"                     rho, phi, rel or img"<<endl
            <<"-S<no. da estacao> : no. da estacao, conforme a ordem do modelo rebocc"<<endl
            <<"-M<data_inclusion> : so plota os dados incluidos no data inclusion file gerado pelo rebocc"<<endl;
        exit(1);
    }

    ifstream ifile;
    string strdif="";
    opcoes opcao;
    int est;

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
            if(strstr(argv[i]+2,"rho")!=0)
	        opcao=rho;
            else if(strstr(argv[i]+2,"phi")!=0)
	        opcao=phi;
            else if(strstr(argv[i]+2,"rel")!=0)
	        opcao=rel;
            else if(strstr(argv[i]+2,"img")!=0)
	        opcao=img;
            else {
                cerr<<"opcao \""<<argv[i]+2<<"\" desconhecida"<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='S') {
            est=atoi(argv[i]+2);
        }
        else if(argv[i][0]=='-' && argv[i][1]=='M') {
            strdif=argv[i]+2;
        }
        else {
            cerr<<"nao reconhece a opcao "<<argv[i]<<endl
	        <<"rode o programa sem opcoes para ver o modo de usar"<<endl;
	    exit(1);
        }
    }

    string line;
    while(getline(ifile,line))
        if(line.find("NUMBER_OF_PERIOD")!=string::npos||line.find("NUMBER_OF_FREQUENCY")!=string::npos)
            break;
    // le periodos
    int nT;
    vector<float> T;
    if(line.find("NUMBER_OF_PERIOD")!=string::npos) {
        nT=atoi(line.substr(line.find("PERIOD")+6,string::npos).c_str());
        //cout<<nT<<endl;
        T.resize(nT);
        for(int t=0;t<nT;t++)
        ifile>>T[t];
    }
    else if(line.find("NUMBER_OF_FREQUENCY")!=string::npos) {
        nT=atoi(line.substr(line.find("FREQUENCY")+9,string::npos).c_str());
        //cout<<nT<<endl;
        T.resize(nT);
        for(int t=0;t<nT;t++) {
            ifile>>T[t];
            T[t]=1./T[t];
        }
    }
    else {
        cerr<<"nao achou no. de periodos ou frequencia"<<endl;
        exit(1);
    }

    while(getline(ifile,line))
        if(line.find("NUMBER_OF_STATION")!=string::npos)
            break;
    int nS=atoi(line.substr(line.find("STATION")+7,string::npos).c_str());
    if(est>nS) {
        cerr<<"no. de estacoes="<<nS<<endl
            <<"logo, nao se pode selecionar a estacao "<<est<<endl;
        exit(1);
    }

    vector<bool> it(nT,true);  // vetor com os periodos a serem plotados
    if(strdif.size()!=0)       // verifica os periodos a serem plotados no data inclusion file
        it=read_dif(strdif,nT,opcao,est);

    cout.setf(ios::fixed);
    cout.precision(4);
    switch(opcao) {
    case rho:
        while(getline(ifile,line))
            if(line.find("app")!=string::npos)
                break;

        if(line.find("app")==string::npos) {
	    cerr<<"nao achou rho (app)"<<endl;
            exit(1);
        }
        else {
	    bool log=false;
            if(line.find("log")!=string::npos)
	        log=true;
	    int si=0;
            if(line.find("*")!=string::npos)
	        si=1;
            for(int t=0;t<nT;t++) {
                for(int s=si;s<=nS;s++) {
       	            double data;
                    ifile>>data;
                    if(s==est && it[t]) {
		        if(log)
			    cout<<setw(12)<<T[t]<<setw(12)<<pow(10.,data)<<endl;
                        else
       	                    cout<<setw(12)<<T[t]<<setw(12)<<data<<endl;
                    }
                }
            }
        }
        break;
    case phi:
        while(getline(ifile,line))
            if(line.find("phs")!=string::npos)
                break;

        if(line.find("phs")==string::npos) {
	    cerr<<"nao achou phi (phs)"<<endl;
            exit(1);
        }
        else {
	    bool rad=false;
            if(line.find("rad")!=string::npos)
	        rad=true;
	    int si=0;
            if(line.find("*")!=string::npos)
	        si=1;
            for(int t=0;t<nT;t++) {
                for(int s=si;s<=nS;s++) {
       	            float data;
                    ifile>>data;
                    if(s==est && it[t]) {
		        if(rad)
			  cout<<setw(12)<<T[t]<<setw(12)<<data*(180./M_PI)<<endl;
                        else
       	                    cout<<setw(12)<<T[t]<<setw(12)<<data<<endl;
                    }
                }
            }
        }
        break;
    case rel:
        while(getline(ifile,line))
            if(line.find("rel")!=string::npos)
                break;

        if(line.find("rel")==string::npos) {
	    cerr<<"nao achou tipper parte real (rel)"<<endl;
            exit(1);
        }
        else {
	    int si=0;
            if(line.find("*")!=string::npos)
	        si=1;
            for(int t=0;t<nT;t++) {
                for(int s=si;s<=nS;s++) {
       	            float data;
                    ifile>>data;
                    if(s==est && it[t])
       	                cout<<setw(12)<<T[t]<<setw(12)<<data<<endl;
                }
            }
        }
        break;
    case img:
        while(getline(ifile,line))
            if(line.find("img")!=string::npos)
                break;

        if(line.find("img")==string::npos) {
	    cerr<<"nao achou tipper parte imaginaria (img)"<<endl;
            exit(1);
        }
        else {
	    int si=0;
            if(line.find("*")!=string::npos)
	        si=1;
            for(int t=0;t<nT;t++) {
                for(int s=si;s<=nS;s++) {
       	            float data;
                    ifile>>data;
                    if(s==est && it[t])
       	                cout<<setw(12)<<T[t]<<setw(12)<<data<<endl;
                }
            }
	}
    }
    ifile.close();
}

vector<bool> read_dif(string strdif,int nT,opcoes opcao,int est) {
    ifstream ifile(strdif.c_str());
    if(!ifile) {
        cerr<<"Nao pode abrir o arquivo \""<<strdif<<'\"'<<endl;
        exit(1);
    }

    string line;
    // procura o inicio da matriz de inclusao de interesse
    if(opcao==rho || opcao==rel) {
        while(getline(ifile,line))
            if(line.find("DATA_INCLUSION_NO_1")!=string::npos)
                break;
        if(line.find("DATA_INCLUSION_NO_1")==string::npos) {
	    cerr<<"nao achou DATA_INCLUSION_NO_1"<<endl;
            exit(1);
        }
    }
    else {
        while(getline(ifile,line))
            if(line.find("DATA_INCLUSION_NO_2")!=string::npos)
                break;
        if(line.find("DATA_INCLUSION_NO_2")==string::npos) {
	    cerr<<"nao achou DATA_INCLUSION_NO_2"<<endl;
            exit(1);
        }
    }
    vector<bool> it(nT);
    for(int t=0;t<nT;t++) {
        getline(ifile,line);
        istringstream sline(line);
        string word;
        sline>>word;
        if(word[est-1]=='0')
            it[t]=false;
        else if(word[est-1]=='1')
            it[t]=true;
        else {
            cerr<<"data inclusion file com valores diferentes de zero ou um"<<endl;
            exit(1);
        }
    }

    return it;
}
