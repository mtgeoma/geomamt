#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "dcmp2j.h"

int main(int argc, char *argv[])
{
    // verifica se o no. de argumentos e' correto, se nao, informa o modo de usar
    if(argc<3 || argc>4) {
        cerr<<argv[0]<<": transforma um arquivo *.dcmp"
                     <<" no formato Jones."<<endl<<endl;
        cerr<<"uso: "<<argv[0]<<" -I<arquivo_entrada> -O<arquivo_saida> [-A<anisotropia>]"<<endl
            <<"-I<arquivo_entrada>: arquivo *.dcmp"<<endl
            <<"-O<arquivo_saida>  : arquivo formato Jones"<<endl
            <<"-A<anisotropia>    : valor para correcao da anisotropia"<<endl
            <<"                     usando apenas -A, sera feito o calculo interativo como o dcmp2j do strike"<<endl;
        exit(1);
    }

    string   ifstr;
    string   ofstr;
    ifstream ifile;
    ofstream ofile;
    double anis=-1.; // anis<0 => nao corrige anisotropia
                     // anis=0 => calculo interativo da correcao da anisotropia
                     // anis>0 => usa o valor de entrada para correcao da anisotropia

    // le os argumentos da linha de comando
    for(int i=1; i<argc; i++) {
        if(argv[i][0]=='-' && argv[i][1]=='I') {
	    ifstr=argv[i]+2;
            ifile.open(argv[i]+2);
            if(!ifile) {
                cerr<<"Nao pode abrir o arquivo \""<<argv[i]+2<<'\"'<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='O') {
	    ofstr=argv[i]+2;
            ofile.open(argv[i]+2);
            if(!ofile) {
                cerr<<"Nao pode abrir o arquivo \""<<argv[i]+2<<'\"'<<endl;
                exit(1);
            }
	}
        else if(argv[i][0]=='-' && argv[i][1]=='A') {
	    string temp(argv[i]+2);
            if(temp.size()==0)
	        anis=0.;
            else {
	        anis=atof(temp.c_str());
                if(anis<=0.) {
		    cerr<<"valor da correcao da anisitropia deve ser maior que zero"<<endl;
                    exit(1);
                }
            }
        }
        else {
            cerr<<"nao reconhece a opcao "<<argv[i]<<endl
	        <<"rode o programa sem opcoes para ver o modo de usar"<<endl;
	    exit(1);
        }
    }

    vector<string> cab=le_cab(ifile);
    vector<datum> azim=le_dados(ifile,AZIM);
    vector<datum> sh=le_dados(ifile,SH);
    vector<datum> ch=le_dados(ifile,CH);
    vector<datum> tw=le_dados(ifile,TW);
    vector<datum> rhoa=le_dados(ifile,RHOA);
    vector<datum> rhob=le_dados(ifile,RHOB);
    vector<datum> phia=le_dados(ifile,PHIA);
    vector<datum> phib=le_dados(ifile,PHIB);
    vector<datum> err=le_dados(ifile,ERR);
    ifile.close();

    ofile.setf(ios::fixed);

    // correct anisotropy
    if( anis==0.) {
        ofile.precision(3);
        cout<<"  #           Period        Anisotropy"<<endl;
        for(int i=0;i<rhoa.size();i++) {
            cout<<setw(5)<<i+1<<setw(15)<<rhoa[i].T<<setw(15)<<sqrt(rhoa[i].val/rhob[i].val)<<endl;
	}

        int n_an1=0, n_an2;
        do {
            cout<<"Give period range to average (e.g. 1 5)"<<endl;
            cin>>n_an1>>n_an2;
            if(n_an1<=0 || n_an1>n_an2 || n_an2 > rhoa.size() ) {
	        cout<<"Inappropriate values entered: re-enter"<<endl;
                n_an1=0;
            }
        } while(n_an1==0);
        n_an1--;
        n_an2--;

        for(int i=n_an1;i<=n_an2;i++)
            anis+=sqrt(rhoa[i].val/rhob[i].val);

        anis/=(float(n_an2-n_an1+1));
        cout<<"Average anisotropy to be applied: "<<anis<<endl;
    }

    if( anis>0.) {
        for(int i=0;i<rhoa.size();i++) {
	    rhoa[i].val/=anis;
            rhoa[i].min/=anis;
            rhoa[i].max/=anis;
	    rhob[i].val*=anis;
            rhob[i].min*=anis;
            rhob[i].max*=anis;
        }
    }


    // calcula as medias
    double Mazim=0.;
    for(int i=0;i<azim.size();i++)
        Mazim+=azim[i].val;
    Mazim/=azim.size();

    double Msh=0.;
    for(int i=0;i<sh.size();i++)
        Msh+=sh[i].val;
    Msh/=sh.size();

    double Mch=0.;
    for(int i=0;i<ch.size();i++)
        Mch+=ch[i].val;
    Mch/=ch.size();

    double Mtw=0.;
    for(int i=0;i<tw.size();i++)
        Mtw+=tw[i].val;
    Mtw/=tw.size();

    double Merr=0.;
    for(int i=0;i<err.size();i++)
        Merr+=err[i].val;
    Merr/=err.size();

    // procura azimute original
    string azim0="oops";
    for(int i=0;i<cab.size();i++) {
        if(cab[i].find(">AZIMUTH")!=string::npos) {
	    azim0=cab[i].substr(cab[i].find("=")+1,string::npos);
            break;
        }
    }

    // escreve cabecalho
    ofile.precision(7);
    ofile<<"# Written by dcmp2j: input file >"<<ifstr<<endl;
    ofile<<"# date: "<<endl;
    ofile<<"#"<<endl;
    ofile<<"#"<<endl;
    ofile<<"# azimuth of original data:"<<azim0<<endl;
    ofile<<"#"<<endl;
    ofile<<"# The azimuth listed under AZIMUTH is the average value"<<endl;
    ofile<<"#"<<endl;
    ofile<<"# average distortion parameters:"<<endl;
    ofile<<"# azimuth:"<<setw(15)<<Mazim<<endl;
    ofile<<"# shear  :"<<setw(15)<<Msh<<endl;
    ofile<<"# twist  :"<<setw(15)<<Mtw<<endl;
    ofile<<"# chann  :"<<setw(15)<<Mch<<endl;
    ofile<<"# error  :"<<setw(15)<<Merr<<endl;
    ofile<<"#"<<endl;
    ofile<<"#     Period       Azim    Shear   Twist   Chann     Error"<<endl;
    for(int i=0;i<azim.size();i++) {
        ofile<<"#"<<setw(15)<<azim[i].T;
        ofile.precision(1);
        ofile<<setw(8)<<azim[i].val;
        ofile<<setw(8)<<sh[i].val;
        ofile<<setw(8)<<tw[i].val;
        ofile<<setw(8)<<ch[i].val;
        ofile.precision(7);
        ofile<<setw(15)<<err[i].val<<endl;
    }
    ofile<<"#"<<endl;

    // acrescenta cabecalho original
    for(int i=0;i<cab.size();i++) {
        if(cab[i].find(">AZIMUTH")!=string::npos) {
            ofile.precision(4);
            ofile<<cab[i].substr(0,cab[i].find("=")+1);
            ofile<<setw(12)<<Mazim<<endl;
        }
        else
            ofile<<cab[i]<<endl;
    }
    ofile.precision(2);
    ofile<<ofstr
         <<"  AZIM"<<setw(8)<<Mazim
         <<"  SHEAR"<<setw(8)<<Msh
         <<"  CHANN"<<setw(8)<<Mch
         <<"  TWIST"<<setw(8)<<Mtw<<endl;

    // escreve os dados
    //ofile.setf(ios::fmtflags(0));
    //ofile.setf(ios::scientific);
    //ofile.precision(6);
    ofile<<"RXY"<<endl;
    ofile<<rhoa.size()<<endl;
    for(int i=0;i<rhoa.size();i++) {
        ofile.setf(ios::scientific);
        ofile.precision(6);
        ofile<<setw(14)<<rhoa[i].T<<setw(14)<<rhoa[i].val<<setw(14)<<phia[i].val;
        if(rhoa[i].max==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<rhoa[i].max;
        if(rhoa[i].min==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<rhoa[i].min;


        if(phia[i].max==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<phia[i].max;
        if(phia[i].min==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<phia[i].min;
        ofile<<"  1.00  1.00"<<endl;
    }

    ofile<<"RYX"<<endl;
    ofile<<rhob.size()<<endl;
    for(int i=0;i<rhob.size();i++) {
        ofile<<setw(14)<<rhob[i].T
             <<setw(14)<<rhob[i].val
             <<setw(14)<<phib[i].val-180.;
        if(rhob[i].max==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<rhob[i].max;
        if(rhob[i].min==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<rhob[i].min;


        if(phib[i].max==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<phib[i].max-180.;
        if(phib[i].min==0.)
            ofile<<setw(14)<<"-999.";
        else
            ofile<<setw(14)<<phib[i].min-180.;
        ofile<<"  1.00  1.00"<<endl;
    }

    ofile.close();
}
