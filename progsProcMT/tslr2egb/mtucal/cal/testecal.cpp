#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<iomanip>
// os arquivos de calibracao do sistema fazem a correcao
//  para [T]^-1 e [V/m]^-1 e de calibracao da caixa fazem a correcao
// para [V]^-1

// este programa recupera (em amplitude) a correcao para [nT] e [uV]
// e o fator de calibracao 
// obs: [uV/m]=[mV/km]

// o programa deve ser invocado seguido de duas letras AB, onde:
// A=H ou E indicando se e' um canal magnetico ou eletrico respectivamente
// B=S ou B indicando se e' um calibracao do sistema ou da caixa
// ao usar a opcao E, deve-se seguir o comprimento do dipolo
// ao usar a opcao H, deve-se seguir o ganho relativo

// entrada padrao: frequencia, parte real e imaginaria da calibracao
// saida padrao: frequencia, amplitude da calibracao, correcao Egbert, calfac

const double NS=6.4;
const long FS=0x7FFFFF;
const int HGN=2;
const int EGN=3;

int main (int argc, char* argv[]) {
    double f, re, im;

    if(argv[1][0]=='H') {   // se for canal magnetico
        double gr=atof(argv[2]);
        while(cin>>f>>re>>im) {
            if(argv[1][1]=='S')   // se for cal. sistema
                if(double a=sqrt(re*re+im*im))
                    cout<<setw(14)<<f
                        <<setw(14)<<a
                        <<setw(14)<<1.e9/(a*FS)
                        <<setw(14)<<(a*NS*(250.e-9/gr))/HGN<<endl;
            if(argv[1][1]=='B')   // se for cal. caixa
                if(double a=sqrt(re*re+im*im))
                    cout<<setw(14)<<f
                        <<setw(14)<<a
                        <<setw(14)<<(250./gr)/(a*FS)
                        <<setw(14)<<(a*NS)/HGN<<endl;
        }
    }
    else if(argv[1][0]=='E') {  // se for canal eletrico
        double eln=atof(argv[2]);
        while(cin>>f>>re>>im) {
            if(argv[1][1]=='S')   // se for cal. sistema
                if(double a=sqrt(re*re+im*im))
                    cout<<setw(14)<<f
                        <<setw(14)<<a
                        <<setw(14)<<1.e6*eln/(a*FS)
                        <<setw(14)<<(a*NS)/(EGN*eln)<<endl;
            if(argv[1][1]=='B')   // se for cal. caixa
                if(double a=sqrt(re*re+im*im))
                    cout<<setw(14)<<f
                        <<setw(14)<<a
                        <<setw(14)<<1.e6/(a*FS)
                        <<setw(14)<<(a*NS)/EGN<<endl;
        }
    }
}
