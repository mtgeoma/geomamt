#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include "fwd1d.h"

using namespace std;

vector<camada> le_modelo(string Mstr);
vector<double> le_periodos(string Tstr);

int main(int argc, char *argv[]) {
  // verifica no. de argumentos:
  if(!(argc==3 || argc==4)) {
    cerr<<"mt1d: calcula rho e phi para um modelo 1D"<<endl
        <<"uso: "<<argv[0]<<" -M<modelo> -T<periodos> [digitos]"<<endl
        <<"<modelo>  : arquivo com os dados do modelo 1D que se deseja calcular"<<endl
        <<"    cada linha deve ter o valor de rho e espessura de uma camada "<<endl
        <<"    sendo que a ultima linha, correspondente ao semi-espaco infinito,"<<endl
        <<"    deve ter apenas o valor de rho"<<endl<<endl
        <<"<periodos>: arquivo com os periodos alvos"<<endl
        <<"    um periodo por linha"<<endl
        <<"    alternativamente, pode-se usar o comando na forma"<<endl
        <<"    -Tpo/nt/nd sendo po a potencia de 10 do periodo inicial"<<endl
        <<"       nt o no. de periodos por decada e"<<endl
        <<"       nd o no. de decadas para se calcular os periodos"<<endl
	<<"<digitos>:no. casas decimais [4]"<<endl;
        exit(1);
  }

  string Mstr, Tstr;
  int ndig=4;
  // le os argumentos
  for(int i=1;i<argc;i++) {
    string cmd(argv[i]);
    if(cmd.substr(0,2)=="-M")
      Mstr=cmd.substr(2,string::npos);
    else if(cmd.substr(0,2)=="-T")
      Tstr=cmd.substr(2,string::npos);
    else if(ndig=boost::lexical_cast<int>(cmd));
    else {
      cout<<cmd<<": opcao nao conhecida"<<endl;
      exit(1);
    }
  }

  vector<camada> C=le_modelo(Mstr);
  vector<double> T=le_periodos(Tstr);

  cout.setf(ios::scientific);
  cout.precision(ndig);
 	for(int i=0; i<T.size(); i++) {
    complex<double> Z=z1d(T[i],C);
    cout<<setw(10+ndig)<<T[i]
        <<setw(10+ndig)<<norm(Z)*T[i]/(2.*M_PI*u0)
        <<setw(10+ndig)<<arg(Z)*(180./M_PI)
        <<setw(10+ndig)<<real(Z)
        <<setw(10+ndig)<<imag(Z)
        <<endl;
  }
}
