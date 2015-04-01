#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

  /* Faz a inversao Niblett-Bostick e Weidelt-Bostick a partir da */
  /* entrada padrao no formato: T rho phi                         */

void print_result(double h, double rhon, double rhob);

int main ()
{
    const double uo=4.E-7*M_PI;

    double Ta, rhoa, phia; // ponto anterior
    double T, rho, phi;    // ponto sendo calculado
    double Tp, rhop, phip; // ponto posterior

    // le os dois primeiros dados
    cin>>T>>rho>>phi;
    cin>>Tp>>rhop>>phip;

    // converte para radianos
    phi*=(M_PI/180.);
    phip*=(M_PI/180.);

    // Calcula em T usando a derivada do lado positivo
    double h=sqrt(rho*T/(2.*M_PI*uo));
    double m=(T/rho)*(rhop-rho)/(Tp-T);
    double rhon=abs(rho*(1.+m)/(1.-m));
    double rhob=abs(rho*(M_PI/(2.*phi)-1.));
    print_result(h,rhon,rhob);

    // atualiza as variaveis
    Ta=T;
    rhoa=rho;
    phia=phi;
    T=Tp;
    rho=rhop;
    phi=phip;

    while(cin>>Tp>>rhop>>phip)
    {
        phip*=(M_PI/180.);

        // Calcula em T usando a media das derivadas dos dois lados
        h=sqrt(rho*T/(2.*M_PI*uo));
        m=(T/rho)*(rhop-rhoa)/(Tp-Ta);
        rhon=abs(rho*(1.+m)/(1.-m));
        rhob=abs(rho*(M_PI/(2.*phi)-1.));
        print_result(h,rhon,rhob);

        // atualiza as variaveis
        Ta=T;
        rhoa=rho;
        phia=phi;
        T=Tp;
        rho=rhop;
        phi=phip;
    }

    // Calcula em T usando a derivada do lado negativo
    h=sqrt(rho*T/(2.*M_PI*uo));
    m=(T/rho)*(rho-rhoa)/(T-Ta);
    rhon=abs(rho*(1.+m)/(1.-m));
    rhob=abs(rho*(M_PI/(2.*phi)-1.));
    print_result(h,rhon,rhob);
}

void print_result(double h, double rhon, double rhob) {
    cout<<setw(12)<<setprecision(4)<<h
        <<setw(12)<<setprecision(4)<<rhon
        <<setw(12)<<setprecision(4)<<rhob
        <<endl;
}
