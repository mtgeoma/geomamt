#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<iomanip>

// faz a transformacao: f re im => f amp fase

int main () {
    double f, re, im;
    while(cin>>f>>re>>im)
        cout<<setw(14)<<f
            <<setw(14)<<sqrt(re*re+im*im)
            <<setw(14)<<atan2(im,re)*(180/M_PI)<<endl;
}
