#include<iostream>
#include<fstream>
#include<strstream>
#include<string>
#include<cmath>
#include <iomanip>

string le_dados(string line);

int main(int argc, char *argv[]) {
    if(argc != 3) {
        cerr<<argv[0]<<": transforma o arquivo *.cts de Re. Im. para Amp. fase"<<endl
            <<"Uso: "<<argv[0]<<" <arquivo de entrada> <arquivo de saida>"<<endl;
        exit(1);
    }
    
    ifstream in(argv[1]);
    if(!in) {
        cerr<<"Nao pode abrir o arquivo \""<<argv[1]<<'\"'<<endl;
        exit(1);
    }

    ofstream out(argv[2]);
    if(!out) {
        cerr<<"Nao pode abrir o arquivo \""<<argv[2]<<'\"'<<endl;
        exit(1);
    }
    
    string line;
    // pula a 1o. linha
    //getline(in,line);
    
    while(getline(in,line)) {
        string saida=le_dados(line);
        out<<saida;
    }
}

string le_dados(string line) {
    istrstream istr(line.c_str());
    double freq, tmp, ch[5][2];
    
    istr>>freq>>tmp;
    for(int i=0; i<5; i++)
        for(int j=0; j<2; j++)
            istr>>ch[i][j];
    
    ostrstream ostr;
    
    ostr.precision(4);
    ostr<<setw(12)<<freq;
    
    ostr.setf(ios::scientific,ios::floatfield);
    
    for(int i=0; i<5; i++) {
        double amp=sqrt(pow(ch[i][0],2)+pow(ch[i][1],2));
        double fase=atan2(ch[i][1],ch[i][0])*(180./M_PI);
        ostr<<"  "<<setw(12)<<amp<<"  "<<setw(12)<<fase;
    }
    
    ostr<<endl;        
    return ostr.str();
}
