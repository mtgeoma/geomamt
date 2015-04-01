#include "dcmp2j.h"

vector<datum> le_dados(ifstream& ifile, opcoes opcao) {

    vector<string> target(11);
    target[AZIM]="regional azimuth";
    target[SH]="shear angle";
    target[CH]="channelling angle";
    target[TW]="twist angle";
    target[RHOA]="app rho a";
    target[RHOB]="app rho b";
    target[PHIA]="imped phase a";
    target[PHIB]="imped phase b";
    target[ERR]="error";
    target[SKEW]="skew";
    target[ANIS]="anis";

    string line;
    bool found=false;
    int n;

    // procura o alvo e le o numero de periodos a ser analizado

    ifile.seekg(0); // vai para o inicio do arquivo
    while(getline(ifile,line)) {
        if(line.find(target[opcao])!=string::npos) {
            istringstream istr(line.c_str());
            istr>>n;
            found=true;
            break;
        }
    }

    if(!found) {
        cout<<"nao achou a frase \""<<target[opcao]<<"\" no arquivo de entrada"<<endl;
        exit(1);
    }

    vector<datum> data(n);
    for(int i=0; i<n;i++) {
        getline(ifile, line);
        datum d;
        istringstream istr(line.c_str());
        istr>>d.T>>d.val>>d.max>>d.min;

        data[i]=d;
    }

    return data;
}
