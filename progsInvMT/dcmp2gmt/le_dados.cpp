#include "dcmp2gmt.h"

vector<datum> le_dados(ifstream& ifile, opcoes opcao) {

    vector<string> target(12);
    target[azim]="regional azimuth";
    target[sh]="shear angle";
    target[ch]="channelling angle";
    target[tw]="twist angle";
    target[rhoa]="app rho a";
    target[rhob]="app rho b";
    target[phia]="imped phase a";
    target[phib]="imped phase b";
    target[err]="error";
    target[skew]="skew";
    target[anis]="anis";
    target[phidif]="phadif";

    string line;
    bool found=false;
    int n;

    // procura o alvo e le o numero de periodos a ser analizado
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

        // para channelling e skew, max e min e' sempre zero
        // para ser coerente com as outras saidas, eles recebem o valor central
        if(opcao==ch || opcao==skew) {
            d.max=d.val;
            d.min=d.val;
        }

        data[i]=d;
    }

    return data;
}
