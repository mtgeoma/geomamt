#include "dcmp2j.h"

vector<string> le_cab(ifstream& ifile) {

    string line;
    vector<string> cab;

    ifile.seekg(0); // vai para o inicio do arquivo
    getline(ifile,line);
    while(line.find("#")==0 || line.find(">")==0) {
        cab.push_back(line);
        getline(ifile,line);
    }

    return cab;
}
