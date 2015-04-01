#ifndef DCMP2J_H
#define DCMP2J_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

typedef struct {
    double T, val, min, max;
} datum;

typedef enum {AZIM, SH, CH, TW, RHOA, RHOB, PHIA, PHIB, ERR, SKEW, ANIS} opcoes;

vector<datum> le_dados (ifstream& ifile, opcoes opcao);
vector<string> le_cab (ifstream& ifile);

#endif // DCMP2J_H
