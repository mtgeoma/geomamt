#ifndef DCMP2GMT_H
#define DCMP2GMT_H

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

typedef enum {azim, sh, ch, tw, rhoa, rhob, phia, phib, err, skew, anis, phidif} opcoes;


vector<datum> le_dados (ifstream& ifile, opcoes opcao);
vector<datum> tofirst (vector<datum> data);

#endif // DCMP2GMT_H
