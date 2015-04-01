#ifndef MODELOS_H
#define MODELOS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

vector< vector<float> > le_mdl_rebocc(string arqE, float xo);
void escreve_mdl_rebocc(string arqE,string arqS, vector< vector<float> > M);

#endif // MODELOS_H
