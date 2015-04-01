#ifndef GEOMETRIA_H
#define GEOMETRIA_H
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

//const float meu_NAN=-999999999;

typedef struct {
  float x,y;
} ponto;

typedef struct {
  float rho;
  vector<ponto> P;
} poligono;

typedef struct {
  float a,b,c;
} reta;

vector<poligono> le_poligono(string arqP, float xo);
bool entre_P1P2(ponto P1, ponto P2, ponto P);
vector<float> cruza_P1P2_X(ponto P1, ponto P2, float x);
vector<float> cruza_P1P2_Y(ponto P1, ponto P2, float y);
vector<float> sort(vector<float> v);
vector<float> cruza_em_x(poligono pol, float x);
vector<float> cruza_em_y(poligono pol, float y);
vector< vector<float> > modelo(vector< vector<float> > M,vector<poligono> pol);

#endif // GEOMETRIA_H
