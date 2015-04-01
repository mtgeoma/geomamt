#include "geometria.h"

bool entre_P1P2(ponto P1, ponto P2, ponto P) {
    ponto Min,Max;

    if(P1.x<=P2.x) {
        Min.x=P1.x;
        Max.x=P2.x;
    }
    else {
        Min.x=P2.x;
        Max.x=P1.x;
    }

    if(P1.y<=P2.y) {
        Min.y=P1.y;
        Max.y=P2.y;
    }
    else {
        Min.y=P2.y;
        Max.y=P1.y;
    }

    if(P.x>=Min.x && P.x<=Max.x && P.y>=Min.y && P.y<=Max.y)
        return true;
    else
        return false;
}

vector<float> cruza_P1P2_X(ponto P1, ponto P2, float x) {
    // verifica se segmento P1P2 é cruzado pela grade vertical em X
    // o tamanho do vetor esta entre 0 (não cruza) ate 2 (grade coincidente com segmento)
    float dx=P2.x-P1.x;
    float dy=P2.y-P1.y;
    vector<float> y;

    if(dx==0) { // segmento vertical
        if(P1.x==x) { // segmento coincidente com a grade, "cruza" nos dois extremos
            y.push_back(P1.y);
            y.push_back(P2.y);
        }
        return y;
    }
    else if(dy==0) { // segmento horizontal
        ponto P={x,P1.y};
        if(entre_P1P2(P1,P2,P)) {
            y.push_back(P1.y);
        }
        return y;
    }
    else { // no caso de dx e dy diferentes de zero:
        ponto P;
        P.x=x;
        P.y=(x-P1.x)*dy/dx+P1.y;
        if(entre_P1P2(P1,P2,P)) {
            y.push_back(P.y);
        }
        return y;
    }
}

vector<float> cruza_P1P2_Y(ponto P1, ponto P2, float y) {
    // verifica se segmento P1P2 é cruzado pela grade horizontal em Y
    // o tamanho do vetor esta entre 0 (não cruza) ate 2 (grade coincidente com segmento)
    float dx=P2.x-P1.x;
    float dy=P2.y-P1.y;
    vector<float> x;

    if(dy==0) { // segmento horizontal
        if(P1.y==y) { // segmento coincidente com a grade, "cruza" nos dois extremos
            x.push_back(P1.x);
            x.push_back(P2.x);
        }
        return x;
    }
    else if(dx==0) { // segmento vertical
        ponto P={P1.x,y};
        if(entre_P1P2(P1,P2,P)) {
            x.push_back(P1.x);
        }
        return x;
    }
    else { // no caso de dx e dy diferentes de zero:
        ponto P;
        P.y=y;
        P.x=(y-P1.y)*dx/dy+P1.x;
        if(entre_P1P2(P1,P2,P)) {
            x.push_back(P.x);
        }
        return x;
    }
}

vector<float> sort(vector<float> v) {
    // ordena em ordem crescente o vetor v ...
    for(int i=0;i<v.size()-1;i++) {
        int min=i;
        for(int j=i+1;j<v.size();j++) {
            if(v[j]<v[min]) {
                min=j;
            }
        }
	float temp=v[min];
	v[min]=v[i];
	v[i]=temp;
    }

    // ... e retira repeticoes
    vector<float> u;
    u.push_back(v[0]);
    int iu=0;
    for(int i=1;i<v.size();i++) {
        if(v[i]!=u[iu]) {
            u.push_back(v[i]);
            iu++;
        }
    }

    return u;
}

vector<float> cruza_em_x(poligono pol, float x) {
    // retorna um vetor em ordem crescente com os valores de y em que o poligono é cruzado na ordenada x
    // se o pologono não é cruzado em x, vetor retornado é de tamanho zero
    vector<float> y;
    for(int p=0;p<pol.P.size()-1;p++) {
        ponto P1=pol.P[p];
        ponto P2=pol.P[p+1];
        vector<float> tmp=cruza_P1P2_X(P1, P2, x);
        for(int i=0;i<tmp.size();i++)
            y.push_back(tmp[i]);
    }

    if(y.size()==0) {
        return y;
    }
    else {
        y=sort(y);
        return y;
    }
}

vector<float> cruza_em_y(poligono pol, float y) {
    // retorna um vetor em ordem crescente com os valores de x em que o poligono é cruzado na ordenada y
    // se o pologono não é cruzado em y, vetor retornado é de tamanho zero
    vector<float> x;
    for(int p=0;p<pol.P.size()-1;p++) {
        ponto P1=pol.P[p];
        ponto P2=pol.P[p+1];
        vector<float> tmp=cruza_P1P2_Y(P1, P2, y);
        for(int i=0;i<tmp.size();i++)
            x.push_back(tmp[i]);
    }

    if(x.size()==0) {
        return x;
    }
    else {
        x=sort(x);
        return x;
    }
}
