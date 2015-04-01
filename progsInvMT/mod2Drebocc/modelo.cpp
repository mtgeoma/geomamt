#include "geometria.h"

vector< vector<float> > modelo(vector< vector<float> > M,vector<poligono> pol) {
    // guarda em co a coluna em que se encontra a origem
    int co=0;
    for (int c=1;c<=M[0].size();c++) {
        if(M[0][0]==M[0][c])
            co=c;
    }
    if(co==0) {
        cerr<<"a coordenada do centro da origem xo="<<M[0][0]<<" deve coincidir com algum centro da grade horizontal"<<endl;
        exit(1);
    }

    for(int p=0;p<pol.size();p++) { // para cada poligono...
        vector< vector<float> > cruza_X(M[0].size()-co);
        for(int c=co;c<M[0].size()-co;c++)
            cruza_X[c-co]=cruza_em_x(pol[p],M[0][c]);

        vector< vector<float> > cruza_Y(M.size()-1);
        for(int r=1;r<M.size();r++)
            cruza_Y[r-1]=cruza_em_y(pol[p],M[r][0]);

        for(int r=1;r<M.size();r++) {
            if(cruza_Y[r-1].size()!=0) {
                for(int c=co;c<M[0].size()-co;c++) { // para cada celula
                    if(cruza_X[c-co].size()!=0) {
                        int ng=0; // no. de vezes que a grade cruza o poligono ate chegar ao centro da celula
	                while(M[r][0]>=cruza_X[c-co][ng]) ng++;
                        if(ng%2==1) { // centro dentro
			    M[r][c]=pol[p].rho;
                        }
                    }
                }
            }
        }
    }
    return M;
}
