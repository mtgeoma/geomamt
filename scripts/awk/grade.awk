# le os arquivos x e y (assume que estam nesta ordem!) com o tamanho das grades em metros
# e calcula as coordenadas xy da grade
# ATENÇÃO! Inverte o sentido e origem do eixo Y (Y'=Ymax-Y) mas não muda a numeração das células e a ordem em que são apresentadas!
# usos:
# 1) para obter a coordenada do centro da celula ao longo da coluna 3 do modelo:
# awk -f grade.awk -v x=3 x.dat y.dat
# 2) para obter a coordenada do centro da celula ao longo da linha 5 do modelo:
# awk -f grade.awk -v y=5 x.dat y.dat
# 3) para obter a coordenada do centro da celula da coluna 3 e linha 5 modelo:
# awk -f grade.awk -v x=3 -v y=5 x.dat y.dat
# 4) para obter a coordenada do centro de todas as celulas do modelo:
# awk -f grade.awk x.dat y.dat
# 5) para obter grade do modelo, ja separada com ">" no começo de cada looping interno
# awk -f grade.awk -v saida=grade x.dat y.dat
# para as opções 4 e 5, o looping interno é do y, para ser o contrario, use "-v LI=x"
BEGIN{
  X[1]=0
  Y[1]=0
}
{
  # guarda o distância acumulativa da grade do 1o. arquivo (assume ser o x)
  if(FNR==NR) {
    X[FNR+1]=$1+X[FNR]
    nX=FNR+1
  }
  else { # mudou de arquivo (assume ser o y) se entrar mais que dois arquivos, o script falha!
    Y[FNR+1]=$1+Y[FNR]
    nY=FNR+1
  }
}
END{
  if(x!=0) {
    if(y!=0) { # 3) para obter a coordenada do centro da celula da coluna x e linha y modelo:
      printf "%f %f %d %d\n",(X[x]+X[x+1])/2,(Y[nY]-(Y[y]+Y[y+1])/2),x,y
    }
    else { # 1) para obter a coordenada do centro da celula ao longo da coluna x do modelo:
      for(y=1;y<=nY-1;y++) {
        printf "%f %f %d %d\n",(X[x]+X[x+1])/2,(Y[nY]-(Y[y]+Y[y+1])/2),x,y
      }
    }
  }
  else {
    if(y!=0) { # 2) para obter a coordenada do centro da celula ao longo da linha y do modelo:
      for(x=1;x<=nX-1;x++) {
        printf "%f %f %d %d\n",(X[x]+X[x+1])/2,(Y[nY]-(Y[y]+Y[y+1])/2),x,y
      }
    }
    else {
      if(LI=="x") { # 1o. looping interno é do x
        if(saida=="grade") { # 5) para obter grade do modelo, ja separada com ">" no começo de cada looping interno
          for(y=1;y<=nY;y++) { # grade ao longo dos paralelos
            printf "> y=%d\n",y
            for(x=1;x<=nX;x++) {
              printf "%f %f %d %d\n",X[x],Y[nY]-Y[y],x,y
            }
          }
          for(x=1;x<=nX;x++) { # grade ao longo dos meridianos
            printf "> x=%d\n",x
            for(y=1;y<=nY;y++) {
              printf "%f %f %d %d\n",X[x],Y[nY]-Y[y],x,y
            }
          }
        }
        else { # 4) para obter a coordenada do centro de todas as celulas do modelo:
          for(y=1;y<=nY-1;y++) {
            for(x=1;x<=nX-1;x++) {
	      printf "%f %f %d %d\n",(X[x]+X[x+1])/2,(Y[nY]-(Y[y]+Y[y+1])/2),x,y
            }
          }
        }
      } # fim do looping interno do x
      else { # 1o. looping interno é do y
        if(saida=="grade") { # 5) para obter grade do modelo, ja separada com ">" no começo de cada looping interno
          for(x=1;x<=nX;x++) { # grade ao longo dos meridianos
            printf "> x=%d\n",x
            for(y=1;y<=nY;y++) {
              printf "%f %f %d %d\n",X[x],Y[nY]-Y[y],x,y
            }
          }
          for(y=1;y<=nY;y++) { # grade ao longo dos paralelos
            printf "> y=%d\n",y
            for(x=1;x<=nX;x++) {
              printf "%f %f %d %d\n",X[x],Y[nY]-Y[y],x,y
            }
          }
        }
        else { # 4) para obter a coordenada do centro de todas as celulas do modelo:
          for(x=1;x<=nX-1;x++) {
            for(y=1;y<=nY-1;y++) {
              printf "%f %f %d %d\n",(X[x]+X[x+1])/2,(Y[nY]-(Y[y]+Y[y+1])/2),x,y
            }
          }
        }
      } # fim do looping interno do y
    }
  }
}
