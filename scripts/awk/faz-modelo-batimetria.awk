# uso: awk -f faz-modelo-batimetria.awk -v nx=69 -v ny=77 -v m=9 -v f=8 z.dat topo.grade
# nx=número de colunas
# ny=número de linhas
# m=character que marca o mar
# f=character que marca o fundo oceânico
# assume que z.dat é um arquivo de duas colulas sendo:
# 1a. coluna=espessura da grade
# 2a. coluna=character que marca a qual camada a grade pertence
BEGIN{
    h=0;
    j=1;
}
{
  # assume que primeiro arquivo é z
  # guarda a profundidade acumulativa (1a. coluna)
  # e a camada a qual a grade pertence (2a. coluna)
  if(FNR==NR) {
      Z[FNR]=$1+h;
      h=Z[FNR];
      nz++;
      C[FNR]=$2;
  }
  else { # mudou de arquivo (assume ser o topo.xyz) se entrar mais que dois arquivos, o script falha!
      if(FNR%nx!=0) {
	  i=FNR%nx;
	  B[i,j]=$3
      }
      else {
	  i=nx;
	  B[i,j]=$3;
	  j++;
      }
  }
}
END{
    if (ny!=j-1) {
	printf "ERRO! ny=%d; j=%d\n",ny,j-1
    }
    else {
	for(k=1;k<=nz;k++) {
	    printf "> camada %d (h=%f)\n",k,Z[k];
	    for(j=1;j<=ny;j++) {
		printf " ";
		for(i=1;i<=nx;i++) {
		    if(B[i,j]>0) {
			printf "%s",C[k];
		    }
		    else if(Z[k]<=-B[i,j]) {
			printf "%s",m;
		    }
		    else {
			printf "%s",f;
		    }
		}
		printf "\n";
	    }
	}
    }
}
