# uso: awk -f D3fwd2XY.awk -v layer=1 mod3dSE3 -v saida=texto
# saida=texto: x y "texto da celula"
BEGIN{
  # no. de camadas
  getline
  nx=$1
  ny=$2
  nz=$3
  # le o comprimento e guarda a posição central da ordenada x da celula
  Xi_1=0
  ci_1=0
  i=0
  while(i!=nx) {
  getline
    for(c=1;c<=NF;c++) {
      i++
      X[i]=Xi_1+($c+ci_1)/2
      Xi_1=X[i]
      ci_1=$c
    }
  }
  # le o comprimento e guarda a distância acumulativa da ordenada y da celula
  Yj_1=0
  j=0
  while(j!=ny) {
  getline
    for(c=1;c<=NF;c++) {
      j++
      Y[j]=$c+Yj_1
      Yj_1=Y[j]
    }
  }
  # calcula a posição central da ordenada y da celula,
  # Inverte o sentido e origem do eixo Y (Y'=Ymax-Y)
  # mas não muda a numeração das células
  Yj_1=0
  Ymax=Y[ny]
  for(j=1;j<=ny;j++) {
    Yj=Y[j]
    Y[j]=(Ymax-(Y[j]+Yj_1)/2)
    Yj_1=Yj
  }

  # guarda as profundidades das camadas em Z
  # atualmente não há uso
  k=0
  Zk_1=0
  while(k!=nz) {
  getline
    for(c=1;c<=NF;c++) {
      k++
      Z[k]=$c+Zk_1
      Zk_1=$c
    }
  }
  # procura a camada que se quer visualizar
  while($1!=layer||$2!="layer") {
    getline
  }
  # guarda os valores da camada
  for(j=1;j<=ny;j++) {
    i=0
    while(i!=nx) {
      getline
      for(c=1;c<=NF;c++) {
        i++
        L[j,i]=$c
      }
    }
  }
}
END{
  for(j=1;j<=ny;j++) {
    for(i=1;i<=nx;i++) {
      printf "%f %f %s\n",X[i],Y[j],L[j,i]
    }
  }
}
