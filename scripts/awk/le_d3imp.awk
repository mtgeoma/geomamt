# le as impedancias
# com x e y constantes, saida será no formato Jones neste ponto
# awk -f le_d3imp.awk -v x=1 -v y=1
# com x, y e T constantes, saida será o vetor de indução neste ponto
# awk -f le_d3imp.awk -v x=1 -v y=1 -v NT=1
# com x e T constantes, saida será vetores de induçao nesta linha e período
# awk -f le_d3imp.awk -v x=1 -v NT=1
# com y e T constantes, saida será vetores de induçao nesta linha e período
# awk -f le_d3imp.awk -v y=1 -v NT=1
# onde
# x e/ou y marca o ponto/linha que se quer recuperar
# NT numero do periodo que se quer analizar
# use -v e=2 para multiplicar por 2 o modulo do vetor de indução
# ATENÇÃO! Inverte o sentido e origem do eixo Y (Y'=Ymax-Y) mas não muda a numeração das células e a ordem em que são apresentadas!
BEGIN{
  t=0
  nP=0
  P[1,1]=0
  P[1,2]=0
}
{
  # guarda as posicoes de X
  while($1!="X-positions") getline
  getline
  nX=0
  while($1!="Y-positions") {
    for(i=1;i<=NF;i++) {
      nX++
      X[nX]=$i
    }
    getline
  }

  # guarda as posicoes de Y
  if($1=="Y-positions") {
    getline
    nY=0
    while($1!="ZXX") {
      for(i=1;i<=NF;i++) {
        nY++
        Y[nY]=$i
      }
      getline
    }
  }

  do {
    t++
    # guarda ZXX
    getline
    f[t]=$1
    for(cy=1;cy<=nY;cy++) {
      cx=0
      getline
      while(NF!=0) {
        if(NF%2!=0) {
          printf "ERRO! Linha %d com numero de campos impar!\n",NR
          nextfile
        }
        for(c=1;c<=NF;c++) {
          if(c%2==1) {
            cx++
            RZXX[cx,cy,t]=$c
          }
          else {
            IZXX[cx,cy,t]=$c
          }
        }
        getline
      }
    }
    # guarda ZXY
    getline
    if($1!="ZXY") {
      printf "ERRO! Linha %d sem rotulo ZXY esperado!\n",NR
      nextfile
    }
    getline
    if(f[t]!=$1) {
      printf "ERRO! Linha %d com frequencia diferente do bloco anterior!\n",NR
      nextfile
    }
    for(cy=1;cy<=nY;cy++) {
      cx=0
      getline
      while(NF!=0) {
        if(NF%2!=0) {
          printf "ERRO! Linha %d com numero de campos impar!\n",NR
          nextfile
        }
        for(c=1;c<=NF;c++) {
          if(c%2==1) {
            cx++
            RZXY[cx,cy,t]=$c
          }
          else {
            IZXY[cx,cy,t]=$c
          }
        }
        getline
      }
    }
    # guarda ZYX
    getline
    if($1!="ZYX") {
      printf "ERRO! Linha %d sem rotulo ZYX esperado!\n",NR
      nextfile
    }
    getline
    if(f[t]!=$1) {
      printf "ERRO! Linha %d com frequencia diferente do bloco anterior!\n",NR
      nextfile
    }
    for(cy=1;cy<=nY;cy++) {
      cx=0
      getline
      while(NF!=0) {
        if(NF%2!=0) {
          printf "ERRO! Linha %d com numero de campos impar!\n",NR
          nextfile
        }
        for(c=1;c<=NF;c++) {
          if(c%2==1) {
            cx++
            RZYX[cx,cy,t]=$c
          }
          else {
            IZYX[cx,cy,t]=$c
          }
        }
        getline
      }
    }
    # guarda ZYY
    getline
    if($1!="ZYY") {
      printf "ERRO! Linha %d sem rotulo ZYY esperado!\n",NR
      nextfile
    }
    getline
    if(f[t]!=$1) {
      printf "ERRO! Linha %d com frequencia diferente do bloco anterior!\n",NR
      nextfile
    }
    for(cy=1;cy<=nY;cy++) {
      cx=0
      getline
      while(NF!=0) {
        if(NF%2!=0) {
          printf "ERRO! Linha %d com numero de campos impar!\n",NR
          nextfile
        }
        for(c=1;c<=NF;c++) {
          if(c%2==1) {
            cx++
            RZYY[cx,cy,t]=$c
          }
          else {
            IZYY[cx,cy,t]=$c
          }
        }
        getline
      }
    }
    # guarda TZX
    getline
    if($1!="ZZX") {
      printf "ERRO! Linha %d sem rotulo ZZX esperado!\n",NR
      nextfile
    }
    getline
    if(f[t]!=$1) {
      printf "ERRO! Linha %d com frequencia diferente do bloco anterior!\n",NR
      nextfile
    }
    for(cy=1;cy<=nY;cy++) {
      cx=0
      getline
      while(NF!=0) {
        if(NF%2!=0) {
          printf "ERRO! Linha %d com numero de campos impar!\n",NR
          nextfile
        }
        for(c=1;c<=NF;c++) {
          if(c%2==1) {
            cx++
            RTZX[cx,cy,t]=$c
          }
          else {
            ITZX[cx,cy,t]=$c
          }
        }
        getline
      }
    }
    # guarda TZY
    getline
    if($1!="ZZY") {
      printf "ERRO! Linha %d sem rotulo ZZY esperado!\n",NR
      nextfile
    }
    getline
    if(f[t]!=$1) {
      printf "ERRO! Linha %d com frequencia diferente do bloco anterior!\n",NR
      nextfile
    }
    for(cy=1;cy<=nY;cy++) {
      cx=0
      getline
      while(NF!=0) {
        if(NF%2!=0) {
          printf "ERRO! Linha %d com numero de campos impar!\n",NR
          nextfile
        }
        for(c=1;c<=NF;c++) {
          if(c%2==1) {
            cx++
            RTZY[cx,cy,t]=$c
          }
          else {
            ITZY[cx,cy,t]=$c
          }
        }
        getline
      }
    }
    getline
  } while ($1=="ZXX")
  while(FNR==nP+1) {
    nP++
    P[nP,1]=$1
    P[nP,2]=$2
    getline
  }
}
END{
  if(nP==0&&P[1,1]==0&&P[1,2]==0) {
    nP++
    P[1,1]=x
    P[1,2]=y
  }
  # calcula o valor da metade da ultima camada em Y
  y0=Y[1]
  for(j=2;j<=nY;j++) {
    y0=Y[j]-(Y[j-1]+y0)
  }
  # acrescenta metade da ultima camada em Y para obter referencia correta
  nY++
  Y[nY]=Y[nY-1]+y0
  if(NT!=0) { # vetores de inducao
    if(e==0) { # fator de multiplicacao do vetor
      e=1
    }
    for(p=1;p<=nP;p++) {
      if(P[p,1]==0) { # vetores de inducao ao longo da linha y
        j=P[p,2]
        for(i=1;i<=nX;i++) {
          azim=(180.0/3.1415927)*atan2(RTZY[i,j,NT],RTZX[i,j,NT])+180+90
          if(azim>=360){azim-=360}
          if(azim<0){azim+=360}
          mod=sqrt(RTZX[i,j,NT]*RTZX[i,j,NT]+RTZY[i,j,NT]*RTZY[i,j,NT])*e
          printf "%s %s %f %f\n",X[i],Y[nY]-Y[j],azim,mod
        }
      }
      else if(P[p,2]==0) { # vetores de inducao ao longo da linha x
        i=P[p,1]
        for(j=1;j<=nY;j++) {
          azim=(180.0/3.1415927)*atan2(RTZY[i,j,NT],RTZX[i,j,NT])+180+90
          if(azim>=360){azim-=360}
          if(azim<0){azim+=360}
          mod=sqrt(RTZX[i,j,NT]*RTZX[i,j,NT]+RTZY[i,j,NT]*RTZY[i,j,NT])*e
          printf "%s %s %f %f\n",X[i],Y[nY]-Y[j],azim,mod
        }
      }
      else { # vetores de inducao no ponto xy
        i=P[p,1]
        j=P[p,2]
        azim=(180.0/3.1415927)*atan2(RTZY[i,j,NT],RTZX[i,j,NT])+180+90
        if(azim>=360){azim-=360}
        if(azim<0){azim+=360}
        mod=sqrt(RTZX[i,j,NT]*RTZX[i,j,NT]+RTZY[i,j,NT]*RTZY[i,j,NT])*e
        printf "%s %s %f %f\n",X[i],Y[nY]-Y[j],azim,mod
      }
    }
  }
  else if(x!=0&&y!=0) { # formato Jones em x y
    i=x
    j=y
    printf "# dados sinteticos do modelo %s\n",substr(FILENAME,1,index(FILENAME,".")-1)
    printf ">STATION   :stn%02d%02d\n",x,y
    printf ">AZIMUTH   =    %8s\n",90
    printf ">LATITUDE  =    %8s\n",Y[nY]-Y[j]
    printf ">LONGITUDE =    %8s\n",X[i]
    printf ">ELEVATION =    %8s\n\n",0

    # imprime ZXX
    printf "ZXX SI units (ohms)\n"
    printf "%d\n",t
    for(k=1;k<=t;k++) {
      real=RZXX[i,j,k]
      imag=-1*IZXX[i,j,k]
      erro=0.05*sqrt(RZXX[i,j,k]*RZXX[i,j,k]+IZXX[i,j,k]*IZXX[i,j,k])
      printf "%14.4e%14.4e%14.4e%14.4e  1\n",1/f[k],real,imag,erro
    }

    # imprime ZXY
    printf "ZXY SI units (ohms)\n"
    printf "%d\n",t
    for(k=1;k<=t;k++) {
      real=RZXY[i,j,k]
      imag=-1*IZXY[i,j,k]
      erro=0.05*sqrt(RZXY[i,j,k]*RZXY[i,j,k]+IZXY[i,j,k]*IZXY[i,j,k])
      printf "%14.4e%14.4e%14.4e%14.4e  1\n",1/f[k],real,imag,erro
    }

    # imprime ZYX
    printf "ZYX SI units (ohms)\n"
    printf "%d\n",t
    for(k=1;k<=t;k++) {
      real=RZYX[i,j,k]
      imag=-1*IZYX[i,j,k]
      erro=0.05*sqrt(RZYX[i,j,k]*RZYX[i,j,k]+IZYX[i,j,k]*IZYX[i,j,k])
      printf "%14.4e%14.4e%14.4e%14.4e  1\n",1/f[k],real,imag,erro
    }

    # imprime ZYY
    printf "ZYY SI units (ohms)\n"
    printf "%d\n",t
    for(k=1;k<=t;k++) {
      real=RZYY[i,j,k]
      imag=-1*IZYY[i,j,k]
      erro=0.05*sqrt(RZYY[i,j,k]*RZYY[i,j,k]+IZYY[i,j,k]*IZYY[i,j,k])
      printf "%14.4e%14.4e%14.4e%14.4e  1\n",1/f[k],real,imag,erro
    }

    # imprime TZX
    printf "TZX\n"
    printf "%d\n",t
    for(k=1;k<=t;k++) {
      real=RTZX[i,j,k]
      imag=-1*ITZX[i,j,k]
      erro=0.05*sqrt(RTZX[i,j,k]*RTZX[i,j,k]+ITZX[i,j,k]*ITZX[i,j,k])
      printf "%14.4e%14.4e%14.4e%14.4e  1\n",1/f[k],real,imag,erro
    }

    # imprime TZY
    printf "TZY\n"
    printf "%d\n",t
    for(k=1;k<=t;k++) {
      real=RTZY[i,j,k]
      imag=-1*ITZY[i,j,k]
      erro=0.05*sqrt(RTZY[i,j,k]*RTZY[i,j,k]+ITZY[i,j,k]*ITZY[i,j,k])
      printf "%14.4e%14.4e%14.4e%14.4e  1\n",1/f[k],real,imag,erro
    }
  }
}
