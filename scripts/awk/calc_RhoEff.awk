# uso: awk -f calc_RhoEff.awk -v n=406 bor.raw
# 'n' e' o bloco a ser calculado (ver relaccao do no. do bloco com o site em bor.txt)
#
BEGIN{
  bloco=sprintf("%04s",n)
}
{
  if($1==bloco){
    getline
    data=$3
    hora=$4
    getline
    getline
    rr=$1
    unit=$2
    printf "> repetition rate: %s %s; data: %s %s\n",rr,unit,data,hora
    getline
    getline
    getline
    while(NF==3){
      if((n=index($1,"M"))!=0){t=substr($1,1,length($1)-1)*1e+6}
      else if((n=index($1,"k"))!=0){t=substr($1,1,length($1)-1)*1e+3}
      else if((n=index($1,"m"))!=0){t=substr($1,1,length($1)-1)*1e-3}
      else if((n=index($1,"u"))!=0){t=substr($1,1,length($1)-1)*1e-6}
      else {t=$1}

      if((n=index($3,"M"))!=0){r=substr($3,1,length($3)-1)*1e+6}
      else if((n=index($3,"k"))!=0){r=substr($3,1,length($3)-1)*1e+3}
      else if((n=index($3,"m"))!=0){r=substr($3,1,length($3)-1)*1e-3}
      else if((n=index($3,"u"))!=0){r=substr($3,1,length($3)-1)*1e-6}
      else {r=$3}
      printf "%e %e\n",(3.9*t),r/(2.3*exp(-0.85))
      getline
    }
    #nextfile
  }
}
