BEGIN{
  n=1
#  printf "it  rms       model norm\n"
}
{
  if(NR==1) {
      arq=FILENAME
  }
  else if(arq!=FILENAME) {
      arq=FILENAME
      n=1
  }

  if($1=="ITERATION" && $3==n) {
    rms=$6
    getline
    if($1=="SMOOTHING") {
	norm=$7
    }
    else {
	norm=$4
    }
    printf "%4d%10s%10s\n",n,rms,norm
    n++
  }
}
