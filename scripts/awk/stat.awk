BEGIN {
 sx=0; n=0;
}
{
  sx+=$1
  n++
  x[n]=$1
}
END {
  med=sx/n
  # Calculo da variancia em dois passos
  # CHAN TF, GOLUB GH, LEVEQUE RJ. AMERICAN STATISTICIAN 37 (3): 242-247 1983
  var=0.
  Sdelta=0.
  for(i=1;i<=n;i++) {
    delta=x[i]-med
    Sdelta+=delta
    var+=delta*delta
  }
  var=(var-Sdelta*Sdelta/n)/(n-1)
  sd=sqrt(var)

    printf "%.11f %.11f %.11f\n", med, sd, sd/med
}
