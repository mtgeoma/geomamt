{
  if(NF==7) {
    T=$1
    rho=$2
    erho=$3
    phi=$4
    ephi=$5

    if(rho>0&&erho/rho<0.01) {
      erho=rho*0.01
    }
    if(phi>0&&ephi/phi<0.01) {
      ephi=phi*0.01
    }

    if($4>90||$4<0) {
      if(NR==24)
        printf "%11.5E %11.5E %11.5E %11.5E %11.5E 1 0\n",T,rho,erho,phi,ephi
      else
        printf "%11.5E %11.5E %11.5E %11.5E %11.5E 0 0\n",T,rho,erho,phi,ephi
    }
    else if(NR==24)
        printf "%11.5E %11.5E %11.5E %11.5E %11.5E 1 1\n",T,rho,erho,phi,ephi
    else
      printf "%11.5E %11.5E %11.5E %11.5E %11.5E %s %s\n",T,rho,erho,phi,ephi,$6,$7
  }
  else {
    print $0
  }
}
