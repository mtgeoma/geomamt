# uso: awk -f put_selecao.awk bor004a/selecao.xy bor004a/xy_bor004a_NM > bor004a/xy_bor004a_NM2
BEGIN{
  n=0
  i=0
  getline
  while($1!="data") {
    n++
    r[n]=$1
    p[n]=$2
    getline
  }
  print $0
}
{
  if(NF==7) {
    i++
    printf "%s %s %s %s %s %s %s\n",$1,$2,$3,$4,$5,r[i],p[i]
  }
  else {
    print $0
  }
}
