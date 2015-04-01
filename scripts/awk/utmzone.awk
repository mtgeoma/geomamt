BEGIN{
 z=1
 l0=-177
}
{
  while($1<l0-3||$1>l0+3) {
    z++
    l0+=6
  }
}
END{
  print z, l0
}
