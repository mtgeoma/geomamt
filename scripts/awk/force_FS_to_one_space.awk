#!/bin/awk -f
# force FS be only one <space>
{
  if(NF>1) {
      printf "%s",$1
    for(c=2;c<=NF;c++) {
      printf " %s",$c
    }
  }
  printf "\n"
}
