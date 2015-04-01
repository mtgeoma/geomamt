#!/bin/awk -f
BEGIN{
  FS=";"
  j=1
  if(taxa==0) {
    taxa=1
  }
}
{
  j*=(1+$NF*taxa/100)
}
END{
  printf "%f\n", (j-1)*100
}
