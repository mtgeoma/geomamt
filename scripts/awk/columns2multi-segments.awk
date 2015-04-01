#!/bin/awk -f
# usage: awk -f columns2multi-segments.awk file.txt
# or
# usage: columns2multi-segments.awk file.txt
#
# where each line in file have the format:
# x y1 ... yn
# assume data columns are ordered.
# If in some line a column data is omitted,
# it will be omitted in output.
# something like:
# x y1 y2 y3
# x y1 y2
# x y1 y2 y3
# x y1    y2 # still y2, besides aligned with y3
# but using: columns2multi-segments.awk -v FS='\t' file.txt
# x\ty1\ty2\ty3
# x\ty1\ty2
# x\ty1\ty2\ty3
# x\ty1\t  \ty3 # now y2 was been omitted (<space> will be removed).
#
# output will be:
# > -Z1
# x y1 1
# : :  :
# > -Zn
# x yn n
# : :  :
#
# value n could be used to colorize data using a cpt file.
# Example:
# makecpt -Cno_green -T1/(n+1)/1 > palette.cpt
# range=$(awk -f columns2multi-segments.awk file.txt | minmax -m -I1)
# columns2multi-segments.awk file.txt | psxy -m -R$range -JX10c -Cpalette.cpt -Wthicker -Ba2f1 -P -K > plot.ps
# columns2multi-segments.awk file.txt | psxy -m -N -R -J -Ss5p -Cpalette.cpt -Wthinner,black -O >> plot.ps
{
  LN=NR; # current line number
  NC[LN]=NF; # number of columns in current line number
  # collect data
  for(c=1;c<=NC[LN];c++) {
    gsub(" ","",$c) # remove any <space>, useful if playing with \t field separator
    d[LN,c]=$c
  }
  # check maximum number of columns
  if(LN==1) {
    MNC=NC[LN]
  }
  else if (NC[LN]>MNC) {
    MNC=NC[LN]
  }
}
END{
  # print output
  for(j=2;j<=MNC;j++) {
    printf "> -Z%d\n",j-1
    for(i=1;i<=LN;i++) {
      if(j<=NC[i]&&length(d[i,j])>0) {
        print d[i,1],d[i,j],j-1
      }
    }
  }
}
