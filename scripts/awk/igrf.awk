# uso:
# awk -f igrf.awk -v k=0 a=2008.0 arquivo.txt
# onde arquivo.txt tem os dados de longitude, latitude e um rótulo
{
    printf "geomag60 /usr/local/share/IGRF10.cof %s D K%s %s %s | awk \47{if(NR==17)print \"%s\",$0}\047\n",a,k,$2,$1,$3;
}
