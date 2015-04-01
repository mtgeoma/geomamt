# uso: awk -f torhoplus+E.awk -v form=edi -v c=3118 -v s=conducting -v d=xy -v begin=10/3 -v end=31/20 -v err=3 -v file=sjr01a_a+10Tg.dat/sjr001_a+10Tg.dat
# form: formato (edi, jones, egbert)
# c: criterion
# s: surface
# d: data (xy | yx)
# b: begin (no. do 1o. periodo)
# e: end   (no. do ultimo periodo)
# f: file
BEGIN{
  n=split(file,f,"/")
  split(begin,b,"/")
  split(end,e,"/")
  printf "data *\n"
  if(c!="")
    printf "task *\n"
#  else if(n==2)
#    printf "bands %d %d %d %d -1\n",1,e[1]-b[1]+1,e[1]-b[1]+2,e[1]-b[1]+2+e[2]-b[2]
  printf "surface %s\n",s
  printf "root %s_%s\n",d,substr(f[1],1,index(f[1],".")-1)
  printf "model\n"
  printf "matrix\n"
  if(c!="")
    printf "criterion %s\n",c
  printf "execute\n"
  if(d=="xy") shift=0.
  else shift=180.
  for(i=1;i<=n;i++) {
    cmd1=sprintf("transferfunction %s -F%s -Dr%s -Egamble | ",f[i],form,d)
    cmd2=sprintf("awk \047{if(NR>=%s&&NR<=%s)printf \042%s\042,1./$1,$2,($2-$3)*%s,$5+%s,($5-$6)*%s}\047",b[i],e[i],"%11.5E %11.5E %11.5E %11.5E %11.5E 0 1\\n",err,shift,err)
    cmd=sprintf("%s%s\n",cmd1,cmd2)
    system(cmd)
  }
  printf "0\n"
  if(c!="") {
    printf "\n"
    for(i=1;i<=n;i++) {
      cmd1=sprintf("transferfunction %s -F%s -Dr%s -Egamble | ",f[i],form,d)
      cmd2=sprintf("awk \047{if(NR>=%s&&NR<=%s)printf \042%s\042,1./$1}\047",b[i],e[i],"%11.5E rho phi\\n")
      cmd=sprintf("%s%s\n",cmd1,cmd2)
      system(cmd)
    }
    printf "0\n"
  }
}
