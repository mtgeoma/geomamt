# uso: ls -l *.ats | awk -f lista_conjuntos.awk
{
    n=split($NF,a,"/")
    if(length(a[n])==12) { # 081E01ZF.ats
	run= substr(a[n],5,2)
	bnd=substr(a[n],8,1)
    }
    else { # 021_V01_C00_R008_TEx_BL_4H.ats
	split(a[n],b,"_")
	run= substr(b[4],2)
	bnd=substr(b[7],1,index(b[7],".ats")-1)
    }
    printf "%s%s %s %s\n",bnd,run,$NF,$5
}
