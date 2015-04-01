# duas primeiras colunas com coordenadas geográficas
# se for no formato dd:mm:ss, converte para dd.ddddd
# awk -f gmt_dms2deg.awk arquivo_gmt_dms.dat
{
    for(c=1;c<=2;c++) {
	n=split($c,d,":")
	if(n==3){ # formato dd:mm:ss
	    deg=int(d[1])
	    min=int(d[2])
	    sec=d[3]
	    if(deg<0) {
		printf " %10.5f",-1*(-1*deg+min/60.0+sec/3600.0)
	    }
	    else {
		printf " %10.5f",(deg+min/60.0+sec/3600.0)
	    }
	}
	else { # não faz nada
	    printf " %s",$c
	}
    }
    for (c=3;c<=NF;c++) {
	printf " %s",$c
    }
    printf "\n"
}
