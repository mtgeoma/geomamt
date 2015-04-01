{
  for (c=1;c<=NF;c+=2) { # pula colunas pares de seleção (não afeta o processo)
    split($c,B,"/") # algo como MT00256/apa046a081_02C.zss
    win=sprintf("%d",substr(B[1],3)) # janela de processamento
    ext=sprintf("%s",substr(B[2],length(B[2])-3)) # guarda a extensão
    B[2]=substr(B[2],1,length(B[2])-4) # remove a extensão
    if(ext==".zrr") { # assume que não há nada além de RR
      split(B[2],C,"RR")
      printf "%s.asc %5s rr;%s.asc\n",C[1],win,C[2]
    }
    else {
      n=split(B[2],C,"_")
      if(n==2) { # ADU06
	op=substr(C[2],4)
	C[2]=substr(C[2],1,3)
      }
      else { # ADU07
	op=substr(C[3],index(C[3],"H")+1)
	C[3]=substr(C[3],1,index(C[3],"H"))
      }
      bad=""
      if(index(op,"bad")!=0) {
	gsub("bad","",op)
	bad=";bad"
      }
      if(n==2) {
	asc=sprintf("%s_%s.asc%s",C[1],C[2],bad)
      }
      else {
	asc=sprintf("%s_%s_%s.asc%s",C[1],C[2],C[3],bad)
	asc=sprintf("%-25s",asc)
      }
      if(length(op)==0) {
	printf "%s %5s ss\n",asc,win
      }
      else {
	printf "%s %5s ss;%s\n",asc,win,op
      }
    }
  }
}
