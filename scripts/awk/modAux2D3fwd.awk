# uso: awk -f modAux2D3fwd.awk arquivo.aux
# transforma o codigo adensado de arquivo.aux para o formato indicial do modelo 3D do d3fwd
BEGIN{
S["1"]="01"
S["2"]="02"
S["3"]="03"
S["4"]="04"
S["5"]="05"
S["6"]="06"
S["7"]="07"
S["8"]="08"
S["9"]="09"
S["a"]="10"
S["b"]="11"
S["c"]="12"
S["d"]="13"
S["e"]="14"
S["f"]="15"
S["g"]="16"
S["h"]="17"
S["i"]="18"
S["j"]="19"
S["k"]="20"
}
{
    if ($1==">") {
	printf " %3d layer #\n",$3;
    }
    else {
	for(i=1;i<=length($1);i++) {
	    printf " %s",S[substr($1,i,1)];
	}
	printf "\n";
    }
}
