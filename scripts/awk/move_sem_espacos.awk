# uso: ls | awk -f ~/awk/move_sem_espacos.awk | awk '{system($0)}'
{
    if(NF>1) {
	printf "mv \"%s\" ",$0;
	for(i=1;i<=NF;i++)printf "%s",$i;
	printf "\n";
    }
}
