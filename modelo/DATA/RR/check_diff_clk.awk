# uso: ls ../*B.clk | awk -f check_diff_clk.awk | awk '{system($0)}' | grep identical
{
    clk[NR]=$1
    N=NR
}
END{
    for(i=1;i<N;i++) {
	for(j=i+1;j<=N;j++) {
	    printf "diff -s %s %s\n",clk[i],clk[j]
	}
    }
}
