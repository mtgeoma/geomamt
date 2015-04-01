{
  if($1=="Final" && $2=="Refinement:" && $3=="Chisq") chisq=$5
}
END {
  printf "%s %s\n",substr(FILENAME,4,index(FILENAME,".")-4),chisq
}
