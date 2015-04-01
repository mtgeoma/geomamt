BEGIN{
  type=""
  charSize=""
  memA=""
}
{
  if(NF>0) {
    if ($3=="//") {
      type=$1
      label=substr($2,1,length($2)-1)
      charSize=""
      mem="0x"substr($4,1,3)
    }
    else if ($2=="//") {
      label=substr($1,1,length($1)-1)
      charSize=""
      mem="0x"substr($3,1,3)
    }
    else if ($4=="//") {
      type=$1
      label=$2
      charSize=substr($3,1,length($3)-1)
      mem="0x"substr($5,1,3)
    }

    if($1=="//") {
      print substr($0,2)
    }
    else {
      if(length(memA)>0) {
	printf "%-7s %-17s %s %3d",typeA,labelA,memA,strtonum(mem)-strtonum(memA)
	if(length(charSizeA)>0) {
	  printf " %5s",charSizeA
	}
	printf "\n"
      }
      typeA=type
      i=match(label,/[[:upper:]]/) # first upper in label
      labelA=substr(label,i)
      memA=mem
      charSizeA=charSize
    }
  }
}
END{
  printf "%-7s %-17s %s %3d",typeA,labelA,memA,strtonum("0x400")-strtonum(memA)
  if(length(charSizeA)>0) {
    printf " %5s",charSizeA
  }
  printf "\n"
}
