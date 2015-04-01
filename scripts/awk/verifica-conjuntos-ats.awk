{
    if(NR==1){
	id=$1
	normal=1
	c=1
	a[c]=$2
	b[c]=$3
    }
    else {
	if(id==$1) {
	    c++
	    a[c]=$2
	    b[c]=$3
	    if(b[1]!=b[c]) {
		normal=0
	    }
	}
	else{
	    if(normal!=1||c!=5) {
		for(i=1;i<=c;i++) {
		    printf "%s %s %s\n",id,a[i],b[i]
		}
	    }
	    id=$1
	    normal=1
	    c=1
	    a[c]=$2
	    b[c]=$3
	}
    }
}
END{
    if(normal!=1||c!=5) {
	for(i=1;i<=c;i++) {
	    printf "%s %s %s\n",id,a[i],b[i]
	}
    }
}
