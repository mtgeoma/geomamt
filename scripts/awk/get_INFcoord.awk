{
    if($1=="Input") {
	stn=$NF
    }
    if($1=="Latitude:") {
	n=split($NF,a,",")
	lat=substr(a[1],1,2)+substr(a[1],3)/60
	if(substr(a[2],1,1)=="S") {
	    lat*=-1
	}
    }
    if($1=="Longitude:") {
	n=split($NF,a,",")
	lon=substr(a[1],1,3)+substr(a[1],4)/60
	if(substr(a[2],1,1)=="W") {
	    lon*=-1
	}
    }
}
END{
    printf "%.5f %.5f %s\n",lon,lat,stn
}
