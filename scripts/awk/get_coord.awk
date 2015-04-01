{
    if(index($0,"LATITUDE")!=0) {
	lat=$NF;
    }
    if(index($0,"LONGITUDE")!=0) {
	lon=$NF;
	print lon,lat,FILENAME;
    }
}
