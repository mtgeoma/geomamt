# uso: awk -v dp=0.10 -v ds=0.05 -v tp=0.04 -f ~/awk/modEM-ErrorFloor.awk data.dat > dataEF.dat
{
  if(($8=="ZXX"||$8=="ZYY")&&(dp>0)) {
    if($11/sqrt($10**2+$9**2)<dp) {
      $NF=sprintf("%.6E",sqrt($10**2+$9**2)*dp)
    }
  }
  else if(($8=="ZXY"||$8=="ZYX")&&(ds>0)) {
    if($11/sqrt($10**2+$9**2)<ds) {
      $NF=sprintf("%.6E",sqrt($10**2+$9**2)*ds)
    }
  }
  else if(($8=="TX"||$8=="TY")&&(tp>0)) {
    if($11<tp) {
      $NF=sprintf("%.6E",tp)
    }
  }
  print $0
}
