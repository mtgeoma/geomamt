c
c******************************************
c
      subroutine savtf(s,nch,tfunc)

      complex s(*),tfunc(3,*)

      do 5 i = 1,nch-2
      ii = ((i+1)*(i+2))/2+1
      tfunc(1,i) = s(ii)
      tfunc(2,i) = s(ii + 1)
      ii = ii+i+1
      tfunc(3,i) = s(ii)
5     continue
      return
      end
