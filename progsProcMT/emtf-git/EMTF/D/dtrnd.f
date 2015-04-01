c
c******************
c
        subroutine dtrnd(x,nch,npts,lfd,id)
 
        real x(nch,npts)

        logical lfd(nch,*)
 
        do 20 i=1,nch
 
        if(.not. lfd(i,id) ) then
           xx = 0.
           xy = 0.
           yb = (npts+1)/2.

              do 5 j=1,npts,5
              xx = xx + (j - yb) * (j - yb)
              xy = xy + (j - yb) * x(i,j)
5             continue

           b = xy/xx

              do 10 j = 1,npts
              x(i,j) = x(i,j) - b * (j-yb)
10            continue

        end if
 
20      continue
        return
        end
