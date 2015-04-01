c
c*********************************
c
        subroutine ptdist(i1,i2,nmax,nacc,lstart)
 
        logical lstart
 
c       find distance between pointers
 
        if(i2 .eq. i1) then
           nacc = 0
        elseif(i2 .gt. i1) then
           lstart = .true.
           nacc = i2-i1
        else
           nacc = nmax - i1 + i2
        end if
 
        if( .not. lstart ) nacc = 0
        return
        end
