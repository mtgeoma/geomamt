
        subroutine cleanrec(ix,n,nchp1,ich,scale,ixout,ixplt,pscale
     1  ,nmad,iuscale)
 
c       point cleaning routine; uses running meadians to look for
c       large isolated glitches (e.g. parity errors); input is nchp1
c       channels of data (first channel is assumed to be record number)
c       stored in array ix; the number of points to be cleaned is n
c       (this is also the length of the running median window) the
c       points to be cleaned should have (n-1)/2 adjacent points on either
c       side: the first record to be cleaned should be ix(.,1) - i.e. in
c       the call to this routine if records k to k+n-1 in array ixx are
c       to be cleaned (here ix is array in calling program) then the
c       first argument in the call statement should be ixx(1,k));  scale
c       is input and is used to find the minimum deviation from the median
c       required for a point to be considered bad (algorithm will adjust
c       scale if the median absolute deviation (mad) is larger); the output
c       is (1) ixout - array of n cleaned records; (2) xplt  -  2 x nchp1 array
c       for plotting - for each channel the two points are the meadian
c       minus and plus one mad (first channel is median record number)
 
c       if iuscale = 0 works as described above; if not the mad is replaced
c       by the ith order statistic of the absolute deviations; for use with
c       cleaning of first differences use iuscale = n-1
 
        parameter (nmax =20)
 
c       nmax is maximum n which can be used with routine; nmad (multiplied
c       times maximum of scale or mad) is maximum deviation allowed before
c       point is fixed
 
        integer ix(nchp1,*),ixout(nchp1,n),iwrk1(nmax),iwrk2(nmax)
     1       ,ixplt(2)
        real scale(nchp1),pscale(nchp1)
 
        nmed = (n+1)/2
        nmed1 = nmed - 1
        if(iuscale.eq.0) iuscale = nmed
 
c       copy data into array iwrk1, sort and check
 
           do 10 j = 1,n
10         iwrk1(j) =  ix(ich,j)
 
        call isort(iwrk1,n)
 
        med = iwrk1(nmed)
 
           do 15 j = 1,n
15         iwrk2(j) = iabs(iwrk1(j) - med)
 
        call isort(iwrk2,n)
 
        mad = iwrk2(iuscale)
        maxd = iwrk2(n)
        madmax = max(mad,nmad)
 
c       compute plotting coordinates
 
        pmad = pscale(ich)*mad*2.
        ixplt(1) = med
        ixplt(2) = pmad
 
        if(maxd.le.madmax) then
 
c       maximum deviation is small enough to consisder all points ok;
c       compute outputs for this channel
 
           do 20 j = 1,n
20         ixout(ich,j) = ix(ich,j)
 
           go to 100
           end if
 
c       else need to check each point; first center point
 
        rscale = max(scale(ich),float(mad))*nmad
        itemp = ix(ich,nmed)-med
        if(iabs(itemp).gt.rscale) then
c     fix point
 
        iii = isign(1,itemp)
        ixout(ich,nmed) = (iwrk1(nmed)+iwrk1(nmed-iii))/2
cpin    print*, 'point fixed; ch. = ',ich,'rec. =',ix(1,nmed)
cpin    print*,'rscale',rscale,'deviation = ',itemp
cpin    print*,'iwrk2',iwrk2
 
       else
           ixout(ich,nmed) = ix(ich,nmed)
       end if
 
c       now do other points
 
        do 50 l =1,nmed1
           l1 = 1-l
           l2 = n-l
              do 30 j = l1,l2
30            iwrk1(j+l) = ix(ich,j)
 
           itemp = iwrk1(nmed)
           call isort(iwrk1,n)
 
           med = iwrk1(nmed)
 
              do 35 j = 1,n
35            iwrk2(j) = iabs(iwrk1(j) - med)
 
           call isort(iwrk2,n)
 
           mad = iwrk2(iuscale)
 
           rscale = max(scale(ich),float(mad))*nmad
           itemp = itemp-med
           if(iabs(itemp).gt.rscale) then
 
c     fix point
 
              iii = isign(1,itemp)
cpin          print*, 'point fixed; ch. = ',ich,'rec.=',ix(1,nmed-l)
cpin          print*,'rscale',rscale,'deviation = ',itemp
cpin          print*,'iwrk2',iwrk2
 
              ixout(ich,nmed-l) = (iwrk1(nmed)+iwrk1(nmed-iii))/2
          else
              ixout(ich,nmed-l) = ix(ich,nmed-l)
          end if
 
           l1 = 1+l
           l2 = n+l
              do 40 j = l1,l2
40            iwrk1(j-l) = ix(ich,j)
 
           itemp = iwrk1(nmed)
           call isort(iwrk1,n)
 
           med = iwrk1(nmed)
 
              do 45 j = 1,n
45            iwrk2(j) = iabs(iwrk1(j) - med)
 
           call isort(iwrk2,n)
 
           mad = iwrk2(iuscale)
 
           rscale = max(scale(ich),float(mad))*nmad
           itemp = itemp-med
           if(iabs(itemp).gt.rscale) then
 
c     fix point
 
cpin          print*, 'point fixed; ch. = ',ich,'rec.=',ix(1,nmed+l)
cpin          print*,'rscale',rscale,'deviation = ',itemp
cpin          print*,'iwrk2',iwrk2
    
              iii = isign(1,itemp)
              ixout(ich,nmed+l) = (iwrk1(nmed)+iwrk1(nmed-iii))/2
          else
              ixout(ich,nmed+l) = ix(ich,nmed+l)
          end if
50      continue
 
100     continue
 
        return
        end
