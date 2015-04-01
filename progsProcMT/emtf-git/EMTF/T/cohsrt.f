      subroutine cohsrt(x,ncht,lrr,nfb,ns,ibtf,cohm,num,coht,cohp,nu
     &  ,coh,lpcoh,cfpcoh,xo,per,w,ntfmax)

c    takes wide coherence band array x,
c        computes coherences for each set, outputs only those
c        in the narrow TF band which achieve requisite minimum
c         coherence; nu is number of good sets which are output in xo;

      integer nchsmx
      parameter (nchmx = 20,nchsmx = (nchmx*(nchmx+1))/2) 
      integer ibtf(2),ns,nfb,ncht,i,j,k,ii,num,
     &       nu,nchs,i1,ib,icohhist(0:100,3),nut
      real coh(3,ns),cohm,cohmin,cohp(2),w(nfb)
      complex s(nchsmx),xxinv(2,2),x(ncht,ns,nfb),xo(ncht,*),
     & xx(2,2),xy(2,2),b(2,2),det,yy(2,2)
      logical lrr,lpcoh,ltemp
      character*20 cfpcoh

c      num is minimum number of points desired for tf estiamtion
c      nut is target number of points for tf est.
c            this is determined from total available
      nbtf = ibtf(2) - ibtf(1) + 1
      nut = nint(cohp(1)*float(ns*nbtf)**cohp(2))
      nut = min(ns*nbtf,nut)
      nut = max(nut,num)

c............compute average power vs freq. to allow rough prewhitening
      call pwrvf(x,ncht,ns,nfb,w,xo)

      if(lpcoh) then
         open(unit=97,file=cfpcoh)
      end if
      nchs = (ncht*(ncht+1))/2
      do 100 i = 1,ns
c..........compute broad band sdm for each set
         do 5 j = 1,nchs
5        s(j) = (0.,0.)
         do 20 j = 1,nfb
20       call wstack(x(1,i,j),w(j),ncht,s)

c********************BLOCK FOR SINGLE STATION DATA ******************
       ltemp = .true.
       if(ltemp) then
c      if(.not.lrr) then
c................straight coherence between E and H
c       here we are finding (1,2) coherence of each channel and 
c         (3) multivariate coherence of both together (i.e.,
c         fraction of total signal which is predicted
         
c........form full complex cross product matrices in a form
c        suitable for simple computation
         xx(1,1) = s(1)
         xx(2,1) = conjg(s(2))
         xx(1,2) = s(2)
         xx(2,2) = s(3)
         xy(1,1) = s(7)
         xy(1,2) = s(11)
         xy(2,1) = s(8)
         xy(2,2) = s(12)
         yy(1,1) = s(10)
         yy(1,2) = s(14)
         yy(2,1) = conjg(s(14))
         yy(2,2) = s(15)
      
c.......solve ls problem using Kramers rule to invert
c          x*x matrix
         det=xx(1,1)*xx(2,2)-xx(1,2)*xx(2,1)
         xxinv(1,1)=xx(2,2)/det
         xxinv(2,2)=xx(1,1)/det
         xxinv(1,2)=-xx(1,2)/det
         xxinv(2,1)=-xx(2,1)/det
         call matmultc(xxinv,xy,b,2,2,2)

c........find coherences
         call mcohc(yy(1,1),xx,xy(1,1),b(1,1),1,2,coh(1,i))
         call mcohc(yy(2,2),xx,xy(1,2),b(1,2),1,2,coh(2,i))
         call mcohc(yy,xx,xy,b,2,2,coh(3,i))
      else
c       for remote ref :
c.......coherence between R and H
c       here we are finding multivariate coherence of both
c       refeerence and local mag channels together (i.e.,
c         fraction of total signal which is predicted)
         
c..........form full complex cross product matrices in a form
c          suitable for simple computation
         xx(1,1) = s(1)
         xx(2,1) = conjg(s(2))
         xx(1,2) = s(2)
         xx(2,2) = s(3)
         xy(1,1) = s(16)
         xy(1,2) = s(22)
         xy(2,1) = s(17)
         xy(2,2) = s(23)
         yy(1,1) = s(21)
         yy(1,2) = s(27)
         yy(2,1) = conjg(s(27))
         yy(2,2) = s(28)
         
c..........solve ls problem using Kramers rule to invert
c          x*x matrix
         det=xx(1,1)*xx(2,2)-xx(1,2)*xx(2,1)
         xxinv(1,1)=xx(2,2)/det
         xxinv(2,2)=xx(1,1)/det
         xxinv(1,2)=-xx(1,2)/det
         xxinv(2,1)=-xx(2,1)/det
         call matmultc(xxinv,xy,b,2,2,2)
 
c............find multivariate coherence
         call mcohc(yy,xx,xy,b,2,2,coh(3,i))

      end if
100   continue

c*************** OUTPUT COHERENCE ETC TO FILE ************************
      if(lpcoh) then
c.........print out coherence vs. set #
         write(97,*) 'period =   ',per,ibtf(1),ibtf(2)
         write(97,'(i5,3f8.3)') (i,(coh(k,i),k=1,3),i=1,ns)
      end if

c************* FIND QUANTILES OF COHERENCE ***************************

         do 105 i = 0,100
         do 105 j = 1,3
105      icohhist(i,j) = 0

         do 110 i = 1,ns
            do 108 k = 1,3
            ii = int(coh(k,i)*100)
            if(( ii.lt.0).or.(ii.gt.100)) then
               print*,'!!!!!!!!! ii out of bounds !!!!!!!! coh= '
     &               ,coh(k,i),'  i= ',i,' k= ',k
               ii = 0
               stop
            end if
            icohhist(ii,k) = icohhist(ii,k) + 1
108         continue
110      continue

         do 115 j = 1,3
         do 115 i = 100,1,-1
115      icohhist(i-1,j) = icohhist(i-1,j) + icohhist(i,j)

      if(lpcoh) then
         write(97,'(i4,3i5)') (i,(icohhist(i,j),j=1,3),i=1,100)
      end if

c*********DETERMINE MINIMUM COHERENCE TO USE****************************
      i0 = 100*coht
         print*,'nbtf',nbtf
         print*,'i0,coht,icohhist(i0,3)',i0,coht,icohhist(i0,3)
      if(icohhist(i0,3)*nbtf.ge.nut) then
c........minimum coherence is target coherence
         cohmin = coht
      else 
         i1 = 100*cohm
          print*,'i1,cohm,icohhist(i1,3)',i1,cohm,icohhist(i1,3)
         if(icohhist(i1,3)*nbtf.ge.nut) then
c...........determine minimum coherence from target number of points
               do 118 i = i0,i1,-1
               if(nbtf*(icohhist(i,3)).ge.nut) then 
                  cohmin = float(i)/100.
                  go to 119
               end if
118            continue
119         continue
         else
            if(icohhist(i1,3)*nbtf.ge.num) then
c..............minimum coherence is input parameter cohm
               cohmin = cohm
            else
c..............determine minimum coherence from minimum number of data num
                  do 120 i = i1-1,1,-1
                  if(nbtf*icohhist(i,3).ge.num) then
                     cohmin = float(i)/100.
                     go to 121
                  end if
                  cohmin = 0.
120               continue
121            continue
            end if
         end if
      end if

      print*,'cohmin,ns,ibt ',cohmin,ns,ibtf

c************* PREPARE FC ARRAY XO FOR OUTPUT **********************
      nu = 0
         do 200 i = 1,ns
         if(coh(3,i).ge.cohmin) then
            do 150 ib =  ibtf(1),ibtf(2)
            nu = nu+1
               do 140 j = 1,ncht
               xo(j,nu) = x(j,i,ib)
140            continue 
150         continue
         end if
         if(nu.eq.ntfmax) return
200      continue
         print*,'nu   ',nu
      return
      end
