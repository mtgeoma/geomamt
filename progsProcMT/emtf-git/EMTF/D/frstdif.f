c
c**************************************
c 
        subroutine frstdif(x,nch,npts,lfd,id,areg,npw,r,aa,bb,work)
 
c    computes first difference of nch channels of data, each of length
c       npts

c       routine heavily modified by ping zhang to do adaptive
c      prewhitening using ar approach

c       fixed by egbert to stabilize ar prewhitening coefficients
c
c      still suspect, June 1993 !!!!!
 
        integer npw(*)
        real x(nch,npts),areg(nch,*),r(npts),aa(npts),work(npts)
        real bb(npts)
        logical lfd(nch,*)

        do 10 i=1,nch
        ntem=npw(i)-1
        if( lfd(i,id) ) then
           do 5 j=1,npts
           x(i,j)=x(i,j+1)-x(i,j)
5          continue 
         else if(ntem.ge.1) then
c***           use autoregressive model to find coefficients

c**          find autocrrelation function
         call autocor(x,npts,ntem,nch,i,r,aa)                
c**   solve system equations for coefficients   

         do 30 k=1,ntem
         bb(k)=r(k+1)
 30      continue
         call chdec(aa,work,ntem,ier)
         if(ier.lt.0) then
c            print*,'cholesky decomp failed in prewhitener'
c            print*,'channel = ',i
c            print*,aa(1)
c            print*,aa(2),aa(3)
c            print*,aa(4),aa(5),aa(6)
            do 35 k = 2,npw(i)
35          areg(i,k) = 0.
            areg(i,1) = 1.
            go to 10
         end if

         call ltslv(work,ntem,bb,1)
         ns=ntem/2
         do 100 k1=1,ns
         temp=bb(k1)
         bb(k1)=bb(ntem-k1+1)
         bb(ntem-k1+1)=temp
 100     continue
         ntol=ntem*(ntem+1)/2
         ns=ntol/2
         do 200 k1=1,ns
         temp=work(k1)
         work(k1)=work(ntol-k1+1)
         work(ntol-k1+1)=temp
 200     continue                
         call ltslv(work,ntem,bb,1)
         do 40 j=1,ntem
         areg(i,j+1)=-bb(ntem-j+1)
 40      continue                     
         areg(i,1)=1.   

c**   construct new data series
         do 50 j=1,npts
         xtem=0.
         do 60 k=1,npw(i)
         xtem=xtem+areg(i,k)*x(i,j+npw(i)-k) 
 60      continue
         x(i,j)=xtem
 50      continue    
        end if
 
10      continue
        return
        end
