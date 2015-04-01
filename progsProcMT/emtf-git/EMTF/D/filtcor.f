c
c***********************************************
c
        subroutine filtcor(x,nch,nfrq,idl,rnrmt,nd,areg,npw,nn,lfd)
c    divides fourier coefficients by transfer function of filters
c    uses table calculated by routine fcorsu

c       replaces routine rnormc
        integer npw(*)
        real x(nch,2,nfrq),areg(nch,*)
        logical lfd(nch,*)
        complex t,rnrmt(nch,nd,*),xcor,wr   !w,
c        print*,'a1= ',areg(1,1),'  a2= ',areg(1,2)
c        print*,'a3= ',areg(1,3),'  a4= ',areg(1,4)
        pian=2.*3.141592654/nn
        do 10 i=1,nfrq
           do 5 j=1,nch   
             if(lfd(j,idl)) then
                xcor=(1.,0.)
             else
                xcor=(0.,0.)
                   do 3 k=1,npw(j)                 
                   wr=cmplx(0.,(npw(j)-k)*pian*i)   
                   xcor=xcor+areg(j,k)*cexp(wr)
 3                 continue
             end if
           t=rnrmt(j,idl,i)*cmplx(x(j,1,i),x(j,2,i))/xcor
           x(j,1,i)=real(t)
           x(j,2,i)=aimag(t)
5          continue

10      continue
        return
        end
