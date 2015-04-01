c
c*****************************************
c
      subroutine edfwts(x,nch,nf,p1,p2,h,w)

      complex h(3),x(nch,nf)
      real edf,p1,p2,w(nf)

         do 10 i = 1,nf
         edf = x(1,i)*conjg(x(1,i))*h(1) + x(2,i)*conjg(x(2,i))*h(3)
     &             +2.*real(conjg(x(2,i))*x(1,i)*h(2))

         if(edf.gt.p2) then
            w(i) = 0.
         else if(edf.gt.p1) then
            w(i) = sqrt(p1/edf)
         else
            w(i) = 1.0
         end if
10       continue

      return
      end
