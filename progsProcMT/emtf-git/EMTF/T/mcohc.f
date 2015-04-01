       subroutine mcohc(yy,xx,xy,b,iq,ip,coh)
c
c       computes multivariate covariance for the multivariate
c         complex linear model
c               Y   =    X     B
c             n x q    n x p  p x q
c              
c            input: xx is x*x, yy is y*y (q times q and p times p
c                hermitian matrices in full storage mode, xy
c                is x*y, B is estimated prediction coefficents,
c            output:  coh is overall coherence of predicted and
c                 observed data (Y)
c            note that this will give an answer for any estimate
c            B, not just the ls estimate; the result returned is the
c             real part of the complex coherence between observed
c             and predicted

      complex xx(ip,ip),yy(iq,iq),xy(ip,iq),b(ip,iq)
      real coh,bxxb,yxb,yt

      yt = 0.
      bxxb = 0.
      do 10 i = 1,iq
         do 5 j = 1,ip
         do 5  k = 1,ip
5        bxxb = bxxb + real(conjg(b(j,i))*xx(j,k)*b(k,i))
      yt = yt + real(yy(i,i))
10    continue

      yxb = 0.
      do 20 i = 1,iq
         do 15 j = 1,ip
         yxb = yxb + conjg(xy(j,i))*b(j,i)
15       continue
20    continue
      coh = yxb*yxb/(yt*bxxb)
 
      if(coh.gt.1) then
         print*,'coh',coh
         print*,'xx',xx
         print*,'yy',yy
         print*,'xy',xy
         print*,'b',b
         print*,'yt,bxxb,yxb',yt,bxxb,yxb
      end if
      return
      end 
