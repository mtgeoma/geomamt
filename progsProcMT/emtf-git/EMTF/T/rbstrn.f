c_________________________________________________________
c 
      subroutine rbstrn(z,nbt,nch,ldf,rdf,lrr,nsmx,ibandlim,
     & leref,ntmx)
c      this routine calls tranls and then corrects error estimates
c       to account for the fraction of downweighted data; also corrects
c       errors to account for correlation between adjacent frequencies
c       caused by data window; currently assumes a pi-prolate window (could
c       be changed  by changing parameter vcorps1 (see egbert and booker,
c       1986 - appendix b for calculation of correction factor)
 
c       z is array of nbt spectral density matrices, nch is number of
c       channels of data at local site
c          , nch2=nch-2, nz = (nch+1)*nch/2 ( = # of  elements
c       of spectral density matrix stored in symmetric stoarage mode (lower
c       triangular half of matrix)), ldf = number of (complex) data points
c       used, rdf = actual degrees of freedom, ldfc = number of adjacent
c       pairs among frequencies stacked - all for each of nbt frequencies
c       (or stations)
 
      complex z(nsmx,*)
      real rdf(ntmx,nbt)
      integer ldf(nbt),ibandlim(2,*)
      logical lrr,leref
      parameter (vcorps1=.17)

      nch2 = nch-2
      nz = ((nch+1)*nch)/2

      write(*,*) 'nsmx,nbt,nch = ',nsmx,nbt,nch
      do k = 1,21
        write(*,*) z(k,1)
      enddo
      if(lrr) then
         do 10 i = 1,nbt
10       call trlsrr(z(1,i),nch,ldf(i)) 
      else if(leref) then
         call rtreref(z,nbt,nz,ldf,nch,nsmx)
      else
         do 20 i = 1,nbt
20       call tranls(z(1,i),nch,ldf(i))
      end if     

c           for each frequency band ...
      do 140 j=1,nbt
      if(ldf(j).gt.2) then

c.........adjust for correlation of adjacent frequencies
         cor =  ibandlim(2,j)-ibandlim(1,j) + 1
         cor = 2.*vcorps1*(cor-1.)/cor
            do 139 k=1,3
            z(k,j)= z(k,j)*(1.+cor)
139         continue
 
c.......adjust errors, coherence to account for fraction of downweighted data
            do 130 i=1,nch2
            ie=(i+2)*(i+3)/2
            if(rdf(i,j).gt.(0.)) then
               cfac=(ldf(j)*2./rdf(i,j))
               err=real(z(ie,j))*cfac*cfac
               coh=aimag(z(ie,j))
               coh = coh/(cfac - (cfac-1.)*coh)
               z(ie,j)=cmplx(err,coh)
            end if
130         continue
         end if
 
140   continue
      return
      end
