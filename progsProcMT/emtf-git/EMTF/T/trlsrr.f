c________________________________________________________________
c
      subroutine trlsrr(s,nch,nu)
 
c      new version; G. Egbert March 1991
c       computes transfer functions with remote reference 
c       technique from spectral density matrix
c
c       input: s = spectral density matrix in symmetric storage mode
c         (complex) for one frequency
c
c               nu = number of frequencies
c
c               nch= # of channels at local station; total channels
c                      used by routine is nch+2
c
c       output: s contains
c            s(1)-s(3) = [RH]-1 [RR] [HR]-1   (i.e., SIG_S for RR case)
c                     ( symmetric storage mode)
c            s(4),s(5)
c            s(7), s(8) = ransfer functions (h first d second)
c            s(11),s(12)         (nch-2 total TFs)

c            s(6), s(10), s(15) = crude error variances (real part)
c                and coherence (imag part)
c         (and the rest of the diagonal block... allow for
c           computation of full error covariance)
c      (note: this is for complex tranfer function; error variance
c      for real and imaginary parts are each one half of this; if
c      iterative robust scheme is used some corrections are needed
c      to account for downweighted data; also some corrections are
c      needed to account for the correlation between adjacent
c      frequencies caused by windowing)
c       s(6)*(s(1-3)) is covariance matrix of transfer function 
c          estimates
c
ccc    NOTE: works for any number nch of local channels (> 2)
c          on input the  sdm looks like the conjugate of:

c    # of columns= 2   nch-2  2 

c         # rows            
c           2     hh
c        nch-2    eh   ee
c           2     rh   re    rr

c         the upper diagonal parts of hh, ee rr are not included
c         hh means avg ( h x conjg(h) transpose )
c          for complex ls we want avg( conjg(h) x h transpose )
c          NEED TO BE CAREFUL ABOUT THIS!!!!!!
  
      parameter (nchmx=20)
      complex s(*),rh(2,2),rr(2,2),re(2,nchmx),det,temp,rt(2,nchmx),
     &    hh(2,2),he(2,nchmx),err(100),wrk(100),cwrk(20)
      integer nu,sindex
 
      if(nu.lt.2) then
c       can't do anyting
          go to 60
      end if

      ns = (nch*(nch+1))/2
      nch2 = nch-2
      n = max(nu-2,1)

cnew      write(*,*) 'n,ns,nch2 = ',n,ns,nch2

c>>>>>>>>>> move cross products into matrices hh,he,rh,rr,re
c.........hh
      hh(1,1) = s(1)
      hh(1,2) = s(2)
      hh(2,1) = conjg(s(2))
      hh(2,2) = s(3)
c.........rh
      rh(1,1) = conjg(s(ns+1))
      rh(1,2) = conjg(s(ns+2))
      rh(2,1) = conjg(s(ns+nch+2))
      rh(2,2) = conjg(s(ns+nch+3))
c.........rr
      rr(1,1) = s(sindex(nch+1,nch+1))
      rr(1,2) = s(sindex(nch+2,nch+1))
      rr(2,1) = conjg(rr(1,2))
      rr(2,2) = s(sindex(nch+2,nch+2))
c........re, he
      do i = 3,nch
         he(1,i-2) = s(sindex(i,1))
         he(2,i-2) = s(sindex(i,2))
         re(1,i-2) = conjg(s(sindex(nch+1,i)))
         re(2,i-2) = conjg(s(sindex(nch+2,i)))
      enddo

c>>>>>>>>>>>> invert RH using Kramers rule:
      det = rh(1,1)*rh(2,2)-rh(2,1)*rh(1,2)
      if(cabs(det).eq.0) then
        write(*,*) 'det = 0'
        write(*,*) 'rh',rh
        go to 60
      endif
      temp = rh(1,1)
      rh(1,1) = rh(2,2)/det
      rh(2,2) = temp/det
      rh(1,2) = -rh(1,2)/det
      rh(2,1) = -rh(2,1)/det

c>>>>>>>>>> compute TFs and put back into eh block
      call matmultc(rh,re,rt,2,2,nch2)
      do i =3,nch
         s(sindex(i,1)) = rt(1,i-2)
         s(sindex(i,2)) = rt(2,i-2)
      enddo

c>>>>>>>>> compute error matrix
c...............error scales
      call mmacb(he,rt,err,nch2,2,nch2)
      do i = 3,nch
         cwrk(i) = s(sindex(i,i))
         do j = 3,i
            i1 = i-2
            j1 = j-2
            s(sindex(i,j)) = s(sindex(i,j))-conjg(err((j1-1)*nch2+i1))
     &                      - err((i1-1)*nch2+j1)
         enddo
      enddo
      call mmacb(rt,hh,wrk,nch2,2,2)
      call matmultc(wrk,rt,err,nch2,2,nch2)

      do i = 3,nch
         do j = 3,i
            i1 = i-2
            j1 = j-2
            s(sindex(i,j)) = s(sindex(i,j))+conjg(err((j1-1)*nch2+i1))
            s(sindex(i,j)) = s(sindex(i,j))/n
         enddo
      enddo
 
      do i = 3,nch
         call mcohc(cwrk(i),hh,he(1,i-2),rt(1,i-2),1,2,coh)
         s(sindex(i,i)) = cmplx(real(s(sindex(i,i))),coh)
      enddo

c..........signal power inverse
      call matmultc(rh,rr,wrk,2,2,2)
      call mmabc(wrk,rh,err,2,2,2)

      do i = 1,2
         do j = 1,i
            s(sindex(i,j)) = err((j-1)*2+i)
         enddo
      enddo
      return

60    do k=1,ns
         s(k)=(-999.,-999.)
      enddo
      return
      end
c___________________________________________________________________
c
      subroutine tranls(s,nch,nu)
 
      parameter (nchmx = 20)
c       computes LS transfer functions from spectral density matrix
c
c       input: s(*) = spectral density matrix in symmetric storage mode
c         (complex)
c
c               nu = number of frequencies in each band
c
c               nch= # of channels
c
c       output: s(*) contains:
c         TFs and error parameters as described in trlsrr
c        (except, s(1)-s(3) contains HH* inv (still SIG_S ... but
cc               now for single station case)

      complex s(*),det,rt(2,nchmx),hh(2,2),he(2,nchmx),err(100)
     &       ,wrk(100),hhinv(2,2)
      integer nu,sindex
 
      if(nu.lt.2) then
c       can't do anyting
          go to 60
      end if

      ns = (nch*(nch+1))/2
      nch2 = nch-2
      n = max(nu-2,1)

cnew      write(*,*) '%*^&**$%$^   nch2 = ',nch2

c>>>>>>>>>> move cross products into matrices hh,he,rh,rr,re
c.........hh
      hh(1,1) = s(1)
      hh(1,2) = s(2)
      hh(2,1) = conjg(s(2))
      hh(2,2) = s(3)

c........he
      do i = 3,nch
         he(1,i-2) = s(sindex(i,1))
         he(2,i-2) = s(sindex(i,2))
      enddo

c>>>>>>>>>>>> invert HH using Kramers rule:
      det = hh(1,1)*hh(2,2)-hh(2,1)*hh(1,2)
      if(cabs(det).eq.0) go to 60
      hhinv(1,1) = hh(2,2)/det
      hhinv(2,2) = hh(1,1)/det
      hhinv(1,2) = -hh(1,2)/det
      hhinv(2,1) = -hh(2,1)/det

c......... put hhinv into hh block
      s(1) = hhinv(1,1)
      s(2) = hhinv(2,1)
      s(3) = hhinv(2,2)

c>>>>>>>>>> compute TFs and put back into eh block
      call matmultc(hhinv,he,rt,2,2,nch2)
      do i =3,nch
         s(sindex(i,1)) = rt(1,i-2)
         s(sindex(i,2)) = rt(2,i-2)
      enddo

      call mmacb(he,hhinv,wrk,nch2,2,2)
      call matmultc(wrk,he,err,nch2,2,nch2)

      do i = 3,nch
         wrk(i) = s(sindex(i,i))
         do j = 3,i
            i1 = i-2
            j1 = j-2
            s(sindex(i,j)) = s(sindex(i,j))-conjg(err((j1-1)*nch2+i1))
            s(sindex(i,j)) = s(sindex(i,j))/n
         enddo
      enddo

      do i = 3,nch
         call mcohc(wrk(i),hh,he(1,i-2),rt(1,i-2),1,2,coh)
         s(sindex(i,i)) = cmplx(real(s(sindex(i,i))),coh)
      enddo
      return

60    do k=1,ns
         s(k)=(-999.,-999.)
      enddo
C70    continue
      return
      end
