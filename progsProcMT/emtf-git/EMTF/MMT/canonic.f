c______________________________________________________________________
c
      subroutine lrcov(cjob,s,var,nt,nsta,ih,ab,r,wrk)

      include 'nstamx.inc'
      parameter (nsmx = ntmx*ntmx)
      complex s(*),wrk(*),s11(nsmx),s22(nsmx),s21(nsmx),ab(nt,nt),
     &    a(ntmx,nstamx),b(ntmx,nstamx)
      character*1 cjob
      real r(*),sdiag(ntmx),she(2),var(*),scale(ntmx)
      integer ih(*),i1(ntmx),i2(ntmx),sindex,ihe(ntmx)

c      write(*,*) 'in lrcov; cjob = ',cjob
      if(cjob.eq.'H') then
ccc      make a crude normalization between E & H fields (so that both are
ccc      of comparable magnitude
         call diag(s,nt,sdiag)
         sh = 0.0
         se = 0.0
         nh = 0
         ne = 0
         do ista = 1,nsta
            do k = ih(ista),ih(ista)+2
               sh = sh+sdiag(k)
               nh = nh+1
               ihe(k) = 1
            enddo
            do k = ih(ista)+3,ih(ista+1)-1
               se = se+sdiag(k)
               ne = ne+1
               ihe(k) = 2
            enddo
          enddo
          sh = sqrt(sh/nh)
          se = sqrt(se/ne)
          she(1) = 1.0
          she(2) = sh/se

          ij = 0
          do i = 1,nt
             do j = 1,i
                ij = ij + 1
                s(ij) = s(ij)*(she(ihe(i))*she(ihe(j)))
              enddo 
           enddo 
        else if((cjob.eq.'N').or.(cjob.eq.'n')) then
ccc     normalize to non-dimensional noise units before computing
ccc     canonical covariances
          do i = 1,nt 
             scale(i) = 1./sqrt(var(i))
          enddo
          ij = 0
          do i = 1,nt
             do j = 1,i
                ij = ij + 1
                s(ij) = s(ij)*scale(i)*scale(j)
              enddo 
           enddo 
        endif

ccc   outer loop over stations
      do ista = 1,nsta

ccc     compute index sets I1,I2
         n1 = 0
         n2 = 0
         do i =  1,nt
            if( (i.ge.ih(ista)).and.(i.lt.ih(ista+1)) ) then
               n2 = n2 + 1
               i2(n2) = i
            else
               n1 = n1 + 1
               i1(n1) = i
            end if
         enddo

ccc       reorder cross-products and partition matrix
         call mks1s2(s,i1,n1,i2,n2,s11,s21,s22)

         call canonic(cjob,n1,n2,s11,s22,s21,a,b,r(ih(ista)),wrk)
         call a_b_ab(a,b,n1,n2,nt,ab(1,ih(ista)))
      enddo    ! do ista = 1,nsta
      if((cjob.eq.'H').or.(cjob.eq.'h')) then
ccc      change S back to original scaling
         ij = 0
         do i = 1,nt
            do j = 1,i
               ij = ij + 1
               s(ij) = s(ij)/(she(ihe(i))*she(ihe(j)))
            enddo 
         enddo 
         do i = 1,nt
            r(i) = sqrt(abs(r(i)))
         enddo 
      else if((cjob.eq.'N').or.(cjob.eq.'n')) then
ccc      change S back to original scaling
         ij = 0
         do i = 1,nt
            do j = 1,i
               ij = ij + 1
               s(ij) = s(ij)/(scale(i)*scale(j))
            enddo 
         enddo 
         do i = 1,nt
            r(i) = sqrt(abs(r(i)))
         enddo 
      endif

      return
      end
c_______________________________________________________________________________
c
      subroutine canonic(cjob,n1,n2,s11,s22,s21,a,b,r,wrk)
      complex s11(*),s22(*),s21(n2,n1),wrk(*),a(n1,*),b(n2,*)
      include 'nstamx.inc'
      real r(*),cnorm,rwrk(ntmx)
      character*1 cjob

ccc   modified from cancor; will do:
ccc   canonical covariance (if cjob = 'V') or correlation (cjob='R')
      
      ns1 = (n1*(n1+1))/2
      ns2 = (n2*(n2+1))/2

      if((cjob.eq.'R').or.(cjob.eq.'r')) then
ccc      compute "normalized off-diagonal block" S11-h S12 S21-h
         call cchdec(s11,wrk,n1,ier)
ccc      bomb proofing for normalization by S22inv :
         call diag(s22,n1,rwrk)
ccc      find average
         sum = 0.0
         do i=1,n1
            sum = sum+rwrk(i)
         enddo
         sum = sum/(n1*10000.)
ccc      sum is now .0001 times avg diagonal element
ccc       add this to diagonal ... protects against null channels
         do i=1,n1
            ii = (i*(i+1))/2
            s22(ii) = s22(ii) + sum
         enddo
         call cchdec(s22,wrk(ns1+1),n2,ier)
         do i = 1,n2
            do j = 1,n1
               a(j,i) = conjg(s21(i,j))
            enddo
         enddo
         call cltslv(wrk,n1,a,n2)
         do i = 1,n2
            do j = 1,n1
               s21(i,j) = conjg(a(j,i))
            enddo
         enddo
         call cltslv(wrk(ns1+1),n2,s21,n1)
      endif
 
ccc   Form S12 S12^ ("normalized" or otherwise)  
      ii = 0
      do i = 1,n1
         do j = 1,i
            ii = ii + 1
            s11(ii) = (0.0,0.0)
            do k = 1,n2
               s11(ii) = s11(ii) + conjg(s21(k,i))*s21(k,j)
            enddo
         enddo
      enddo
 
ccc   Form S12^ S12 ("normalized" or otherwise)  
      ii = 0
      do i = 1,n2
         do j = 1,i
            ii = ii + 1
            s22(ii) = (0.0,0.0)
            do k = 1,n1
               s22(ii) = s22(ii) + s21(i,k)*conjg(s21(j,k))
            enddo
         enddo
      enddo

ccc   Do eigenvector decomposition
      nm = min(n1,n2)
      nt = ns1+ns2+1
      call rsp(s11,n1,nm,a,r,0,wrk(nt),dum)
      call rsp(s22,n2,nm,b,r,0,wrk(nt),dum)

      
      if((cjob.eq.'R').or.(cjob.eq.'r')) then
         call cutslv(wrk,n2,a,n1)
         call cutslv(wrk(ns1+1),n2,b,n1)

ccc      renormalize a and b
         do i = 1,nm
            temp = cnorm(a(1,i),n1)
            do j = 1,n1
               a(j,i) = a(j,i)/temp
            enddo
            temp = cnorm(b(1,i),n2)
           do j = 1,n2
              b(j,i) = b(j,i)/temp
           enddo
        enddo
      endif

      return
      end
c______________________________________________________________________________
c
      subroutine ltcmlt(s,a,n,m)
      complex s(*),a(n,m)
      integer sindex
 
      do j = 1,m
         do i = 1,n
            a(i,j) = a(i,j)*s(sindex(i,i))
            do k = i+1,n
               a(i,j) = a(i,j) + a(k,j)*conjg(s(sindex(i,k)))
            enddo
         enddo
      enddo
      return
      end
c______________________________________________________________________________
c
      function cnorm(x,n)
      complex x(n)
      real cnorm
      integer n,i

      cnorm = 0.0
      do i = 1,n
         cnorm = cnorm + x(i)*conjg(x(i))
      enddo
      cnorm = sqrt(cnorm)
      return
      end
ccc______________________________________________________________________
c
      subroutine mks1s2(s,i1,n1,i2,n2,s11,s21,s22)
      complex s(*),s11(*),s22(*),s21(n2,n1)
      integer sindex,i1(n1),i2(n2)
      ii = 0
      do i = 1,n1
         do j = 1,i
            ii = ii + 1
            s11(ii) = s(sindex(i1(i),i1(j)))
         enddo
      enddo

      ii = 0
      do i = 1,n2
         do j = 1,i
            ii = ii + 1
            s22(ii) = s(sindex(i2(i),i2(j)))
         enddo
      enddo

      do i = 1,n2
         do j = 1,n1
            if(i2(i).ge.i1(j)) then
               s21(i,j) = s(sindex(i2(i),i1(j)))
            else
               s21(i,j) = conjg(s(sindex(i1(j),i2(i))))
            end if
         enddo
      enddo
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine a_b_ab(a,b,n1,n2,nt,ab)
      complex a(n1,n1),b(n2,n1),ab(nt,n1)

      do j = 1,n2
         do i = 1,n2
            ab(i,j) = b(i,j)
         enddo 
         do i = 1,n1
            ab(n2+i,j) = a(i,j)
         enddo 
      enddo
      return
      end
