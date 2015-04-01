      subroutine rsp(s,n,p,u,ev,icon,wrk,cov)

c    subroutine returns estimate of response space of dimension p
c    s is input SDM (n x n); resposne space is spanned by columns
c    of u; eigenvalues are returned in ev; icon is an integer variable
c    which is =0 for isotrpic errors; =1 for diagonal covariance
c    =2 for general covariance;  parameters for error covariance are
c    stored in cov

c           this version does not return orthogonal response space vectors
ccc    ALSO:  not normalized any more (evecs are in physical units of E and H)
      integer p,icon,n,i,j,nj1,ierr
      complex s(*),u(n,p)
      real cov(*),ev(n),wrk(n,n,6)

      if( icon .le. 1) then
         call mkarai(s,n,wrk(1,1,1),wrk(1,1,2))

         if(icon .eq.1) then
            do 3 i = 1,n
3           wrk(i,1,5) = sqrt(cov(i))

            do 5 i = 1,n
            do 5 j = 1,n
            wrk(i,j,1) = wrk(i,j,1)/(wrk(i,1,5)*wrk(j,1,5))
5           wrk(i,j,2) = wrk(i,j,2)/(wrk(i,1,5)*wrk(j,1,5))
         end if
      else
c       compute transformation matrix w1 = (cov)**-.5
         call cchdeci(cov,n,ierr)
         if(ierr.eq.0) then
c             transform s to s' = w1 s w1^
            call mkicf(wrk,n)
            call cltslv(cov,n,wrk,n)
            call usuc1(s,wrk,n,n,wrk(1,1,3),wrk(1,1,5))
            call mkarai(wrk(1,1,5),n,wrk(1,1,1),wrk(1,1,2))
         else
            write(*,*) 'error in cholesky decomposition',
     &             'error covariance is numerically singular'
            call mkarai(s,n,wrk(1,1,1),wrk(1,1,2))
         end if
      end if

      call cheig(n,wrk(1,1,1),wrk(1,1,2),wrk(1,2,5),wrk(1,3,5),
     &     wrk(1,1,3),wrk(1,1,4),ev,ierr)
      if(ierr.ne.0) print*,'error in svd'
c      write(*,*) 'ev',ev

      do i = 1,n/2
         temp = ev(i)
         ev(i) = ev(n-i+1)
         ev(n-i+1) = temp
      enddo

      if(icon .le. 1) then
         do 10 i = 1,n
         do 10 j = 1,p
         nj1 = n-j+1
         u(i,j) = cmplx(wrk(i,nj1,3),wrk(i,nj1,4))
10       continue

         if(icon.eq.1) then
            do 15 i = 1,n
            do 15  j = 1,p
15          u(i,j) = u(i,j)*wrk(i,1,5)
         end if
      else
         call pwrsptrn(wrk(1,1,3),wrk(1,1,4),n,p,cov,u)
      end if

ccc   output vectors are not normalized any more

      return
      end 
c_______________________________________________________________
c
        subroutine pwrsptrn(zr,zi,n,neig,t1,pwrsp)
 
        real zr(n,n),zi(n,n)
        complex pwrsp(n,neig),t1(*),temp
 
c       gets neig dominant eigenvectors out of arrays zr,zi and pre-
c       multiplies by t1; puts result in pwrsp

        kk = 0
        do 20 j = 1,neig
        jj = n-j+1
        ii = 0
           do 20 i = 1,n
           temp = (0.,0.)
              do 10 k =1,i
              ii = ii + 1
10            temp = temp + t1(ii) * cmplx(zr(k,jj),zi(k,jj))
           pwrsp(i,j) = temp
20      continue
        return
        end
c______________________________________________________________________
c
        subroutine mkicf(s,n)
 
c       makes an n x n complex identity matrix in full storage mode
 
        complex s(n,n)
 
        do 5 i = 1,n
        do 5 j = 1,n
5       s(i,j) = (0.,0.)
 
        do 10 i = 1,n
10      s(i,i) = (1.,0.)
 
        return
        end
