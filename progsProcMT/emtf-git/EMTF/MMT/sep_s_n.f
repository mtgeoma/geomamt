      subroutine sep_s_n(cjob,s,n,u1,nu1,u2,nu2,c,work,s0)
      include 'nstamx.inc'
      complex s(*),work(n,n,*),u1(n,nu1),u2(n,nu2)
     &  ,s0(n,n),c(nu2,nu1)
      real eval(ntmx)
      character*1 cjob

ccc   orthonormalize u1  (current best guess for "plane wave response vectors")
      call cgrsch(u1,n,nu1,work)

ccc   project SDM (onto orthogonal complement of u1) 
      ij = 0
      do i = 1,n
         do j = 1,i-1
            ij = ij + 1
            work(i,j,1) = s(ij)
            work(j,i,1) = conjg(s(ij))
         enddo
         ij = ij + 1
         work(i,i,1) = s(ij)
      enddo

ccc   temporary for debuging
      call rsp(s,n,nu2,u2,eval,0,work(1,1,2),dum)
c      write(*,*) 'raw eigenvalues of singularized SDM'
c      write(*,*) (eval(k),k=1,n)
ccccccccc

      ns = (n*(n+1))/2
      call s_proj(work,s0,n,u1,nu1,work(1,1,2))

      ij = 0
      do i = 1,n
         do j = 1,i
            ij = ij + 1
            s(ij) = s0(i,j)
         enddo
      enddo

ccc   get domiant eigenvectors from projected SDM
      call rsp(s,n,nu2,u2,eval,0,work(1,1,2),dum)
ccccccccc
c      write(*,*) 'residual eigenvalues of projected, singularized SDM'
c      write(*,*) (eval(k),k=1,n)
      call ccmult(u1,u1,work(1,1,2),nu1,n,nu1)
c      write(*,*) 'U1^U1'
c      write(*,*) (work(k,1,2),k=1,4)
      call ccmult(u2,u2,work(1,1,2),nu2,n,nu2)
c      write(*,*) 'U2^U2'
c      write(*,*) (work(k,1,2),k=1,9)
      call ccmult(u2,u1,work(1,1,2),nu2,n,nu1)
c      write(*,*) 'U2^U1'
c      write(*,*) (work(k,1,2),k=1,6)
ccccccccc
ccc    (how is n2 chosen ???)
      if((cjob.eq.'C').or.(cjob.eq.'c')) then
ccc      not assuming secondary sources are uncorrelated with PW
ccc      just use orthogonal complement of U1 for non-PW source estimate
          write(*,*) 'assuming coherent noise correlated with PW'
         return
      endif
ccc   else assume that secondary sources are uncorrelated; can then
cc     refine non-PW source estimate
ccc 
ccc   U2^S 
      call ccmult(u2,work,work(1,1,2),nu2,n,n)
ccc   U2^SU1
      call cmult(work(1,1,2),u1,c,nu2,n,nu1)
c         write(*,*) 'U2^SU1'
c         do i = 1,nu2
c            write(*,'(2e12.4,2x,2e12.4)') (c(i,j),j=1,nu1)
c         enddo
ccc   U2^SU2
      call cmult(work(1,1,2),u2,work,nu2,n,nu2)
c         write(*,*) 'U2^SU2'
c         write(*,'(2e12.4,2x,2e12.4,2x,2e12.4)') (work(k,1,1),k=1,9)
c
      call sym_stor(work,s,nu2)
         
ccc   compute U2^SU2 inv U2^SU1
      call cchdeci(s,nu2,ier)
      call cltslv(s,nu2,c,nu1)
      call cutslv(s,nu2,c,nu1)
c         write(*,*) 'Coefficients C'
c         do i = 1,nu2
c            write(*,'(2e12.4,2x,2e12.4)') (c(i,j),j=1,nu1)
c         enddo

ccc   compute U2 = U2 + U1C^
      call cmultc(u1,c,work,n,nu1,nu2)
      c1 = (1.,0.)
      call cmadd(u2,c1,work,n,nu2)
ccc   orthonormalize u2
      call cgrsch(u2,n,nu2,work)

ccc   compute covariances for U1, U2
ccc   (refinements for later ...)

      return
      end
c______________________________________________________________________
c
      subroutine s_proj(s,s0,n,u,m,work)
ccc   s is nxn hermitian SDM, regular storage mode
ccc   u is nxm orthonormal matrix
ccc   output is  s0 = [ I - UU^ ] S [ I - UU^ ] (in S)
      complex s(n,n),u(n,m),work(n,n,2),c1m,c1,s0(n,n)

c      write(*,*) 's'
ccc   copy s to s0
      do i = 1,n
         do j = 1,n
            s0(i,j) = s(i,j) 
         enddo
c         write(*,'(15e8.2)') (real(s0(i,j)),j=1,n)
c         write(*,'(15e8.2)') (aimag(s0(i,j)),j=1,n)
      enddo

ccc   form UU^S
      c1m = (-1.,0.)
      call ccmult(u,s,work,m,n,n)
      call cmult(u,work,work(1,1,2),n,m,n)

c      write(*,*) 'UU^S'
c      do i = 1,n
c         write(*,'(15e8.2)') (real(work(i,j,2)),j=1,n)
c         write(*,'(15e8.2)') (aimag(work(i,j,2)),j=1,n)
c      enddo

ccc   form S-UU^S-SUU^
      call cmadd(s0,c1m,work(1,1,2),n,n)
      call cmaddc(s0,c1m,work(1,1,2),n,n)

c      write(*,*) 'S-UU^S-SUU^'
c      do i = 1,n
c         write(*,'(15e8.2)') (real(s0(i,j)),j=1,n)
c         write(*,'(15e8.2)') (aimag(s0(i,j)),j=1,n)
c      enddo

ccc   form UU^SUU^
      call cmult(work(1,1,2),u,work,n,n,m)
      call cmultc(work,u,work(1,1,2),n,m,n)
      c1 = (1.,0.)
      call cmadd(s0,c1,work(1,1,2),n,n)

c      write(*,*) 'S-UU^S-SUU^+UU^SU^U'
c      do i = 1,n
c         write(*,'(15e8.2)') (real(s0(i,j)),j=1,n)
c         write(*,'(15e8.2)') (aimag(s0(i,j)),j=1,n)
c      enddo

      return
      end 
c______________________________________________________________________
c
      subroutine cmadd(a,c,b,n,m)
ccc   form sum a+c*b, c a scalar, a and b nxm.  Return in a.
      complex a(n,m),b(n,m),c
      do i = 1,n
         do j = 1,m
            a(i,j) = a(i,j) + c*b(i,j)
         enddo
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine cmaddc(a,c,b,n,m)
ccc   form sum a+c*b^, c a scalar, a and b nxm.  Return in a.
      complex a(n,m),b(m,n),c
      do i = 1,n
         do j = 1,m
            a(i,j) = a(i,j) + c*conjg(b(j,i))
         enddo
      enddo
      return
      end
c______________________________________________________________________
c
        subroutine cmult(a,b,c,nn,np,nq)

c       multiplies n x p matrix a by p x q matrix b and puts result in c

        complex a(nn,np),b(np,nq),c(nn,nq),temp

        do 10 in=1,nn
        do 10 iq=1,nq
        temp= (0.,0.)

           do  5 ip=1,np
5          temp=temp+a(in,ip)*b(ip,iq)

10      c(in,iq)=temp
        return
        end
c______________________________________________________________________
c
      subroutine ccmult(a,b,c,nn,np,nq)
ccc   multiplies conjugate  transpose of n x p matrix a by p x q matrix b and puts result in c
      complex a(np,nn),b(np,nq),c(nn,nq),temp
      do in=1,nn
         do iq=1,nq
         temp= (0.,0.)
            do ip=1,np
               temp=temp+conjg(a(ip,in))*b(ip,iq)
            enddo
            c(in,iq)=temp
         enddo
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine cmultc(a,b,c,nn,np,nq)
ccc   multiplies n x p matrix a by conjugate transpose of q x p matrix b and puts result in c
      complex a(nn,np),b(nq,np),c(nn,nq),temp
      do in=1,nn
         do iq=1,nq
            temp= (0.,0.)
            do ip=1,np
               temp=temp+a(in,ip)*conjg(b(iq,ip))
            enddo
            c(in,iq)=temp
         enddo
      enddo
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine sc_u(var,u,nt,nu,ijob)

ccc   if ijob = 1 multiply by sqrt(var)
ccc   if ijopb = -1 divide
      complex u(nt,nu)
      real var(nt)
      do i = 1,nt
         if(ijob.eq.1) then
            scale = sqrt(var(i))
         else
            scale = 1./sqrt(var(i))
         endif
         do j = 1,nu
            u(i,j) = u(i,j)*scale
         enddo
      enddo
      return
      end
ccc_____________________________________________________________________
c
      subroutine sym_stor(work,s,nu2)
      complex work(nu2,nu2),s(*)

      ij = 0
      do i = 1,nu2
         do j = 1,i
            ij = ij + 1
            s(ij) = work(i,j)
         enddo
      enddo
      return
      end
