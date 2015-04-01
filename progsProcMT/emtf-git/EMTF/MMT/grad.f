      subroutine pwgslv(c,uvar,u,v,nt,nvec)
      parameter (pwfac=100.,nvmx=20,ndamp=30,nvnv=(nvmx*(nvmx+1))/2)
      complex c(6,nvec),u(nt,5),v(nt,nvec),bp(nvmx),cmiss(nvmx)
     &  , cpc(nvnv),cpcd(nvnv),bt(nvmx)
      real uvar(*),p(6,5),m(6,5),eta(6,5),rlim(2),Jv(0:ndamp),
     &  Jc(0:ndamp)
      data rlim/1.e-5,1.e+3/
      steplog=alog(rlim(2)/rlim(1))/ndamp

      print*,'uvar',(uvar(k),k=1,5)

      do ipwg = 1,5
         do j = 1,6
            p(j,ipwg) = 1.0
            m(j,ipwg) = 1.0
            eta(6,5) = 0.0
         enddo
         eta(ipwg,ipwg) = 1.0
      enddo

      p(1,1) = pwfac
      p(2,1) = pwfac
      p(1,2) = pwfac
      p(2,2) = pwfac
      p(6,1) = .2
      p(6,2) = .2
      p(6,3) = 0.
      p(6,4) = 0.
      p(6,5) = 0.

      do i=1,5
         m(1,i) = 0.
         m(2,i) = 0.
      enddo
      m(6,1) = 1.
      m(6,2) = 1.
      do i =3,5
         m(6,i) = 0.
      enddo
       
      do ipwg=1,5

ccc       form C^PC
         ii = 0
         do i = 1,nvec
            do j = 1,i
               ii = ii + 1
               cpc(ii) = (0.,0.)
               do k = 1,6
                  cpc(ii) = cpc(ii)+conjg(c(k,i))*p(k,ipwg)*c(k,j)
               enddo
            enddo
         enddo

ccc    form C^Pn
          do i = 1,nvec
          bt(i) = 0.0
             do k = 1,6
                bt(i) = bt(i) + conjg(c(k,i))*p(k,ipwg)*eta(k,ipwg)
             enddo
          enddo

ccc       scale for range of trade-off parameters tried
      scale = 0.0
         do k = 3,nvec
            scale= scale + uvar(k)
         enddo
         temp = 0.
         do i = 3,nvec
            ii = (i*(i+1))/2
            temp = temp + cpc(ii)
         enddo
         scale = temp/scale

ccc      solve system for ndamp+1 damping factors
         do idamp = 0,ndamp
            damp = rlim(1)*scale*(exp(steplog*idamp))

            ij = 0
            do i = 1,nvec
               bp(i) = bt(i)
               do j = 1,i
                  ij = ij + 1
                  cpcd(ij) = cpc(ij)
               enddo
               cpcd(ij) = cpcd(ij) + damp*uvar(i)
            enddo

            call cchdeci(cpcd,nvec,ier)
            call cltslv(cpcd,nvec,bp,1)
            call cutslv(cpcd,nvec,bp,1)
ccc          bp now contains soln vector for one damping parameter

ccc        calculate constraint misfit
            do i = 1,6
               cmiss(i) = eta(i,ipwg)
               do k = 1,nvec
                  cmiss(i) = cmiss(i) - c(i,k)*bp(k)
               enddo
            enddo
            Jc(idamp) = 0.0
            do i = 1,6
               Jc(idamp) = Jc(idamp) + m(i,ipwg)*abs(cmiss(i))**2.
            enddo

ccc         calculate error variance of linear combination
            Jv(idamp) = 0
            do k = 1,nvec
               Jv(idamp) = Jv(idamp) + uvar(k)*abs(bp(k))**2.
            enddo
         write(*,*) idamp,damp,Jv(idamp),Jc(idamp)
         if(mod(idamp,5).eq.0) write(*,'(5(2f7.4,2x))')(bp(k),k=1,5)
         enddo !  do idamp = 1,ndamp

ccc  chose damping parameter   (Let's just rty something; adjust later)
         do idamp = 0,ndamp
            if(100.*Jv(idamp).le.Jc(idamp) ) go to 20
         enddo
20       damp = rlim(1)*scale*(exp(steplog*idamp))
         write(*,*) idamp,damp,Jv(idamp),Jc(idamp)
ccc     and solve system for this damping parameter
         ij = 0
         do i = 1,nvec
            bp(i) = bt(i)
            do j = 1,i
               ij = ij + 1
               cpcd(ij) = cpc(ij)
            enddo
            cpcd(ij) = cpcd(ij) + damp*uvar(i)
         enddo
   
         call cchdeci(cpcd,nvec,ier)
         call cltslv(cpcd,nvec,bp,1)
         call cutslv(cpcd,nvec,bp,1)

c        form output vectors as linear combinations of eigenvectors
         write(*,'(5(2f7.4,2x))') (bp(k),k=1,5)
         do j = 1,nt
            v(j,ipwg) = (0.,0.)
            do k = 1,nvec
               v(j,ipwg) = v(j,ipwg) + u(j,k)*bp(k)
            enddo
         enddo
      enddo  ! do ipwg = 1,5

      return
      end
c___________________________________________________________
c
      subroutine mkuvar(ev,nt,pwg,nvec,evar,uvar,nf)
      real ev(*),pwg(11,nvec),uvar(nvec),evar

ccc      Crude estimate of variance
ccc        of each eigenvector ....

      uvar(1) = 1./(nf*(ev(1)/evar-1.))
      uvar(2) = 1./(nf*(ev(2)/evar-1.))
      do j = 3,nvec
         uvar(j) = (1.- pwg(2,j) )/(nt-5)
      enddo
c      do j = 1,nvec
c         uvar(j) = 1./(nf*abs(ev(j)/evar-1.))
c      enddo

      return
      end
c______________________________________________________________________
c
      subroutine pwgfit(u,w,c,nt,nvec,nh,ih,pwg)
      complex u(nt,nvec),w(nt,7),c(6,20),cdot,ctemp
      real pwg(11,20)
      integer ih(nh)

ccc     calculate dot product of eigenvectors [ U ] with
ccc     idealized gradient unit vectors [ W ] 
ccc        store in C, with PWG (j,i), j= 1,5 containing
ccc        an estimate of how well various idealized
ccc        terms are approximated by eigenvectors ...

c      print*,'IN PWGFIT'
c      print*,'nt,nvec,nh,ih',nt,nvec,nh,ih
c      write(*,*) 'u'
c      do i = 1,nt
c         write(*,'(10e12.4)') (u(i,j),j=1,5)
c      enddo
c      write(*,*) 'w'
c      do i = 1,nt
c         write(*,'(14e10.3)') (w(i,j),j=1,7)
c      enddo
      iounit = 0
      do i = 1,nvec
         do j = 1,5
            c(j,i) = cdot(w(1,j),u(1,i),nt)
            pwg(j+2,i) = real(c(j,i)*conjg(c(j,i)))
         enddo
         c(6,i) = cdot(w(1,7),u(1,i),nt)
      enddo
c      write(*,*) 'C'
c      do i = 1,6
c         write(*,'(12e10.3)') (c(i,j),j=1,5)
c      enddo

ccc       N-S gradients of H and E-W gradients of D
      do i = 1,nvec
         pwg(10,i) = abs((c(3,i)+c(4,i))/2.)**2.
         pwg(11,i) = abs((c(3,i)-c(4,i))/2.)**2.
      enddo

ccc       calculate magnitude in eigenvectors of
ccc        horizontal magnetic components [pwg(1,*)], and
ccc        vertical magnetics    [ pwg(8,*) ]
      do i =  1,nvec
         pwg(1,i) = 0.
         pwg(8,i) = 0.
         do j = 1,nh
            pwg(1,i) =pwg(1,i) + u(ih(j),i)*conjg(u(ih(j),i)) +
     &         u(ih(j)+1,i)*conjg(u(ih(j)+1,i))
            pwg(8,i) = pwg(8,i) + u(ih(j)+2,i)*conjg(u(ih(j)+2,i))
         enddo
      enddo

ccc    pwg(2,i) is an estimate of how well the horizontal components
ccc     of each eigenvector can be approximated as a linear combination
ccc      of the idealized source plane wave and gradient vectors
      do i = 1,nvec
         pwg(2,i) = 0.
         do k = 1,2
            pwg(2,i) = pwg(2,i) + pwg(k+2,i) + pwg(k+5,i)
         enddo
         ctemp = cdot(w(1,6),u(1,i),nt)
         pwg(2,i) = pwg(2,i) + real(ctemp*conjg(ctemp))
         ctemp = cdot(w(1,7),u(1,i),nt)
         pwg(9,i) = real(ctemp*conjg(ctemp))
         pwg(2,i) = pwg(2,i)/pwg(1,i)
         pwg(9,i) = pwg(9,i)/pwg(8,i)
      enddo

ccc     Why do this ????? .... so first two "plane wave" constraints
ccc     correpond to avg=normal
      do i = 1,nvec
         do j = 1,2
            c(j,i) = c(j,i)/sqrt(float(nh))
         enddo
      enddo

      return
      end
c______________________________________________________________________
c
      subroutine llkm(stcor,nsta,stkm,ier)
c       transform lat and long to kilometers using a simple projection

      integer nsta,ier
      real stcor(2,nsta),stkm(2,nsta),xb,yb,syb

      ier = 0

ccc     compute array center (in lat-lon)
      xb = 0.
      yb = 0.
      do  i=1,nsta
         if(stcor(1,i) .eq. 0.) then
            print*,'WARNING!!!!!! - station coordinates are (0.0,0.0)'
            ier = -1
            return
         end if
         xb = xb+stcor(1,i)
         yb = yb+stcor(2,i)
      enddo
      xb = xb/nsta
      yb = yb/nsta

      do i=1,nsta
         stkm(1,i) = (stcor(1,i)-xb)*111.04
         stkm(2,i) = (yb-stcor(2,i))*111.04*cos(stcor(1,i)*.01745)
      enddo

      syb = 0.
      do i = 1,nsta
         syb = syb + stkm(2,i)
      enddo
      syb = syb/nsta

      do i = 1,nsta
         stkm(2,i) = stkm(2,i) - syb
      enddo

      return
      end

c______________________________________________________________________
c
      subroutine mkw(stcc,nt,ih,nsta,w,theta,wrk)
      complex w(nt,7),wrk(*),cdot
      integer ih(nsta)
      real stcc(2,nsta),gstcc(2,100)

ccc     make idealized plane wave and gradient vectors
c       rotate station coordinates into new coordinate system
c          (with x-axis theta degrees east of N)
c         local coordinate system for each vector is the same
c         rotated coordinate system
ccc    Need as INPUT: station coordinates in array centered, km coord. s yst.

      c = cos(theta*3.14159/180.)
      s = sin(theta*3.14159/180.)

      do i = 1,nsta
         gstcc(1,i) = c*stcc(1,i)-s*stcc(2,i)
         gstcc(2,i) = s*stcc(1,i)+c*stcc(2,i)
      enddo

      do j = 1,nt
         do k = 1,7
            w(j,k) = (0.,0.)
         enddo
      enddo

c       Contents of array W:
c       1,2 : uniform source vectors (horiz. only)
c       3,4,5 : normalized gradient vectors  (horiz. only)
c       6 : by repalcing 3 with 6, the gradient vectors are orthonormal;
c       7 : vertical field (uniform Hz)
      do i = 1,nsta
         nch = ih(i+1) - ih(i)
         w(ih(i),1) = (1.,0.)
         w(ih(i)+1,2) = (1.,0.)
         w(ih(i),3)  = gstcc(1,i)
         w(ih(i)+1,3)  = 0.
c         w(ih(i)+1,3)  = gstcc(2,i)
c         w(ih(i),4)  = gstcc(1,i)
c         w(ih(i)+1,4)  = -gstcc(2,i)
         w(ih(i),4)  = 0.
         w(ih(i)+1,4)  = gstcc(2,i)
         w(ih(i),5)  = gstcc(2,i)
         w(ih(i)+1,5)  = gstcc(1,i)
         w(ih(i),6)  = gstcc(1,i)
         w(ih(i)+1,6)  = gstcc(2,i)
         if(nch.ne.4) w(ih(i)+2,7) = (1.,0.)
      enddo

c      orthonormalize gradient vectors
ccc     G-S orthogonalization .... 4 and 5 are already orthog., so this
ccc      just makes 6 orthog. to 4 and 5
      call cgrsch(w(1,4),nt,3,wrk)

c      Now normalize all
      do i = 1,7
         wsum = real(cdot(w(1,i),w(1,i),nt))
         wsum = sqrt(wsum)
         do j = 1,nt
            w(j,i) = w(j,i)/wsum
         enddo
      enddo

ccc    temporary test: print out
c      write(*,*) '!!!!! W !!!!!'
c      do i=1,7
c       write(*,*) (w(j,i),j=1,nt)
c      enddo
c
      return
      end
c______________________________________________________________________
c
      subroutine cgrsch(u,n,ip,wrk)

      complex u(n,ip),wrk(n),uv,cdot
      real uu
      integer i,ip,n,j,k

c     complex Gram-Schmidt orthogonalization of vectors in array u;
c      result overwrites u

      do i = 1,ip
         do k = 1,n
            wrk(k) = u(k,i)
         enddo
         do j = 1,i-1
            uv = cdot(u(1,j),u(1,i),n)
            do k = 1,n
               wrk(k) = wrk(k) - uv*u(k,j)
            enddo
         enddo
         uv = cdot(wrk,wrk,n)
         uu = sqrt(real(uv))
         do k = 1,n
            u(k,i) = wrk(k)/uu
         enddo
      enddo
      return
      end
