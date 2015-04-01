      subroutine cn_wt(cjob,p,nt,nd,npc,nfb,xpc,u,ev,pw,v,wt,xv,
     &    var)
      include 'nstamx.inc'
      integer npcmx,nbmx,nsmx,ntmx2,nsmx2,nt,nd,npc,nfb,
     & nu2,nu,nv,nw,ntpc,nt2
      parameter (npcmx=5,nbmx=50,nsmx=(ntmx*(ntmx+1))/2,
     &       ntmx2=ntmx*2,nsmx2=nsmx*2)
      complex u(nt,npc),xpc(nd,npc),pw(nt,2),v(nt,npc),
     &       s0(ntmx,ntmx),s(nsmx),wrk(10,nsmx),xv(nd,npc),
     &        c(nsmx),work(nsmx,20),ut(ntmx,npcmx)
      real ev(nt),p,wt(nd),var(nt)
      character*1 cjob

c      write(91,*) 'In cn_wt'
c      write(91,*) 'p,nt,nd,npc,nfb',p,nt,nd,npc,nfb
c      write(91,*) 'pw'
c      write(91,'(10e12.4)') pw

ccc   express everyting in noise units
      do i = 1,nt
         temp = sqrt(var(i))
         do j = 1,npc
            u(i,j) = u(i,j)/temp
         enddo
         do j = 1,2
            pw(i,j) = pw(i,j)/temp 
         enddo
      enddo

ccc   copy u to ut (ut overwritten in lcpcls
      ntpc = nt*npc
      call ccopy(ntpc,u,1,ut,1)

ccc   change u back to physical units
      do i = 1,nt
         temp = sqrt(var(i))
         do j = 1,npc
            u(i,j) = u(i,j)*temp
         enddo
      enddo

ccc      construct "singularized SDM"
      call sing_s(nt,npc,ev,u,s)

ccc   sep_s_n takes "singularized sdm", pw source vectors
ccc    returns "cnoise" ...
      nu2 = npc-2
ccc      cjob = 'N'    ! signal and noise not correlated 
ccc      cjob = 'C'    ! signal and noise correlated
      call sep_s_n(cjob,s,nt,pw,2,v(1,3),nu2,c,work,s0)
c      write(91,*) 'After sep_s_n'
c      write(91,*) 'pw'
c      write(91,'(10e12.4)') pw
c      write(91,*) 'npc,CJOB = ',npc,cjob,' v'
c      write(91,'(10e12.4)') ((v(i,j),i=1,nt),j=3,npc)

ccc   copy pw into first part of array "v" Isecond part contains
ccc    coherent noise vectors already)
      nt2 = nt*2
      call ccopy(nt2,pw,1,v,1)

      nu = npc
      nv = npc
ccc     lcpcls finds coefficients relating pw-bart vectors to eigenvectors
      call lcpcls(ut,v,nu,nv,nt,c)

ccc   compute PW-CN coefficients
      call xu_xv(nv,nv,nd,c,xpc,xv)

ccc   find weights
      nw = nd/nfb
      call pwb_wts(xv,nv,nw,nfb,nd,wt,p)
c      write(91,*) 'nd,nv,nw,nfb,p,wt',nd,nv,nw,nfb,p
c      write(91,'(80i1)') (nint(wt(i)),i=1,nd)
       
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine pwb_wts(xv,nv,nw,nfb,nd,w,p)
 
      complex xv(nw,nfb,nv)
      real p,w(nd),pw_pwr
      integer nv,nw,nfb,nd,i,j,k,ik
 
      do i = 1,nw
         w(i) = 0.0
         do k = 1,nfb
            do j = 3,nv
               w(i) = w(i)+abs(xv(i,k,j))**2.
            enddo
         enddo
         pw_pwr = 0.0
         do k = 1,nfb
            pw_pwr = pw_pwr + (abs(xv(i,k,1))**2.+abs(xv(i,k,2))**2.)
         enddo
         w(i) = w(i)/pw_pwr
      enddo
c      write(91,'(80i2)') (nint(10*w(i)),i=1,nw)
      do i = 1,nw
         if(w(i).le.p) then
            do k = 1,nfb
               ik = i+(k-1)*nw
               w(ik) = 1.0
            enddo
         else
            do k = 1,nfb
               ik = i+(k-1)*nw
               w(ik) = 0.0
            enddo
         endif
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine lcpcls(u,v,nu,nv,nt,c)
      complex u(nt,nu),v(nt,nv),c(nv,nu),work(2000)
      character*1 side,trans,uplo,diag

ccc   compute [V^V]inv V^U = C [ use QR ]
      lwork = 100*nv
ccc       work(1:nv) = scalar factors of elementary reflectors (save)
ccc       work(np+1:np+lwork)
      iw1 = 1
      iw2 = 2*nv+1
      iw3 = iw2 + 2*lwork

c      write(91,*) 'nu,nv',nu,nv
c      write(91,*) 'u'
c      do i = 1,nt
c         write(91,*) (u(i,k),k=1,nu)
c      enddo
c      write(91,*) 'v'
c      do i = 1,nt
c         write(*,*) (v(i,k),k=1,nv)
c      enddo
ccc   Q-R factorization of V
      call cgeqrf(nt,nv,v,nt,work(iw1),work(iw2),lwork,info)

c         write(*,*) 'after Q-R info = ',info
 
ccc   compute Q^U (remember: Q is the full nt x nt orthogonal matrix

c         write(*,*) 'U before multiplying by Q'
c         do i = 1,nt
c            write(*,'(a1,10f8.4)') 'R',(real(u(i,j)),j=1,nu)
c            write(*,'(a1,10f8.4)') 'I',(aimag(u(i,j)),j=1,nu)
c         enddo

      side = 'L'
      trans = 'C'
      call cunmqr(side,trans,nt,nu,nv,v,nt,work(iw1),u,nt
     &    ,work(iw2),lwork,info)

c         write(*,*) 'info = ',info
c         write(*,*) 'U after multiplying by Q'
c         do i = 1,nt
c            write(*,'(a1,10f8.4)') 'R',(real(u(i,j)),j=1,nu)
c            write(*,'(a1,10f8.4)') 'I',(aimag(u(i,j)),j=1,nu)
c         enddo

ccc   move upper part of Q^U to C  , then compute R inv C
ccc    (LS solution)
      do i = 1,nv
         do j = 1,nu
            c(i,j) = u(i,j)
         enddo
      enddo
      trans = 'N'
      diag = 'N'
      uplo = 'U'
      call ctrtrs(uplo,trans,diag,nv,nu,v,nt,c,nv,info)

         write(*,*) 'info = ',info
         write(*,*) 'C'
         do i = 1,nv
            write(*,'(a1,10f8.4)') 'R',(real(c(i,j)),j=1,nu)
            write(*,'(a1,10f8.4)') 'I',(aimag(c(i,j)),j=1,nu)
         enddo

      return
      end
c______________________________________________________________________ 
c
      subroutine xu_xv(nv,nu,nd,c,xu,xv)
      complex xu(nd,nu),xv(nd,nv),c(nv,nu)
       write(*,*) 'in xu_xv; nd,nv,nu',nd,nv,nu 
       write(*,*) 'c',c
      do i = 1,nd
         do j = 1,nv
            xv(i,j) = (0.0,0.0)
            do k = 1,nu
               xv(i,j) = xv(i,j) + xu(i,k)*c(j,k)
            enddo
         enddo
      enddo
      return
      end
