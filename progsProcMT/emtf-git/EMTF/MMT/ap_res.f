ccc_____________________________________________________________________
ccc   General apparent resistivies and phases

      subroutine ap_res(z,sig_s,sig_e,period,nche,k0,nbt,ap)
      include 'nstamx.h'
      complex z(nchemx,2,*),sig_s(2,2,*),
     &  sig_e(nchemx,nchemx,nbt)
      real ap(nchemx,2,5,*),period(nbt),zerr
      integer k0,nche,nbt

      parameter (rad_deg=57.2958)

c       computes apparent resisitvities, phases, errors etc. from
c       the input arrays z, sig_s, and sig_e; z gives EMAP (MT) TFs
c       SIG_S is HH^(inverse) (or analogue for RR etc);
c       SIG_E is residual covariance matrix
c       transfer functions (including complex impedences);
c       output is given in array ap which contains for i=1,4:
c               (1)   app.resistivity
c               (2)   apparent resistivity error
c               (3)   phase of impedence (not of ap.res.!!!!)
c               (4)   phase error
c               (5)   multiple coherence

ccc   here nche is # of ELECTRIC channels (only ... Hz not included)
ccc      nch = ih(ista+1)-ih(ista)
ccc      nche = ih(ista+1)-ie(ista)
ccc   k0 + 1 is start of electric field channels for this stations
ccc      k0 = nch - nche - 2
      write(*,*) 'k0,nche',k0,nche
      do ib = 1,nbt
        do l = 1,2
          do k = 1,nche 
            k1 = k + k0
ccc         apparent resistivity
            ap(k,l,1,ib) =
     &      (abs(z(k1,l,ib))**2.)*period(ib)/5.
ccc         2 SE for apparent resistivity
            zerr = real(sig_e(k1,k1,ib)*sig_s(l,l,ib))
            ap(k,l,2,ib) = 2.*sqrt(zerr*period(ib)*
     &                                    ap(k,l,1,ib)*2/5.)
ccc         for plotting convenince make maximum 2se = rho !!!!!!!
            ap(k,l,2,ib) = min(ap(k,l,1,ib)-.01,
     &               ap(k,l,2,ib))
            ap(k,l,2,ib) = max(.01,ap(k,l,2,ib))

ccc         phase
            ap(k,l,3,ib) = rad_deg*
     &        atan(aimag(z(k1,l,ib))/real(z(k1,l,ib)))
ccc         2 SE for phase
            ap(k,l,4,ib) = 2.*sqrt(zerr/2.)*
     &                               rad_deg/abs(z(k1,l,ib))

ccc         Multiple coherence
ccc         **** NOT IMPLEMENTED ****
          enddo  ! k
        enddo  ! l
      enddo  ! ib
      return
      end
c______________________________________________________________________
c
      subroutine wrt_rho(nbt,period,ap,csta,arrayid,theta,ldf)

      include 'nstamx.h'

      integer nbt,fid,ldf(ngrpmx,*)
      real ap(nchemx,2,5,*),period(nbt),theta
      character*3 csta
      character*10 arrayid
      character*80 cfile

c      write(*,*) 
      fid = 41
      do length = 3,1,-1
         if(csta(length:length) .ne.' ') go to 5
      enddo
      length = 0
5     continue
      cfile = 'mt_'//csta(1:length)//'.'//arrayid
      write(*,*) 'file ',cfile,' open'
      open(unit=fid,file=cfile,status='unknown')
      write(fid,*) csta,arrayid
      write(fid,6005) theta
      write(fid,6003)
      do j = 1,nbt
         write(fid,6001) period(j),ldf(j),(ap(1,2,k,j),k=1,4),0.0,
     &   (ap(2,1,l,j),l=1,4),0.0
      enddo
      close(fid)
      return
6001  format(f10.3,i6,f8.2,f8.2,2f7.2,1x,f4.3,f8.2,f8.2,2f7.2,1x,f4.3)
6003  format('    per. ',' ndata',13x,'x - y',32x,'y - x'/'   (sec) ',
     1 6x,'   ap. res. ','  2 se ',' phase ',' 2 se','  coh ',
     2  '  ap. res. ','  2 se ',' phase ',' 2 se',' coh')
6005  format('rotation is ',f7.2,' degrees east of geomagnetic north')
      end
ccc_____________________________________________________________________
ccc
      subroutine wrt_emap(nbt,period,ap,csta,arrayid,orient,ldf,
     &    k0,nche)

      include 'nstamx.h'
      integer ldf(nbt),iex(nchemx),iey(nchemx)
      logical lex
      integer nbt,fid
      real period(nbt),ap(nchemx,2,5,*),theta,orient(2,*)
      character*3 csta
      character*10 arrayid
      character*80 cfile,ctitle(2,4)

      ctitle(1,1) = 'TM apparent resistivity'
      ctitle(1,2) = '2 SE TM apparent resistivity'
      ctitle(1,3) = 'TM phase'
      ctitle(1,4) =  '2 SE : TM phase'
      ctitle(2,1) = 'TE apparent resistivity'
      ctitle(2,2) = '2 SE TE apparent resistivity'
      ctitle(2,3) = 'TE phase'
      ctitle(2,4) =  '2 SE : TE phase'

ccc   computes mt interpretation parameters for EMAP configuration:
ccc   a bunch of E channels referenced to a single pair of H.
ccc   The first channel is assumed to be H oriented parallel to the span
ccc   Apparent resistivities corresponding to E fields with the same
ccc   orientation are referred to as "TM"; apparent resistivities for
ccc   orthogonal E fields are referred to as "TE"

      ney = 0
      nex = 0
      do i = k0+1,nche
         lex = (abs(orient(1,i+2)-orient(1,1)).lt.1.)
         if(lex) then
            nex = nex+1
            iex(nex) = i
         else
            ney = ney + 1
            iey(ney) = i
         endif
      enddo
            
      fid = 41
      do length = 3,1,-1
         if(csta(length:length) .ne.' ') go to 5
      enddo
      length = 0
5     continue
      cfile = 'em'//csta(1:length)//'.'//arrayid
      open(unit=fid,file=cfile,status='unknown')
      write(fid,*) csta,arrayid
      write(fid,*) nbt,' = # of periods'
      write(fid,*) nex,ney, '   = # of TM, TE dipoles'
      do l = 1,4
         write(fid,*) ctitle(1,l)
         do ib = 1,nbt
            write(fid,'(f12.5,i6,8f9.2)') period(ib),ldf(ib),
     &               (ap(iex(k),2,l,ib),k=1,nex)
         enddo
      enddo
      do l = 1,4
         write(fid,*) ctitle(2,l)
         do ib = 1,nbt
            write(fid,'(f12.5,i6,8f9.2)') period(ib),ldf(ib),
     &               (ap(iey(k),1,l,ib),k=1,ney)
         enddo
      enddo

      close(fid) 
 
      return
      end
ccc_____________________________________________________________________
ccc
ccc   output apparent resistivities and phases
      subroutine rho_out(z,sig_s,sig_e,period,nbt,nsta,ih,ie,
     &  orient,theta,csta,arrayid,ldf)
      include 'nstamx.h'
      complex z(nchemx,2,nstamx,*),sig_s(2,2,*),
     &  sig_e(nchemx,nchemx,nstamx,*),sig_s1(2,2,nbmax),
     &  z1(nchemx,2,nbmax),sig_e1(nchemx,nchemx,nbmax)
      real ap(nchemx,2,5,nbmax),period(nbt),orient(2,*),theta
      integer ldf(nbmax),ih(*),ie(*)
      character*10 arrayid
      character*3 csta(*)

      write(*,*) 'ih',(ih(k),k=1,nsta+1)
      write(*,*) 'ie',(ie(k),k=1,nsta)
      do ista = 1,nsta
         call z1sta(z,sig_e,sig_s,z1,sig_s1,sig_e1,ista,nbt,
     &                     nbmax,nstamx,nchemx)
        nch = ih(ista+1) - ih(ista)
        nche = ih(ista+1) - ie(ista)
        k0 = nch - nche - 2
        if(nche.eq.2) then
          write(*,*)  'rotating'
ccc       standard MT ... rotate to chosen coordinate system (theta)
ccc       output in standard MT format
          nchet = nche+k0
          call mt_rot(z1,sig_s1,sig_e1,theta,orient(1,ih(ista)),
     &      nchet,nbt)
          call ap_res(z1,sig_s1,sig_e1,period,nche,k0,nbt,ap)
          call wrt_rho(nbt,period,ap,csta(ista),arrayid,theta,ldf)
        else if((nche.gt.2) .and. (k0.ge.0)) then
ccc       EMAP ::  don't rotate, just output TM, TE apparaent resistivities
ccc       NOTE:  the direction of the first magnetic channel is taken to
ccc       define the strike of the EMAP line;  both modes are output however,
ccc       so this only effects accuracy of headings in output files
          write(*,*)  ' not rotating'
          call ap_res(z1,sig_s1,sig_e1,period,nche,k0,nbt,ap)
          write(*,*)  ' done with ap_res'
   
          call wrt_emap(nbt,period,ap,csta(ista),arrayid,
     &      orient(1,ih(ista)),ldf,k0,nche)
          write(*,*)  ' done with wrt_emap'
        endif
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine mt_rot(z,sig_s,sig_e,theta,orient,nche,nbt)
      include 'nstamx.h'
      real orient(2,*),u(nchmx,nchmx)
      complex z(nchemx,2,*),sig_s(2,2,*),sig_e(nchemx,nchemx,*)
     &  ,work(3,3)

ccc   MT_ROT takes TFs in Z for NCHE channels of output data
ccc   ( NCHE = 3 for MT impedance with tipper, = 2 without tipper)

      nch = nche+2
      call mkccmt(theta,u,orient,nch)
c      do k = 1,5
c         write(*,*) (u(k,l),l=1,5)
c      enddo

ccc   transform impedance
      do ib = 1,nbt
        do j = 1,nche
c          write(*,*) (z(j,k,ib),k=1,2)
          do k = 1,2
            work(j,k) = (0.,0.)
            do l = 1,nche
              do m = 1,2
                work(j,k) = work(j,k) + u(j+2,l+2)*z(l,m,ib)*u(k,m)
              enddo
            enddo
          enddo
        enddo
        do j = 1,nche
          do k = 1,2
            z(j,k,ib) = work(j,k)
          enddo
c          write(*,*) (z(j,k,ib),k=1,2)
        enddo

ccc     transform sig_s
        do j = 1,2
          do k = 1,2
            work(j,k) = (0.,0.)
            do l = 1,2
              do m = 1,2
                work(j,k) = work(j,k) + u(j,l)*sig_s(l,m,ib)*u(k,m)
              enddo
            enddo
          enddo
        enddo
        do j = 1,2
          do k = 1,2
            sig_s(j,k,ib) = work(j,k)
          enddo
        enddo

ccc     transform sig_e
        do j = 1,nche
          do k = 1,nche
            work(j,k) = (0.,0.)
            do l = 1,nche
              do m = 1,nche
                work(j,k) = work(j,k) + u(j+2,l+2)*sig_e(l,m,ib)*
     &            u(k+2,m+2)
              enddo
            enddo
          enddo
        enddo
        do j = 1,nche
          do k = 1,nche
            sig_e(j,k,ib) = work(j,k)
          enddo
        enddo
      enddo    !ib
      return
      end
