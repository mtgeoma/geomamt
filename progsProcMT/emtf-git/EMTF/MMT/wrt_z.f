      subroutine wrt_z(z,sig_s,sig_e,nbt,ngrp,igrp,cgrp,period,stcor,
     &   decl,orient,chid_grp,sta,sta_grp,level,ibandlim,arrayid,chead,
     &   ldf,dr)
     
      include 'nstamx.inc'

      complex z(nchemx,2,nstamx,*),sig_s(2,2,nstamx,*),
     &  sig_e(nchemx,nchemx,nstamx,*)
      real period(nbt),stcor(2,*),decl(*),orient(2,*),dr(*)
      character*20 cgrp(ngrp)
      integer sta_grp(nchemx,*)
      character*6 chid_grp(nchemx,*)
      character*3 sta(*)
      character*10 arrayid
      character*80 cfile,chead
      integer nbt,level(nbt),ibandlim(2,nbt),ldf(ngrpmx,nbt),fid,ngrp
     &   ,igrp(ntmx+1,*),lar_id,iclong

      fid = 51
      do k = 1,ngrp
         nch = igrp(1,k)
         if(nch .ge.3) then
ccc         output TF file ...
            nche = nch-2
            do length = 20,1,-1
               if(cgrp(k)(length:length) .ne.' ') go to 5
            enddo
         length = 0
5        continue
         lar_id = iclong(arrayid,10)
         cfile = cgrp(k)(1:length)//'_'//arrayid(1:lar_id)//'.zmm'
         open(unit=fid,file=cfile,status='unknown')

c...  write header information
         ista = sta_grp(1,k)
         write(fid,*) 'TRANSFER FUNCTIONS IN MEASUREMENT COORDINATES'
         write(fid,*) '********** WITH FULL ERROR COVARINCE*********'
         write(fid,'(a80)') chead
         write(fid,'(a20)') cgrp(k)
         write(fid,105) stcor(1,ista),stcor(2,ista),decl(ista)
         write(fid,110) nch,nbt
         write(fid,*) 'orientations and tilts of each channel '
         do l=2,igrp(1,k)+1
            ll = igrp(l,k)
            write(fid,115) ll,orient(1,ll),orient(2,ll),
     &       sta(sta_grp(l-1,k)),chid_grp(l-1,k)
         enddo
         write(fid,*)

         do ib = 1,nbt
            write(fid,120) period(ib),level(ib),ibandlim(1,ib),
     &           ibandlim(2,ib)
            write(fid,125) ldf(k,ib),dr(1)
c...        TF (impedance + tipper ... emap dipole etc)
            write(fid,*) 'Transfer Functions'
            do ii = 1,nche
               write(fid,140) (z(ii,l,k,ib),l=1,2) 
            enddo
c...        SIGMA_S^-1 (inverse coherent signal power matrix)
            write(fid,*) 'Inverse Coherent Signal Power Matrix'
            write(fid,140) sig_s(1,1,k,ib)
            write(fid,140) (sig_s(2,ii,k,ib),ii=1,2)
c...        RESIDUAL COVARIANCE
            write(fid,*) 'Residual Covaraince'
            do ii = 1,nche
               write(fid,140) (sig_e(ii,l,k,ib),l=1,ii)
            enddo
         enddo
         close(fid)
         endif
      enddo
c100   format('station    : ',a20)
105   format('coordinate ',f10.5,1x,f10.5,1x,'declination ',f8.2)
110   format('number of channels ',i3,2x,' number of frequencies',i4)
115   format(i5,1x,f8.2,1x,f8.2,1x,a7,2x,a6)
120   format('period : ',1pe12.6,3x,' decimation level ',i3,3x,
     +       ' freq. band from ',i6,' to ',i6)
125   format('number of data point ',i6,' sampling freq. ',1pe12.6,
     +     ' Hz')
140   format(16es12.4)
      end
ccc_____________________________________________________________________
ccc
      subroutine z1sta(z,sig_e,sig_s,z1,sig_s1,sig_e1,ista,nbt,
     &                     nbmax,nstamx,nchemx)

      complex z(nchemx,2,nstamx,nbmax),sig_s(2,2,nstamx,nbmax),
     &  sig_e(nchemx,nchemx,nstamx,nbmax),z1(nchemx,2,nbmax),
     &  sig_s1(2,2,nbmax),sig_e1(nchemx,nchemx,nbmax)

      do ib = 1,nbt
        do l = 1,2
          do k = 1,2
            sig_s1(k,l,ib) = sig_s(k,l,ista,ib)
          enddo
        enddo

        do k = 1,nchemx
          do l = 1,2
            z1(k,l,ib) = z(k,l,ista,ib)
          enddo
          do l = 1,nchemx
            sig_e1(k,l,ib) = sig_e(k,l,ista,ib)
          enddo
        enddo
      enddo

      return
      end
