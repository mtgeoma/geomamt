       subroutine sdmwrite(nbt,ns,name,cdirout,nsmx,zt,fb,rdf,
     &          nch, chid, orient,bw,ntmx,stcor,decl)

c     writes the spectral density matrix to a file. Format is the intermediate 
c     GEOTOOLS format readable with the routine EDIWRITE to transform into an
c     EDI file. SDM is rotated into a geographic north coordinate system.
c     SDMWRITE only writes 
        character*40    name
        character*30    cdirout
        character*80    cfile
        character*6     chid(*)    
        
        integer         SPEC_UNIT,ERR
        integer         stcor1dms(3),stcor2dms(3)
        integer         ipair(2,3)
        integer         npair, nopair

        complex         zt(nsmx,*), zz(28)
        real            fb(nbt),rdf(ntmx,nbt)
        real            bw(*)
        real            orient(2,*)
        real            xx(20,20)
        real            stcor(2), decl
        real            u(7,7)
        real            orient1(2)
        real            cc(2,2)
        PARAMETER       (SPEC_UNIT = 28)

        if (nch .ne. 5 .and. nch .ne. 7 ) then
           write(*,'(''This number of channels is not supported for '',/,
     &          ''SDM output ( NCH = '',i3,'')'')') nch
           return
        endif


        ll = iclong(cdirout,30)
        ll2 = iclong(name,45)
        if (name(ll2-4:ll2) .eq. '.zss' .or.
     &       name(ll2-4:ll2) .eq. '.zrr') then
           ll2 = ll2-4
        endif
        cfile = cdirout(1:ll)//'/'//name(1:ll2)//'.sdm'

        open (unit = SPEC_UNIT, file = cfile)           

        theta = 0.0
        npair = nch/2
        ipair(1,1) = 1
        ipair(2,1) = 2
        nopair = 3
        ipair(1,2) = 4
        ipair(2,2) = 5
        if (npair .gt. 2) then
           ipair(1,3) = 6
           ipair(2,3) = 7
        endif
        do i = 1,nch
           orient(1,i) = orient(1,i) + decl
        enddo

c     first set up the rotation matrix from the angles of all paired channels
c     this needs to be done only ones for the whole set of frequencies, so lets do
c     it here
c     fill rotation matrix with zeros

      do k = 1,nch
         do l = 1,nch
            u(l,k) = 0.0
         enddo
      enddo

c     if there are HZs make the corresponding diagonal element 1
      if (nopair .gt. 0)  u(nopair,nopair) = 1.0

      do ip = 1, npair
c     get rotation matrices for pairs
         orient1(1) = orient(1,ipair(1,ip))
         orient1(2) = orient(1,ipair(2,ip))
         call ccmat(cc, orient1,theta)
c     map the pair rotation matrices into the full rotation matrix
         do k = 1,2
            do l = 1,2
               u(ipair(k,ip),ipair(l,ip)) = cc(k,l)
            enddo
         enddo
      enddo
c      write(*,'(7f7.3)') ((u(j,k), k = 1,nch),j = 1,nch)
C-EDI     ----WRITE 3 LINE HEADER, WITH SECTID,NFREQS,AND SPEC_ROT
      WRITE(SPEC_UNIT,'(A,2X,A)',IOSTAT=ERR) 'SITE',name(1:ll2)
      WRITE(SPEC_UNIT,'(A,2X,I4)',IOSTAT=ERR) 'NFREQ', nbt
      WRITE(SPEC_UNIT,'(A,A)',IOSTAT=ERR) 'SPEC_ROT','  0'


      DO I=1,nbt

           do j = 1,ns
             zz(j) = zt(j,i)
           enddo

c         print*,' VOR ROTSDM', (zz(k),k=1,ns)
c        call ccsdm(zz, nch, ipair, npair, orient, theta)
         call rotsdm(zz, u, nch)
c         print*,' NACH ROTSDM', (zz(k),k=1,ns)
         kj = 0
         do k = 1, nch
            do j = 1, k
               kj = kj + 1
c               xx(j,k) = real(zz(kj))
c               if (j.ne.k) xx(k,j) = aimag(zz(kj))
               xx(k,j) = real(zz(kj))
               if (j.ne.k) xx(j,k) = -1. * aimag(zz(kj))
            enddo
         enddo
         
c         if (bw(i) .eq .0.0) 
         bw(i) = 0.3*(1./fb(i))       
         
 110     FORMAT(A,1PE10.4,1X,1PE10.4,1X,I5,1X,I5)
         WRITE(SPEC_UNIT,110,IOSTAT=ERR) ' FREQ/BW/STACKS/HARMS ',
     1        1./fb(i),bw(i),int(rdf(3,i)),1
         
 120     FORMAT( 5(1X,1PE15.8))
         WRITE(SPEC_UNIT,120,IOSTAT=ERR) ((xx(J,K),K=1,nch),
     1        J=1,nch)
         
         IF (ERR.NE.0) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'ERROR: Error writing Spectra file data. ',ERR
            GOTO 99
         END IF
      enddo 

      close (SPEC_UNIT)
      call dec2dms(stcor(2),stcor2dms)
      call dec2dms(stcor(1),stcor1dms)
      num_ref = 0
      if (nch.gt.5) num_ref = -1

      print*
      print*,' Crosspower matrix written to ', cfile(1:iclong(cfile,80))

      write(*,199) name(1:3), name(1:ll2),stcor1dms,stcor2dms,
     &     0, 100, 0, 0, num_ref, 1, 999999
 199  format(/,'The follolwing line is a template for the SITES file',/,
     &         'needed by EDIWRITE (GEOTOOLS). Paste it into the ',/,
     &         'SITES file and edit the 3 char site name to a number',/,
     &         'and move the sdm file into the spectra directory as',/,
     &         'a file with this number as the file name.',//,
     &     a3,': ',a8,', ',i3.2,':',i2.2,':',i2.2,', ',i3.2,':',i2.2,
     &     ':',i2.2,',',i3,',',i7,',',i5,',',i5,',',i4,',',i4,',',i8.6)

 99   close (SPEC_UNIT)
      return
      end

c***********************************************************
      subroutine dec2dms(orig,dms)
      real      degs, orig
      integer   dms(3)

      degs = orig
      sig = sign(1.,degs)
      dms(1) = int(abs(degs))
      degs = (abs(degs) - float(dms(1)))*3600.
      dms(1) = sig * dms(1)
      dms(2) = int(degs/60.)
      dms(3) = int(mod(degs,60.))
      return
      end
