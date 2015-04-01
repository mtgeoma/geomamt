      subroutine ediazimuth( ind, azimuth )
c
c     gets the azimuth from an EDI file (must be open already)

      implicit none

      integer ind

      character*80 line
      
      real azimuth, azim, azim2
      
      include 'ctrlblk.inc'

      rewind( ind )
c...look for ZROT in impedance file
      rewind( ind )
      read(ind,'(a)') line
      do while( index(line,'ZROT').eq.0 )
        read(ind,'(a)',end=11) line
      enddo
      read(ind,*) azim, azim2
      if( iprint.ge.2 ) then
        write(*,*) 'edicopyinfo: azims from impedance edi =>', 
     &             azim, azim2
      endif
      if( azim.ne.azim2 ) then
        write(*,*)'Different azimuths for different frequencies'
	write(*,*)'STRIKE cannot handle this case'
	stop
      endif
      goto 12
      
c...ZROT not found, must be an impedance file, look for ROTSPEC entry in SPECTRA line
11    rewind( ind )
      read(ind,'(a)') line
      do while( index(line,'>SPECTRA ').eq.0 )
        read(ind,'(a)',end=9999) line
      enddo
      read(line(index(line,'ROTSPEC')+8:),'(i4)') azim
c...find second azimuth
      do while( index(line,'>SPECTRA ').eq.0 )
        read(ind,'(a)',end=9999) line
      enddo
      read(line(index(line,'ROTSPEC')+8:),'(i4)') azim2
      if( iprint.ge.2 ) then
        write(*,*) 'edicopyinfo: azims from spectra edi =>', azim, azim2
      endif
      if( azim.ne.azim2 ) then
        write(*,*)'Different azimuths for different frequencies'
	write(*,*)'STRIKE cannot handle this case'
	stop
      endif
      
12    azimuth = azim
      if( iprint.ge.1 ) then
        write(*,*) 'ediazimuth: azimuth =', azimuth
      endif

      return
      
9999  write(*,*)'edicopyinfo: Could not find either ZROT or ROTSPEC'//
     &          ' in EDI file'
      return
      end
