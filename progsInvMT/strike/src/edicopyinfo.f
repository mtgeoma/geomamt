      subroutine edicopyinfo( ind, outd )

      implicit none

      integer ind, outd

      character*80 line

      integer leq, lat, latmin, latsec, lon, lonmin, lonsec, llon, lelev

      real xlat, xlon, elev, azim, azim2
      
      include 'ctrlblk.inc'

c...gets inforamtion from EDI input and write to output file 
c   (both must be open already)

      if( iprint.ge.2 ) write(*,*)'Entered edicopyinfo'

c...Latitude
      rewind( ind )
      read(ind,'(a)',end=9901) line
      do while(index(line,'REFLAT').eq.0  )
        read(ind,'(a)',end=9901) line
      enddo
      leq = index(line,'=')
      if( line(leq+3:leq+3).eq.':' ) then
        read(line(leq+1:),'(i2,1x,i2,1x,i2)') lat, latmin, latsec
      else
        read(line(leq+1:),'(i3,1x,i2,1x,i2)') lat, latmin, latsec
      endif
      if( lat.lt.0 ) then
        xlat = float(lat) - float(latmin)/60. - float(latsec)/3600.
      else
        xlat = float(lat) + float(latmin)/60. + float(latsec)/3600.
      endif
      write(outd,'(a12,f12.4)') '>LATITUDE  =', xlat

c...Longitude
      rewind( ind )
      read(ind,'(a)',end=9902) line
      do while(index(line,'REFLONG').eq.0  )
        read(ind,'(a)',end=9902) line
      enddo
      leq = index(line,'=')
      llon = index(line(leq+1:),':') - 1
      if( llon.eq.4 ) then
        read(line(leq+1:),'(i4,1x,i2,1x,i2)') lon, lonmin, lonsec
      elseif( llon.eq.3 ) then
        read(line(leq+1:),'(i3,1x,i2,1x,i2)') lon, lonmin, lonsec
      elseif( llon.eq.2 ) then
        read(line(leq+1:),'(i2,1x,i2,1x,i2)') lon, lonmin, lonsec
      else
        write(*,*)'Unusual length for lon deg =', llon
        stop  
      endif  
      if( lon.lt.0 ) then
        xlon = float(lon) - float(lonmin)/60. - float(lonsec)/3600.
      else
	xlon = float(lon) + float(lonmin)/60. + float(lonsec)/3600.
      endif
      write(outd,'(a12,f12.4)') '>LONGITUDE =', xlon


c...Elevation
      rewind( ind )
      read(ind,'(a)',end=9903) line
      do while( index(line,'REFELEV').eq.0 )
        read(ind,'(a)',end=9903) line
      enddo
      leq = index(line,'=')
      read(line(leq+1:),*) lelev
      elev = float(lelev)
      write(outd,'(a12,f12.4)') '>ELEVATION =', elev


c...Azimuth
c...look for ZROT in impedance file
      call ediazimuth( ind, azim )
      
12    write(outd,'(a12,f12.4)') '>AZIMUTH   =', azim
      
      return

9901  write(*,*)'edicopyinfo: Could not find REFLAT in EDI file'
      return
9902  write(*,*)'edicopyinfo: Could not find REFLONG in EDI file'
      return
9903  write(*,*)'edicopyinfo: Could not find REFELEV in EDI file'
      return
9904  write(*,*)'edicopyinfo: Could not find either ZROT or ROTSPEC'//
     &          ' in EDI file'
      return

      end
