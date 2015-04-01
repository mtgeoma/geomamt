ccc_____________________________________________________________________
ccc
      subroutine asc_rec(nch)

ccc   asc_rec  reads in clock reset information (time of zero record +
ccc    time of universal clock zero), and sets up sample numbering for
ccc    ASCII data files

      include 'iounits.inc'
      include 'input.inc'

      integer mday(12)
      data mday/31,28,31,30,31,30,31,31,30,31,30,31/

      integer irec
      character*80 cform
      character*256 cline

      common/ASCBLK/irec,cform

ccc   set format for data files
ccc   For more specific formats the number of channels might have
ccc   to be hard coded; set here, and remove code which counts channels on
ccc   one line at end of file (this will still work if channels are separated
ccc   by blanks)
      cform = '*'
 
      if(l_clk_hd) then
ccc      read clock reset info from first three lines of data file
ccc      SAMPLING RATE (seconds)
         read(in_unit,*) dr
ccc      INSTRUMENT CLOCK ZERO TIME - year month, day, time (ut)
ccc          [ hour,min,sec] (i.e., the time of the first record in the file)
         read(in_unit,*) iyr,imo,iday,ihr,imin,isec
ccc      UNIVERSAL CLOCK ZERO TIME - year month, day, time (ut) 
ccc          [ hour,min,sec] ( i.e., reference time for synchronizing stations)
         read(in_unit,*) iyru,imou,idayu,ihru,iminu,isecu
      else 
ccc      read header info from clock reset file
         open(unit = clk_unit, file = cfclk)
ccc      SAMPLING RATE (seconds)
         read(clk_unit,*) dr
ccc      INSTRUMENT CLOCK ZERO TIME - year month, day, time (ut)
ccc          [ hour,min,sec] (i.e., the time of the first record in the file)
         read(clk_unit,*) iyr,imo,iday,ihr,imin,isec
ccc      UNIVERSAL CLOCK ZERO TIME - year month, day, time (ut) 
ccc          [ hour,min,sec] ( i.e., reference time for synchronizing stations)
         read(clk_unit,*) iyru,imou,idayu,ihru,iminu,isecu
         close(clk_unit)
      endif
 
ccc   use clock reset info to set sample numbering for data in ASCII file
ccc    compute julian day number for instrument clock zero
      jday = 0
      do i = 1,imo-1
         if( (mod(iyr,4).eq.0) .and. (i.eq.2) ) then
            jday = jday+mday(i) + 1
         else
            jday = jday + mday(i)
         end if
      enddo
      jday = jday + iday
 
ccc   compute julian day number for universal clock zero
      jdayu = 0
      do i = 1,imou-1
         if( (mod(iyru,4).eq.0) .and. (i.eq.2) ) then
            jdayu = jdayu+mday(i) + 1
         else
            jdayu = jdayu + mday(i)
         end if
      enddo
      jdayu = jdayu + idayu
 
ccc   record number (relative to universal clock zero) of first record in file
      jday = jday-jdayu
      ihr = ihr - ihru
      imin = imin - iminu
      isec = isec - isecu

      time_int = (((jday*24+ihr)*60+imin)*60+isec)
      sampfreq = 1./dr
      irec = 0

ccc   count number of channels in one line
      print*, 'Starting sample # = ',irec
      read(in_unit,'(a256)') cline
      m1 = icfirst(cline,256)
      m2 = iclong(cline,256)
      nch = 1
      do m = m1,m2-1
         if( (cline(m:m).eq. ' ').and.(cline(m+1:m+1).ne.' '))
     &     nch = nch+1
      enddo
      write(0,*) 'NCH  = ',nch
ccc   repositionascii input file at beginiing of data
      rewind(in_unit)
      if(l_clk_hd) then
         read(in_unit,*) 
         read(in_unit,*) 
         read(in_unit,*) 
      endif
      return
      end
ccc_____________________________________________________________________
ccc
      subroutine rdasc(npts,ix,ierr,nch)

c      ix = data array (integer)
c      npts : on input = number of points to try reading returned
c           : on retrun =  number of points actually returned
c      ierr = 0 for success, -1 for end, -3 for big gap, -2
c      nch = number of data channels

      include 'iounits.inc'
      include 'input.inc'
      include '../include/datsz.inc'

      integer irec
      character*80 cform
      common/ASCBLK/irec,cform
 
      integer ix(0:nch,*),ierr,nch
         
      npts = nblk
      ierr = 0
      if(l_asc_rec_num) then
         print*, 'This option not yet available'
         print*, ' (ASCII files with sample numbers)'
         print*, ' Stopping'
         stop
      else
ccc      read in consecutive data records with no sample numbers
         if(cform(1:1).eq.'*') then
            do j = 1,npts
               read(in_unit,*,end=20) (ix(k,j),k=1,nch)
               ix(0,j) = irec
               irec = irec + 1
            enddo
         else
            do j = 1,npts
               read(in_unit,cform,end=20) (ix(k,j),k=1,nch)
               ix(0,j) = irec
               irec = irec + 1
            enddo
         endif
c         do j =1,npts
c            write(0,*) (ix(k,j),k=0,nch)
c         enddo
         return
 
20       continue
         ierr = -1
         npts = j-1
         return
      endif
      end
