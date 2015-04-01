      real function ddmmss2deg(line)
      character*(*) line
      integer sign,deg,min,i,j
      real sec

      sign = 1
      i = 1
c     skip leading blanks
      do while( line(i:1).eq.' ' )
         i = i + 1
      enddo

      if(line(i:1).eq.'-') then
         sign = -1
         i = i + 1
      endif

      j = index(line,':')-1
      read(line(i:j),*) deg
      i = j+2
      j = index(line,':',back=.TRUE.)-1
      read(line(i:j),*) min
      i = j+2
      j = lnblnk(line)
      read(line(i:j),*) sec

      ddmmss2deg = sign*(float(deg)+float(min)/60.+sec/3600.)
      return
      end
