      subroutine date(string)
      character*(*) string
      integer iarray(3)
	external idate

      call idate(iarray)
      write(string(1:2),'(i2.2)') iarray(1)
      string(3:3) = '/'
      write(string(4:5),'(i2.2)') iarray(2)
      string(6:6) = '/'
      write(string(7:8),'(i2.2)') iarray(3) - (iarray(3)/100*100)

      return
      end


