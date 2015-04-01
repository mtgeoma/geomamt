      subroutine edicopycomm( ind, outd )

      implicit none

      integer ind, outd

      integer lnblnk
      character*80 line

c...copies comment block from EDI input to output file 
c   (both must be open already)

      write(outd,'(a)') '#'
      write(outd,'(a)') '# Comment (INFO) block from edi file'
      write(outd,'(a)') '#'

      rewind( ind )

      read(ind,'(a)',end=9999) line

      do while(line(1:5).ne.'>INFO'  )
        read(ind,'(a)',end=9999) line
      enddo

      do while( line(1:1).ne.'>' )
        write(outd,'(2a)') '# ', line(1:lnblnk(line))
        read(ind,'(a)',end=9999) line
      enddo


      backspace(ind)

9999  return
      end
