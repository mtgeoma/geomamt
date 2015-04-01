      subroutine jcopyinfo( ind, outd )
      integer ind, outd
      character*80 line

c...copies info block from input to output file (both must be open already)

      read(ind,'(a)',end=9999) line

      do while(line(1:1).eq.'>'  )
        if(outd.ne.0) write(outd,'(a)') line(1:lnblnk(line))
        read(ind,'(a)',end=9999) line
      enddo

      backspace(ind)

9999  return
      end
