      subroutine trunc (line, length)
      character*(*) line
      integer length,start,mid,end

      include 'ctrlblk.inc'

c
c this routine trims leading and trailing blanks from 'line', returning
c its new size in 'length'.  the actual length of the string remains
c unchanged; it is padded on the right with blanks if leading blanks are
c removed, but 'length' is not increased.
c
 
      length = len( line )
c
c      lastbl = len_trm( line )
c              ---- lnblnk is a sun fortran library function
c                   that returns the last non-blank char in a string
      lastbl = lnblnk( line )

      if( iprint.ge.20 ) then
        write(*,*)'TRUNC: entered with line=',line
        write(*,*)' input length & last non-blank=',length, lastbl
      end if

      if( length.eq.0 .or. lastbl.eq.0 ) then
        length = 0
        return
      end if
 
      do 10 start = 1, length
c	write(*,*)'  start,line(start:start)=',start,' ',
c    &            line(start:start)
        if( line(start:start).ne.' ' ) goto 1
   10 continue
      length = 0
      return
 
    1 if( start.gt.1 ) then
        line(1:length-start+1) = line(start:length)
        line(length-start+2:length) = ' '
      end if
 
c
c do a binary search for the last non-blank character in 'line'.
c
      START = 1
      END = LEN(LINE)
      DO WHILE (START .LT. END)
        MID = (START+END)/2
        IF (LINE(MID:END) .EQ. ' ') THEN
           END = MID-1
        ELSE
           START = MID+1
        END IF
      END DO
      IF (LINE(START:START) .EQ. ' ') THEN
        LENGTH = START-1
      ELSE
        LENGTH = START
      END IF

c	write(*,*)'TRUNC: leaving with length=',length
      RETURN
      END
