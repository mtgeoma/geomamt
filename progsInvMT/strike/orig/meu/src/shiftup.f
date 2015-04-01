      subroutine shiftup (string)

      implicit none

      character*(*) string
      integer table(256), ns, i, ic, tc
c converts lowercase letters in string to uppercase

      data table /
     1  1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
     2 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
     3 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
     4 49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,
     5 65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
     6 81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,
     7 65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
     8 81,82,83,84,85,86,87,88,89,90, 123,124,125,126,127,128,
     9 128*0/

      include 'ctrlblk.inc'

      if( iprint.ge.3 ) then
        write(*,*)'SHIFTUP: entered with string=',string
      end if

      if( string.eq.' ' .or. string(1:1).eq.char(0) ) return

      call trunc( string, ns )
      if( iprint.ge.3 ) then
        write(*,*)'SHIFTUP: after TRUNC, ns, string=',
     &            ns, ' ', string(:ns)
      end if

      if( ns.eq.0 ) return

      do 1 i = 1, ns
        ic = ichar(string(i:i))
        tc = table(ic)
        string(i:i) = char(tc)
        if( iprint.ge.10 ) then
          write(*,*)'i, ic, string(i:i)=',i,ic,' ',string(i:i)
        end if
    1 continue

      if( iprint.ge.5 ) then
        write(*,*)'SHIFTUP: leaving with string=',string(:ns)
      end if

      return
      end  
