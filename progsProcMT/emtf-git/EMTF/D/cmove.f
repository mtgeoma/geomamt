
c
c**********************************
c
      subroutine cmove(c1,c2,nmv,i1)
      character*1 c1(*)
      character*120 c2

      do 5 i = 1,nmv
      j = i1 - 1 + i
      c1(j) = c2(i:i)
5     continue

      i1 = i1 + nmv

      return
      end
