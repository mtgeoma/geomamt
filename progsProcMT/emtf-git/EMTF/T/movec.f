c
c**********************************
c
      subroutine movec(c1,c2,nmv,i1)
      character*1 c1(*)
      character*120 c2

      do 5 i = 1,nmv
      j = i1 - 1 + i
      c2(i:i)=c1(j)
5     continue

      i1 = i1 + nmv

      return
      end
