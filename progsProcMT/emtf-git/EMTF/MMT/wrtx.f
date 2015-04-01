      subroutine openx(ioux,cfile)
      character*80 cfile
      write(*,*) 'ioux,cfile',ioux,cfile
      open(unit=ioux,file=cfile,form='unformatted')
      return
      end
c______________________________________________________________________
c
      subroutine wrtx(ioux,x,nt,nseg,nfb,var,comment)
      character*80 comment
      integer nt,nseg,nfb
      real var(nt)
      complex x(nfb,nseg,nt)
      write(ioux) nt,nseg,nfb
      write(ioux) comment
      write(ioux) var
      write(ioux) x
      return
      end
