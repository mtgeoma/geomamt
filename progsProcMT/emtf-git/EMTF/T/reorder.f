      subroutine reorder(x,ncht,nfreq,iin,nchout,iout,lrref,irref)
ccc   sets everything up to allow for any two channels
ccc   to be used as the input reference, any other two
ccc   channels to be used for a remote reference, and any
ccc   subset of the remaining channels to have TF estimated

ccc   2  channels to be used as the input reference are given in iin
ccc   if lrref = .ture., 2 channels to be used for remote are
ccc     given in irref 
ccc   channels to estimate TFs for are given in iout

      integer NCHMX_LOC
      parameter (NCHMX_LOC = 100) 
      complex xtemp(NCHMX_LOC),x(ncht,nfreq)
      integer iout(nchout),iin(2),irref(2)
      logical lrref

      if(lrref) then
         nchstk = nchout+4
      else
         nchstk = nchout + 2
      endif


      do i = 1,nfreq
         xtemp(1) = x(iin(1),i)
         xtemp(2) = x(iin(2),i)
         do k= 1,nchout
            xtemp(k+2) = x(iout(k),i)
         enddo
         if(lrref) then
            xtemp(3+nchout) = x(irref(1),i)
            xtemp(4+nchout) = x(irref(2),i)
         endif
         do k = 1,nchstk
            x(k,i) = xtemp(k)
         enddo 
      enddo
      return
      end
