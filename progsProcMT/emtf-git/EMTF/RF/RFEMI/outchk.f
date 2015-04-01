      subroutine outchk(ix,nch,nout,iout,ir,iseg,nseg)

      integer iout,nout,nch,ir,iseg(2,*),nseg,ix(0:nch,nout)

      nchp1 = nch+1
      it1 = ix(0,1)
      if(it1 .gt. ir) then
ccc      start new segment
         iseg(2,nseg) = ir-1
         nseg = nseg+1
         iseg(1,nseg) = it1
         ir = it1 + nout
      else if(it1 .lt. ir) then
         write(0,*) 'ERROR !!! it1 > ir '
         write(0,*) 'it1,ir = ',it1,ir
         stop
      else
         ir = ir + nout
      endif

      if(ix(0,nout) .ge. ir ) then
ccc      must be gaps in block
ccc      not programed for this case yet
         write(0,*) 'ERROR!!!  Gaps in block'
         stop
      else if(ix(0,nout) .lt. ir-1) then
         write(0,*) 'ERROR!!! record numbers wrong'
         stop
      else
         call outdo(it1,nout,nchp1,ix,iout)
      endif
      return
      end
