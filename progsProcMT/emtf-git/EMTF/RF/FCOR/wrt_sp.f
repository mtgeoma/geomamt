
      subroutine wrt_sp(fid,stname,stcor,decl,nch,samprate,clock,
     &  dipole_length,gain,sensitivity,nf,nfmax,ftype)

ccc   writes out sp******* file in standard format ... 
ccc   dipole_length is in  km
ccc   
      integer fid,nch,nf(nch)
      real stcor(2),decl,samprate,orient(nch),gain(nch),clock(2),
     &     sensitivity(nch)
      character*(*) stname,ftype(nfmax,nch)
      character*1 chid(nch)
      write(fid,'(a20)') stname
      write(fid,'(2f10.4)') stcor
      write(fid,'(f10.4)') decl
      write(fid,'(i2)') nch
      write(fid,'(e12.4)') samprate
      write(fid,'(2f10.4)') clock
      do ich = 1,nch
         write(2,'(a1)') chid(ich)
         if(chid(ich).eq.'E') then
            write(fid,'(4f10.4)') dipole_length,orient(ich),0.,gain(i)
         else
            write(fid,'(2f10.4)') orient(ich),0.
         endif
         write(fid,'(f10.4,i3)') sensitivity(ich),nf(ich)
         do k = 1,nf(ich)
            write(fid,*) ftype(k,ich)
         enddo 
      enddo
      return
      end
