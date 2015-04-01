      subroutine prteig(ev,ev2,u,period,nt,nf,sta,iounit,outname,
     &  id,iband,nev,nsta,ih,lstcc,stcor,declav,lhd,cvecname)
        
ccc   OLD text file output for multi- station eigenvector analysis
ccc   changed 1 Oct, 1996

      integer iband(2),ih(nsta)
      complex u(nt,nev)
      real ev(nt),ev2(nev),stcor(2,*),declav,tpower
      character*10 outname,cvecname
      character*3 sta(*)
      character*80 cformat,chd
      logical lstcc,lhd
ccc   output eigenvalues,vectors

      cformat = '(1x,a3,1x,2f7.2,5(2f8.3,1x)/19x,5(2f8.3,1x))'
      chd = '(1x,4hsta.,3x,4hlat.,6x,4hlon.)'

      if(lhd) then
ccc      output title, eigenvalues
         write(iounit,601) outname,period
601      format('         array  ',a10,
     &         '   period =   ',e12.4,' seconds  ')
        
         write(iounit,*) '# of data vectors used = ', nf
         write(iounit,*) 'Decimation level  = ',id,'; Freqs. = ',iband

         if(nf.eq.0) return

         write(iounit,*) 'Eigenvalues (normalized to local noise):'
         write(iounit,615) (ev(k),k=1,nt)
         write(iounit,*) 'Fraction of Total Power (SNR units)'
         tev = 0.0
         do k = 1,nt
            tev = tev+ev(k)
         enddo 
         write(iounit,614) (ev(k)/tev,k=1,nt)
         write(iounit,*) 'Power in PCs of Sig. Comp. (nT**2/Hz):'
         write(iounit,615) (ev2(k),k=1,nev)
         write(iounit,*) 'Fraction of Total Power (nT**2/Hz)'
         tev = 0.0
         do k = 1,nev
            tev = tev + ev2(k)
         enddo
         write(iounit,614) (ev2(k)/tev,k=1,nev)
615      format(8f10.2)
614      format(8f10.6)
      end if

ccc   write out vectors in U'
      write(iounit,*)
      do j=1,nev
         write(iounit,'(a10,i3)') cvecname,j
         write(iounit,chd) 
         do i=1,nsta
            write(iounit,cformat) sta(i),stcor(1,i),stcor(2,i),
     &                  (u(k,j),k=ih(i),ih(i+1)-1)
         enddo
      enddo
      write(iounit,*) '_______________________________________',
     & '_______________________________________'
      return
      end
