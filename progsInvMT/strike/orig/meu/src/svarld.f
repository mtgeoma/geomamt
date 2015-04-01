      subroutine svarld(data,weight,mean,n,svar)

      real data(n),weight(n),mean,svar,xn,sumw
c 
c     subroutine to calculate the weighted sample variance
c        For linear data
c
c     note: mean derived if set to -999. on input
c           weights unused if data(1) = weight(1)
c           i.e., called  svarld(data,data,...)
c
c 
c MISDAT     - missing data marker; set to -999. in MISDAT.inc
      include 'misdat.inc'
      include 'ctrlblk.inc'

        if( iprint.ge.1 ) then
          write(*,*)'svarld: entered with n, mean =', n, mean
        endif


      if(n.le.1) then
        write(*,*)' svarld: impermissable n = ', n
        svar=0.0
        return
      end if

      xn=Float(n)

      if( mean.eq.-999. ) then
        mean = 0.0
        wtot = 0.0
        do 1 i = 1, n
          if( data(i).eq.MISDAT ) goto 1
          if( data(1).eq.weight(1) ) then
            mean = mean + data(i)
            wtot = wtot + 1.
          else
            mean = mean + data(i)*weight(i)
            wtot = wtot + weight(i)
          end if
    1   continue
        mean = mean / wtot
      end if
      
c 
      svar=0.0
      if(data(1).ne.weight(1))goto 3

      xn = 0.
      do i=1,n
        if( data(i).ne.MISDAT ) then
          xn = xn + 1.
          svar=svar+(data(i)-mean)**2
        endif
      enddo
      svar=svar/(xn-1.0)
      goto 9999
c 
    3 sumw=0.0
      do i=1,n
        if( data(i).ne.MISDAT ) then
          sumw=sumw+weight(i)
          svar=svar+weight(i)*((data(i)-mean)**2)
        endif
      enddo
      svar=(xn/(xn-1.))*svar/sumw

9999  if( iprint.ge.2 ) then
        write(*,*)'svarld: returning with mean, svar =',mean,svar
      endif

      return
      end

