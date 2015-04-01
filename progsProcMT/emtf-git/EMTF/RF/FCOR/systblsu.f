c___________________________________________________
c
        subroutine systblsu(dr,nch,nfil,afparam,iftype,sc,cda,cfsTF)
        
c       computes table of transfer functions for subroutine fctf
c       for uncorrecting transfer functions back to measurement units

ccc         corrections include
c       3) correction of units of fourier coefficients (to counts/(sqrt(hz)))
c       4) analogue low pass (anti-alias) filters of instruments
c       5) analogue high pass filters of instruments (if any)
c        5') system response table
c       6) scale factors for individual channels to convert form counts
c                to physically meaningful units; for magnetic channels
c                this will be nt for electric mv/km; transformation of
c                measured channels to a "nice" coordinate system is not
c                done with this 

c       =========> afparam(nfilmax,nch) analogue filter/system response parameters
c>>>>>>>>>> 28 Feb, 1991 changed to character*80 array
c       =========> iftype(nfilmax,nch) indicator variable for formula for 
c                    filter/response
c       =========> sc(nch) = scale factors to convert counts to physical
c                       units
c       =========> cda clock offset in seconds (+ for fast)
c      
c       <========  systbl(nch,nsystbl) = conmplex table output
      
      include 'fcor.inc'
      parameter (pi=3.141592654)

      real dr,sc(nch),freq(0:nsystbl)
      complex systbl(nchmx,nsystbl),afcor,temp,tc
      integer iftype(nfilmax,nch),nfil(nch)
      character*80 afparam(nfilmax,nch)
      character*40 cfile,ctemp
      character*80 cfsTF

      if(nch.gt.nchmx) then
         write(0,*) 'ERROR IN systblsu !!!!'
         write(0,*) 'NCH = ',nch,'  NCHMX = ',nchmx
         write(0,*) 'STOPPING'
         stop
      endif
c>>>> correct units of fc's & for fixed clock offset
      freq(0) = steplog/(2.*dr)
      do j = 1,nsystbl
         freq(j) = freq(j-1)/steplog
         tc = cexp(cmplx(0.,freq(j)*2.*pi*cda))
         do i = 1,nch
            systbl(i,j) = sc(i)*sqrt(dr)
         enddo
      enddo

c>>>>>   corrections for analogue instrument filters/system response,
      do 40 ich = 1,nch
         do 35 k = 1,nfil(ich)
c         read info for each filter out of character array, set
c               up tables for table look up of system response
         j = iftype(k,ich)
         if((j.eq.3).or.(j.eq.4).or.(j.eq.6)) then
            read(afparam(k,ich),*) gain,t0,alpha
         else if((j.eq.2).or.(j.eq.5)) then
            read(afparam(k,ich),*) gain,t0
         else if(j.eq.7) then
            gain = 1.0
            read(afparam(k,ich),*) cfile,imode,igain
            ctemp = 'sensors/'//cfile
            call resptbl(ctemp,j,imode,igain,freq,temp,0)
         else if(j.eq.8) then
            gain = 1.0
            read(afparam(k,ich),*) cfile
            ctemp = 'sensors/'//cfile
            call resptbl(ctemp,j,imode,igain,freq,temp,0)
         end if

         do i = 1,nsystbl
            period = 1./freq(i)
c        corect for filters, system response
            if(j.le.6) then
               temp = afcor(j,t0,period,alpha)
            else if((j.eq.7).or.(j.eq.8)) then
               call resptbl(cfile,j,imode,igain,freq(i),temp,1)
            else if(j.eq.9) then
               temp = 1.0
               gain = 1.0
            end if
            temp = temp*gain
            systbl(ich,i) = systbl(ich,i)/temp
         enddo ! do i = 1,nsystbl

35       continue 
40    continue
                       
       open(unit=93,file=cfsTF)
       write(93,'(2i5)')  nch,nsystbl
       do i = 1,nsystbl
          write(93,'(11e10.4/10x,10e10.4)') freq(i),(systbl(j,i),
     .          j=1,nch)
       enddo
      close(93)
      return
      end
