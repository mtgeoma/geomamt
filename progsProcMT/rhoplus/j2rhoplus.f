        program j2rhoplus


c  reads a Jones (J-format) or Egbert mt output file and outputs
c     a rho+ data file for both modes
c  Uses standard in. Out files have name xy_<root> and yx_<root>
c   the root is detemined from the input file name

c  Writes out a rho+ file with data and "tasks". Data have exclude
c    flags on if the coherence is less than value prompted for

c  Input insensitive to blank lines

      implicit none

       integer nfmx
       parameter(nfmx=70)
       real PI
       parameter ( PI = 3.141592 )

       real rho(nfmx,2),erho(nfmx,2),
     &      pha(nfmx,2),epha(nfmx,2),
     &      coh(nfmx,4)
       real freq(nfmx)
       real cohlim
       real errmin,errmin_rho,errmin_pha

       integer i,j,k
       integer nfreq, length, iounit
c      integer nstrt(4),nstp(4),incrho(4),incpha(4)

       character*80 infil
       character*80 outfxy,outfyx,bands
c       character*80 bands
       character*6 mssg
       character*20 root
       character*10 answer

       logical egbert, offset

       integer findstr
       external findstr

c      call setieee()

       write(0,*)'Warning: must be run in same directory as data files'
       write(0,*)'Enter y or Y to continue'
       read(5,'(a)')answer
       if(findstr(answer,'y').eq.0)stop

c  what kind of data file?

       egbert=.true.
       write(0,*)'Is this an Egbert (e,E) or Jones (j,J) file?'
       read(5,'(a)')infil
       if(findstr(infil,'e') .ne. 1)then
          if(findstr(infil,'j') .ne. 1)then
             stop 'Answer must start with E or J (case insensitive)'
          endif
          egbert=.false.
       endif

c  get input and root names
       if(egbert)then
           write(0,*)'Enter name of Egbert mt file'
       else
           write(0,*)'Enter name of Jones mt file'
       endif
       read(5,'(a)')infil

c  Following option is not documented and so has been disabled
c  Do you want the rho+ file to enable different gains in of bands?
c
c       read(5,'(a)')answer
c       if(findstr(answer,'y').ne.0)then
c          write(0,*)'Enter band index limits on one line'
c          bands(1:6)='bands '
c          read(5,'(a)')bands(7:80)
c          offset=.true.
c       else
          offset=.false.
c       endif

c  translate any nulls to blanks
       if(index(infil,char(0)).ne.0)then
         do 10 i=1,len(infil)
10         if(infil(i:i) .eq. char(0))infil(i:i)=' '
       endif

c get the data
       if(egbert)then
         call rdmteg(infil,nfmx,nfreq,freq,rho,erho,pha,epha,coh,root)
       else
         call rdmtj(infil,nfmx,nfreq,freq,rho,erho,pha,epha,coh,root)
       endif

       write(0,*)'Give error floor (in % of rho_a). 0. for no floor'
       read(5,*) errmin
       if( errmin.gt.0. ) then
         errmin_rho = errmin/100.
         errmin_pha = (180./PI)*atan(errmin/200.)
	WRITE(*,*)'errmin_rho, errmin_pha = ', errmin_rho, errmin_pha
         do k = 1, 2
           do i = 1, nfreq
             erho(i,k) = max( erho(i,k), rho(i,k)*errmin_rho )
             epha(i,k) = max( epha(i,k), errmin_pha )
           enddo
         enddo
       endif

c limit the bandwidth of the data
       call limband(nfmx,nfreq,freq,rho,erho,pha,epha,coh)

c open output files
       length=index(root,' ')-1
       outfxy(4:length+3)=root
       outfxy(1:3)='xy_'
       outfyx(4:length+3)=root
       outfyx(1:3)='yx_'
       length = length+3
       open(9,file=outfxy,err=996)
       open(10,file=outfyx,err=997)

c write a rho+ control file with data and tasks at end
c (Note, a stop line is inserted between data and tasks that must be
c  deleted to do the tasks)

       write(0,*)'You can use coherence to set the data exclude flags'
       write(0,*)'Enter y or Y if you want to set these flags'
       read(5,*)answer
       if(findstr(answer,'y') .ne. 0)then
         write(0,*)'Enter lower bound on coherence. (=0 if >.999 or <0)'
         read(5,*)cohlim
         if(cohlim .lt. 0. .or. cohlim .gt. .999)then
           write(0,'(f7.3,a)')cohlim,' not acceptable. Setting to 0.'
           write(0,*)'This means ALL data flags will be set to INCLUDE'
           cohlim=0.
         endif
       else
         cohlim=0.
       endif

c      write(0,*)'You can set the exclude data flags in 1-4 bands'
c      write(0,*)'Enter y or Y if you want to set these flags'
c      read(5,*)answer
c      if(findstr(answer,'y') .ne. 0)then
c        write(0,*)'Enter data of form start# stop# IR IP'
c        write(0,*)'Where # is frequency index (1=highest)'
c        write(0,*)'IR = 0 to exc rho; IP = 0 to exc phase'
c        i=1
c175     read(5,*,err=178)nstrt(i),nstp(i),incrho(i),incpha(i)
c        if(nstrt(i).eq.0) go to 180
c        i=1+1
c        if(i.gt.4)then
c           write(0,*)'Only 4 bands allowed'
c           go to 180
c        endif
c        write(0,*)'Enter info for another band or 0 to go on'
c        go to 175
c178     write(0,*)'Bad data. Try again or 0 to go on'
c        go to 175
c      endif
c180   continue
       
       do 200 k=1,2
         iounit=8+k
         write(iounit,'(a)')'data *'
         write(iounit,'(a)')'task *'
         write(iounit,'(a)')'surface insulating'
         if(k .eq. 1)then
           write(iounit,'(a5,a)')'root ',outfxy(1:length)
         else
           write(iounit,'(a5,a)')'root ',outfyx(1:length)
         endif
         if(offset)then
           length=index(bands,'       ')
           write(iounit,*)bands(1:length)
         endif
         write(iounit,'(a)')'model '
         write(iounit,'(a)')'plot '
         write(iounit,'(a)')'matrix '
         write(iounit,'(a)')'execute '

         do 100 j=1,nfreq
           if(coh(j,k) .gt. cohlim)then
               mssg='  1  1' 
           else
               mssg='  0  0' 
           endif
           if(pha(j,k).gt.90.0 .or. pha(j,k).lt.0)mssg='  1  0'
           write(iounit,'(5g12.4,a)')freq(j),rho(j,k),
     &          erho(j,k),pha(j,k),epha(j,k),mssg
100      continue
         write(iounit,*)'0   %end of data'
         write(iounit,*)'0   %remove this line to do tasks'
         do 150 j=1,nfreq
           write(iounit,'(g12.4,a)')freq(j),'  rho  phi'
150      continue
         write(iounit,*)'0   %end of tasks'

200    continue
       close(9)
       close(10)

       stop 'Normal end'

996    write(0,*)'Unable to open file: ',outfxy
       stop 'Abort'
997    write(0,*)'Unable to open file: ',outfyx
       stop 'Abort'
998    write(0,*)'Unable to open file: ',infil
       stop 'Abort'
       end

c***********************************************************************

       subroutine limband(nfmx,nfreq,freq,rho,erho,pha,epha,coh)

c inputs:
       integer nfmx
c inputs and outputs:
       integer nfreq
       real freq(nfmx)
       real rho(nfmx,2),erho(nfmx,2),
     &      pha(nfmx,2),epha(nfmx,2),
     &      coh(nfmx,2)

c prompts for information on how to limit bandwidth and then alters
c arrays and nfreq to reflect the choosen bandwidth

c assumes that data are monotonically ordered high frequency to low

c calls subroutines setbw, setblim, bndlm
c uses integer function findstr

       integer maxf, minf
       
       character*10 answer

       integer findstr
       external findstr
     
       write(0,*)'Enter y or Y if you want to set bandwidth'
       read(5,'(a)')answer

       if(findstr(answer,'y').ne.0)then
          call setbw(nfreq,freq,maxf,minf)
       else
          write(0,*)'Enter y or Y if you want to set data index limits'
          read(5,'(a)')answer
          if(findstr(answer,'y').ne.0)then
             call setblim(nfreq,maxf,minf)
          else
            return
          endif
       endif

       if(maxf.gt.1)then
          call bandlim(nfmx,maxf,minf,freq,rho,erho,pha,epha,coh)
       endif

       nfreq=minf-maxf+1
       write(0,*)nfreq,
     &   ' data included with indices from: ',maxf,minf

       return
       end

c***********************************************************************

       subroutine bandlim(nfmx,maxf,minf,freq,rho,erho,pha,epha,coh)

c inputs:
       integer nfmx
       integer maxf
       integer minf
c inputs and outputs:
       real freq(nfmx)
       real rho(nfmx,2),erho(nfmx,2),
     &      pha(nfmx,2),epha(nfmx,2),
     &      coh(nfmx,2)

c   maxf - index of 1st element used in arrays
c   minf - index of 1st element used in arrays

c shifts arrays freq,rho,erho,pha,epha,coh so that first entries 
c   correspond to the index maxf in original array

       integer k
   
       call Vshft(maxf-1,freq,minf)

       do 100 k=1,2
         call Vshft(maxf-1,rho(1,k),minf)
         call Vshft(maxf-1,erho(1,k),minf)
         call Vshft(maxf-1,pha(1,k),minf)
         call Vshft(maxf-1,epha(1,k),minf)
         call Vshft(maxf-1,coh(1,k),minf)
100    continue

       return
       end

c***********************************************************************

       subroutine Vshft(nshft,x,n)
       integer nshft,n
       real x(n)

c shifts values nshft places in vector x
c 
       integer i

       do 100 i=1,n-nshft
         x(i)=x(i+nshft)
100    continue
       return
       end
       
c***********************************************************************
 
      integer function nindex(zo,z,n)
c        returns the index of the element of z closest to zo.
c        assumes z(j) monotonic.
c
      integer i,n
      real zo,z(n)
 
      if( z(n).gt.z(1) )then
          do 100 i= 1, n-1
              if( .5*(z(i)+z(i+1)) .gt. zo )then
                  nindex = i
                  return
              end if
 100      continue
          nindex = n
      else
          do 200 i= n, 2, -1
              if( .5*(z(i)+z(i-1)) .ge. zo )then
                  nindex = i
                  return
              end if
 200      continue
          nindex = 1
      end if
 
      return
      end

c***********************************************************************

       subroutine setbw(nfreq,freq,maxf,minf)
c inputs: 
       integer nfreq
       real freq(nfreq)
c outputs:
       integer maxf,minf

       real fmax, fmin

       integer  nindex
       external nindex

       write(0,*)'Enter shortest period or NEGATIVE highest frequency'
       read(5,*)fmax
       if(fmax. eq. 0. .or. fmax .lt. -freq(1))then
         maxf=1
       else
         if(fmax.gt.0.) fmax= 1./fmax
         if(fmax.lt.0.) fmax= -fmax
         maxf=nindex(fmax,freq,nfreq)
       endif
 
       write(0,*)'Enter longest period or NEGATIVE lowest frequency'
       read(5,*)fmin
       if(fmin. eq. -0. .or. fmin .lt. 1./freq(nfreq))then
         minf=nfreq
       else
         if(fmin.gt.0.) fmin= 1./fmin
         if(fmin.lt.0.) fmin= -fmin
         minf=nindex(fmin,freq,nfreq)
       endif
 
       if(minf.gt.nfreq)minf=nfreq

       return
       end

c***********************************************************************
 
       subroutine setblim(nfreq,maxf,minf)
c inputs: 
       integer nfreq
c outputs:
       integer maxf,minf
 
         write(0,*)'Enter index of 1st (highest) included frequency'
         read(5,*)maxf
 
         write(0,*)'Enter index of last (lowest) included frequency'
         read(5,*)minf
         if(minf.le.maxf)then
           minf=maxf+1
           write(0,*)'Too small, setting to: ',minf
         endif
         if(minf.gt.nfreq)then
           minf=nfreq
           write(0,*)'Too large, setting to: ',nfreq
         endif

         return
         end

c***********************************************************************
 
       subroutine rdmtj(infil,nfmx,nfreq,freq,rho,erho,
     &                               pha,epha,coh,root)
c inputs:
       character*80 infil
       integer nfmx
c outputs:
       integer nfreq
       real freq(nfmx)
       real rho(nfmx,2),erho(nfmx,2),
     &      pha(nfmx,2),epha(nfmx,2),
     &      coh(nfmx,4)
       character*(*) root

c  Uses integer functions "findstr" and "lnblnk"
c  Uses character function "getline"

c  reads mt data from a Jones J file=infil

c  Input insensitive to blank lines

c  Determines root by stripping off trailing ".dat" if it exists

       real rmax,rmin,pmax,pmin

       integer i,k
       integer findstr,lnblnk

       character*120 line,getline
       character*3 modehd(2)

       external findstr, getline, lnblnk

       k=findstr(infil,'.dat')
       if (k.gt.0) then
          root=infil(1:k-1)
       else
          root=infil
       endif

       open(11,file=infil,status='old',err=998)

       modehd(1)='RXY'
       modehd(2)='RYX'

       do 500 k=1,2

c  find start of data for each mode

         rewind(11)
 50      line=getline(11)
         if(findstr(line,modehd(k)).eq.0)go to 50

         line=getline(11)
         read(line,*)nfreq

         do 100 i=1,nfreq
            line=getline(11)
            read(line,*) freq(i),rho(i,k),pha(i,k),
     &           rmax,rmin,pmax,pmin,coh(i,2*k-1),coh(i,2*k)
c convert periods (>0) to frequencies
            if(freq(i).gt.0.)freq(i) = 1./freq(i)
            if(freq(i).lt.0.)freq(i) = -freq(i)
c move YX phase to 1st quadrant
            if(k.eq.2)pha(i,k)=180.+pha(i,k)
c generate errors from stated response limits
            erho(i,k)=abs(rmax-rmin)/2.
            epha(i,k)=abs(pmax-pmin)/2.
 100     continue
 500   continue

       return

998    write(0,*)'Unable to open file: ',infil
       stop 'Error in rdmtj'

       end

c***********************************************************************
 
       subroutine rdmteg(infil,nfmx,nfreq,freq,
     &                    rho,erho,pha,epha,coh,root)
c inputs:
       character*80 infil
       integer nfmx
c outputs:
       integer nfreq
       real freq(nfmx)
       real rho(nfmx,2),erho(nfmx,2),
     &      pha(nfmx,2),epha(nfmx,2),
     &      coh(nfmx,4)
       character*(*) root

c root is derived from infil and can be used to open output files
c   with same root as the input file

c  Uses integer function "findstr"
c  Uses character function "getline"

c  reads mt data from an Egbert file=infil

c  Input insensitive to blank lines

c  Determines root by stripping off initial "mt" if it exists

       integer i,k
       integer findstr

       character*120 line,getline

       external findstr, getline

       if(infil(1:2) .eq. 'mt')then
          root=infil(3:22)
       else
          root=infil(1:20)
       endif
       open(11,file=infil,status='old',err=998)

c  read data

c find line just before data
5      line=getline(11)
       if(findstr(line,'(sec)').eq.0)go to 5

c read a line and determine whether it contains data

       i=1
10     line=getline(11)
       if(findstr(line,'impedance').eq.0)then
          read(line,*)freq(i),
     &     rho(i,1),erho(i,1),pha(i,1),epha(i,1),coh(i,1),
     &     rho(i,2),erho(i,2),pha(i,2),epha(i,2),coh(i,2)
          i=i+1
          go to 10
       else
          nfreq=i-1
       endif

c convert period to freq
         do 100 i=1,nfreq
           freq(i)=1./freq(i)
100    continue

c convert errors to 1 sd
       do 101 k=1,2
         do 101 i=1,nfreq
           erho(i,k)=erho(i,k)/2.
           epha(i,k)=epha(i,k)/2.
101    continue

       return

998    write(0,*)'Unable to open file: ',infil
       stop 'Abort'

       end

c***********************************************************************
 
       subroutine rdrsp(infil,nfmx,nfreq,freq,
     &                    rho,erho,pha,epha,coh,root)
c inputs:
       character*80 infil
       integer nfmx
c outputs:
       integer nfreq
       real freq(nfmx)
       real rho(nfmx,2),erho(nfmx,2),
     &      pha(nfmx,2),epha(nfmx,2),
     &      coh(nfmx,4)
       character*(*) root
c
c  Reads MT and D+ responses from .rsp files

c  infil must have either xy_ or yx_ prefix.
c  Determines root by stripping off path and xy_ or yx_ prefix
c    and .rsp suffix from infil

c  errors are determined as the greater of the claimed error
c    or the misfit to the D+ model response. Included data
c    have their coherences set to 1.

c  Excluded data (flagged with zero errors) are replaced with
c    D+ responses and errors equal to the misfit to D+
c    and coherences of 0.

c  Uses integer function "findstr"
c  Uses character function "getline"

c  reads mt data from a rho+ .rsp files

       parameter(nmx=100)
       real d_rho(nmx,2), d_pha(nmx,2)

       integer i,k
       integer nfyx
       integer findstr

       logical xy, yx

       character*80 infil1, infil2

       external findstr

       if(nfmx.gt.nmx)stop 'nfmx > nmx in rdrsp. Recompile'

       k=findstr(infil,'.rsp')
       if(k.eq.0)then
         stop 'Input file does not have .rsp suffix'
       else
          j1=findstr(infil,'xy_')
          j2=findstr(infil,'yx_')
          if(j1.gt.0)then
             infil1=infil
             infil2=infil
             infil2(j1:j1+2)='yx_'
          elseif(j2.gt.0)then
             infil2=infil
             infil1=infil
             infil1(j2:j2+2)='xy_'
          else
             stop 'Input file does not begin with xy_ or yx_'
          endif
       endif
       root=infil(j+3:k-1)

c  open files and read data

       xy=.true.
       yx=.true.
       open(11,file=infil1,status='old',err=5)
       go to 6
5      xy=.false.
6      open(12,file=infil2,status='old',err=7)
       go to 8
7      yx=.false.
8      continue
       if(.not. xy .and. .not. yx)stop 'No XY or YX data'
       if(.not. xy )write(0,*)'Warning: no XY data'
       if(.not. yx )write(0,*)'Warning: no YX data'

       if(xy)then
            i=1
10          read(11,*,err=12)freq(i),rho(i,1),erho(i,1),
     &              d_rho(i,1),pha(i,1),epha(i,1),d_pha(i,1)
            i=i+1
            go to 10
12          nfreq=i-1
       endif
       if(yx)then
            i=1
14          read(12,*,err=16)freq(i),rho(i,2),erho(i,2),
     &              d_rho(i,2),pha(i,2),epha(i,2),d_pha(i,2)
            i=i+1
            go to 14
16          nfyx=i-1
       endif
       if(xy .and. yx .and. (nfreq.ne.nfyx))
     &    stop 'Number of XY and YX data not equal'


c Assign errors and set exclude flags in coherences

       if(xy)then
          n=1
       else
          n=2
       endif
       if(yx)then
          m=2
       else
          m=1
       endif

       do 101 k=n,m
         do 101 i=1,nfreq

           d=abs(rho(i,k)-d_rho(i,k))
           if(erho(i,k).gt.0.)then
              erho(i,k)=max(erho(i,k),d)
              coh(i,2*k-1)=1.0
           else
              rho(i,k)=d_rho(i,k)
              erho(i,k)=d
              coh(i,2*k-1)=0.0
           endif

           p=abs(pha(i,k)-d_pha(i,k))
           if(epha(i,k).gt.0.)then
              epha(i,k)=max(epha(i,k),p)
              coh(i,2*k)=1.0
           else
              pha(i,k)=d_pha(i,k)
              epha(i,k)=p
              coh(i,2*k)=0.0
           endif
101    continue

       return
       end

c***********************************************************************
       real function rmsflg(n,x,y,flag)
c inputs:
       integer n
       real x(n), y(n), flag(n)

c  computes the rms difference between vectors x and y
c   EXCLUDING data with flag(i)=0.

       integer i, nuse

       nuse=n
       rmsf=0.

       do 100 i=1,n
         if(flag(i) .ne. 0.)then
           rmsf=rmsf+(x(i)-y(i))**2
         else
           nuse=nuse-1
         endif 
100    continue

       rmsf=sqrt(rmsf/float(nuse))

       return
       end

c***********************************************************************

c***********************************************************************
 
      character*(*) function getline(iin)
      integer iin, lnblnk
      external lnblnk

c  reads a unit until it finds a non-blank line
 
 100  read(iin,'(a)',end=200) getline
      if( lnblnk(getline).eq.0 ) go to 100
 
      return
 200  getline='end'
 
      return
      end

c***********************************************************************
 
      integer function lnblnk(string)
      character*(*) string
c     Returns the position of the last non-blank character in string
c
c     (Actually there is a standard unix function of the same name
c     that does the exactly same thing.  This is provided for
c     compatibility on non-unix systems.)
 
      do 100 lnblnk=len(string), 1, -1
          if(     string(lnblnk:lnblnk).ne.' '
     &     .and.string(lnblnk:lnblnk).ne.'      ' ) return
  100 continue
      lnblnk=0
 
      return
      end

c***********************************************************************
 
      integer function findstr(str1,str2)
      character*(*) str1, str2
c     returns the position of str2 in str1.  Ignores case.
c     returns 0 if str2 not found in str1
 
      integer i, j, capdif
      logical same
 
      capdif= ichar('a')-ichar('A')
 
      do 20 i= 1, len(str1)-len(str2)+1
         do 10 j=1,len(str2)
 
            same= str1(i+j-1:i+j-1) .eq. str2(j:j)        .or.
 
     &       'A'.le.str2(j:j) .and. str2(j:j).le.'Z' .and.
     &       ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j))+capdif .or.
 
     &       'a'.le.str2(j:j) .and. str2(j:j).le.'z' .and.
     &       ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j)) - capdif
 
            if( .not.same ) go to 20
 10      continue
         findstr=i
         return
 20   continue
 
      findstr=0
      return
 
      end


c***********************************************************************

        character*80 function getwrd(line,n)
        character line*(*)
        integer n
 
c       returns the n'th word (delineated by blanks) in line
c       returns 'null' if there are not enough words.

        integer i,ln,nstart,nend

        ln=len(line)

c skip over the first n-1 words
        nstart=1
        do 50 i=1,n-1
           nend=nstart+index(line(nstart:ln),' ')-2
           do 10 j=nend+1,ln
              if(line(j:j).ne.' ')go to 13
10         continue
13         continue
c find start of n'th word
           nstart=j
50      continue

c find end of nth word
        nend=nstart+index(line(nstart:ln),' ')-1

        if(nstart.gt.nend)then
          getwrd='null'
          return
        else
          getwrd=line(nstart:nend)
          return
        endif
        
        end
