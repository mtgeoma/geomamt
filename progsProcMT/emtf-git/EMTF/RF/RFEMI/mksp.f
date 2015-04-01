      subroutine mksp(cfile,chead,nhd,ichd,isfreq,nch,clock,chid,dr)

c       program reads header from EMI files and makes
c       system parameter file (sp* ) for dnff
ccc     VERY UGLY ::: way to complex, written before I really understood
ccc     EMI header formats.    (I'm not sure I do yet!!!)

      include '../../include/nchmx.inc'
      common /FORMBLK/ir0,irelrc,irec0,nt,nbl,iord
      integer iord(nchmx),iordinv(nchmx),stname_length
      character*1 chead(81,50)
      integer ichd(nhd),nhd,isfreq,nbl
      real stcor(2),dr,clock(2),orient(nchmx),
     &      sensitivity(nchmx),electrode(nchmx)
      integer nch,idb(nchmx),ieboxg(nchmx),
     &      ieboxf(nchmx),ieboxm(nchmx),iemode(nchmx)
      logical lefield(nchmx),l_x(nchmx),l_z(nchmx),l_y(nchmx)
      character*12 sysfile(nchmx)
      character *20 stname
      character*10 csearch
      character*1 c1
      character*2 c2,ftype(nchmx),chid(*)
      character*3 c3
      character*10 c10
      character*60 cline
      character*80 cfile

c     find station name
      csearch = 'SITE    1:'
      call findhd(csearch,chead,nhd,iline)
      if(iline.eq.-1) then
         write(0,*) 'ERROR:::: STATION NAME NOT FOUND',nhd
         stname = 'sta'
      else
         istart = ichd(iline)
         call strpend(chead,iline,length,istart)
         if(length.eq.0) then
            write(0,*) 'ERROR:::: STATION NAME NOT FOUND',nhd
            stname = 'sta'
          else
             stname_length= length
             do 7 i = istart,istart + length - 1
7                stname(i-istart+1:i-istart+1) = chead(i,iline)
                 print*,'stname = ',stname,stname_length
          end if
       end if

ccc    open system parameter file
ccc    mr is length of "file root" -- data file name with suffix 
ccc   (begining with a dot (.) and any blanks stripped off
       mr = irlong(cfile,40)

       open(unit=2,file=cfile(1:mr)//'.sp')

c      find number of channels
      csearch = 'NO OF CH :'
      call findhd(csearch,chead,nhd,iline)
      if(iline.eq.-1) then
         print*,'enter number of channels'
         read(5,*) nch
      else
         istart = ichd(iline)
         call strpend(chead,iline,length,istart)
         if(length.eq.1) then
            c1 = chead(istart,iline)
            read(c1,'(i1)') nch
         else
c             (asssuming  length.eq.2)
            c2(1:1) = chead(istart,iline)
            c2(2:2) = chead(istart+1,iline)
            read(c2,'(i2)') nch
         end if
      end if
      write(0,*) 'nch = ',nch

c     sampling rate (in hz)
c      csearch = 'DATA TYPE:'
c      call findhd(csearch,chead,nhd,iline)
      csearch = 'PARAMETER:'
      call findhd(csearch,chead,nhd,iline)
      do 1020 i = 1,60
         cline(i:i) = chead(10+i,iline)
1020     continue
      read(cline,1090) nblks,nsegs,fhp,flp,sfreq,nbl
      write(0,*)'nblks,nsegs,fhp,flp,sfreq,nbl',nblks,nsegs,fhp,flp,
     &   sfreq,nbl
      isfreq = nint(sfreq)

c      get system info for each channel
      do 100 k = 1,nch
         if(k.lt.10) then
            write(csearch,50) k
         else
            write(csearch,51) k
         endif
50       format('CHANEL  ',i1,':')
51       format('CHANEL ',i2,':')
         call findhd(csearch,chead,nhd,iline)
         if(iline.eq.-1) then
            print*,'forget it'
            stop
         end if
         istart = ichd(iline)
c         call strpend(chead,iline,length,istart)
         c2(1:1) = chead(istart-5,iline)
         c2(2:2) = chead(istart-4,iline)
         read(c2,'(i2)') idb(k)
c         print*,'k,idb(k)',k,idb(k)
         istart = istart-6
         call strpend(chead,iline,length,istart)
         istart = istart-1
         call strpend(chead,iline,length,istart)
         istart = istart-1
         call strpend(chead,iline,length,istart)
c           after three pieces of junk we get to the orientation
ccc         which might well be meaningless in general!!!
         istart = istart-1
         call strpend(chead,iline,length,istart)
         iend = istart + length - 1
         c3(1:1) = chead(iend-2,iline)
         c3(2:2) = chead(iend-1,iline)
         c3(3:3) = chead(iend,iline)
         if(c3(2:2).ne.' ') then
            read(c3,'(f3.0)') orient(k)
         else
            read(c3(3:3),'(f1.0)') orient(k)
         endif
          print*,'k,orient(k),',k,orient(k)
c          then to the system parameter file identifier
         istart = istart-1
         call strpend(chead,iline,length,istart)
         do i = 1,length
            sysfile(k)(i:i) = chead(istart+i-1,iline)
         enddo
c         print*,sysfile(k)
ccc      added Feb. 96 ... 
ccc      find out which channels are "x" andy which "y"
         l_x(k) = (chead(12,iline).eq.'x') 
         l_y(k) = (chead(12,iline).eq.'y') 
         l_z(k) = (chead(12,iline).eq.'z') 
         
c         get extra info for electric channels
          write(csearch,'(6hCHAEXT,i3,1h:)') k
          call findhd(csearch,chead,nhd,iline)
          if(iline.eq.-1) then
ccc          not finding a CHAEXT means this is an H 
             lefield(k) = .false.
             sensitivity(k) = 1.0
             go to 100
          else
ccc          finding a CHAEXT means this is an E
             lefield(k) = .true.
    
             istart = 11
             call strpstrt(chead,iline,length,istart)
             write(0,*)  'length,iend',length,iend
             iend = istart
             ie1 = iend-length+1
             ie2 = iend
             do kk = 1,length
                c10(kk:kk) = chead(iend-length+kk,iline)
             enddo
             write(0,*) length
             read(c10(1:length),*) electrode(k)
             print*,'electrode(k) = ',electrode(k)
             istart = istart + 1
             call strpstrt(chead,iline,length,istart)
             c1 = chead(istart-length+1,iline)
             read(c1,'(i1)') ieboxg(k)
ccc             print*,'ieboxg = ',ieboxg(k)
             istart = istart + 1
             call strpstrt(chead,iline,length,istart)
             if(chead(istart-length+1,iline).eq.'L') then
                ieboxm(k) = 0 
             else
                ieboxm(k) = 1
             end if
      
             istart = istart + 1
             call strpstrt(chead,iline,length,istart)
             if(chead(istart+1,iline).eq.'N') then
                ieboxf(k) = 1 
             else
                ieboxf(k) = 0
             end if
             iemode(k) = 1  + ieboxf(k) + 2*ieboxm(k)
             sensitivity(k) = .001
          endif
100      continue

c      print*,'ieboxg',ieboxg
c      print*,'ieboxf',ieboxf
c      print*,'ieboxm',ieboxm
c      print*,'idb',(idb(k),k=1,5)

      dr = 1./isfreq
      do k = 1,nch
         sensitivity(k) = 1./(sensitivity(k)*(10.**(idb(k)/20.)))
         ll = 8
         if(sysfile(k)(8:8) .eq. ' ') ll = 7
         do l = 1,ll
            iic = ichar(sysfile(k)(l:l))
            if((iic.ge.65).and.(iic.le.90)) then
               sysfile(k)(l:l) = char(iic+32)
            endif
         enddo
         sysfile(k) = sysfile(k)(1:ll)//'.rsp'
      enddo

ccc   stcor should give staton lat/lon; decl should be declination
      stcor(1) = 0.0
      stcor(2) = 0.0
      decl = 0.0
      write(2,'(a)') stname(1:stname_length)
      write(2,'(2f10.4)') stcor
      write(2,'(f10.4)') decl
      write(2,'(i2)') nch
      write(2,'(e12.4)') dr
      write(2,'(2f10.4)') clock
      
ccc######################################################################

ccc   find order of channels, ids (will resort so that H channels are first
      call find_ord(nch,lefield,l_x,l_y,l_z,iordinv,chid,ftype)
      do  k = 1,nch
         do l = 1,nch
            if(iordinv(l).eq.k) then
               iord(k) =l
            endif
         enddo
      enddo
      write(0,*) 'iord = ',iord
ccc   Using orientation of 1st channel, and x-y defs change other
ccc   orientations to agree
      orient_x = orient(iordinv(1))
      orient_y = orient_x + 90.
      do k = 1,nch
         if(l_x(k)) orient(k) = orient_x
         if(l_y(k)) orient(k) = orient_y
         if(l_z(k)) orient(k) =  0.0
      enddo

      do k = 1,nch
         l = iordinv(k)
         write(2,'(a2)') chid(k)
         if(lefield(l)) then
            write(2,'(4f10.4)') 1.e-3*electrode(l),orient(l),0.,1.
         else
            write(2,'(2f10.4)') orient(l),0.
         endif
         write(2,'(f10.4,i3)') sensitivity(l),3
         write(2,'(a2)') ftype(k)
         if(lefield(l)) then
            write(2,1010) sysfile(l),iemode(l),ieboxg(l)
         else
            write(2,1009) sysfile(l)
         endif
         write(2,'(a2)') 'LP'
         write(2,'(f10.5)') flp
         write(2,'(a2)') 'HP'
         write(2,'(f10.5)') fhp
      enddo
      close(2)

      do k = 1,nch
         iord(k) = iord(k) + 1
      enddo
      return

1004  format(i2,8x,2i2,4x,4i2,30x,5i2)
1005  format(2e10.4)
1007  format(5a8)
1008  format(f10.4)
1009  format(1h',a12,1h')
1010  format(1h',a12,1h',2i2)
1090  format(2i10,3f10.5,i10)
      end
c______________________________________________________________________
c
      subroutine findhd(csearch,chead,nhd,iline)
      character*10 csearch,ctest
      character*1 chead(81,50)
      do 10 i = 1,nhd
      do 5 j = 1,10
         ctest(j:j) = chead(j,i)
5        continue
         if(ctest.eq.csearch) then
            iline = i
            return
         end if
10       continue
         iline = -1
         return
      end
c______________________________________________________________________
c
      subroutine strpend(chead,iline,length,istart)
      character*1 chead(81,50)
      integer length,iline
c       strips off tail of header line (starting from position istart
c       and going backwards until first a non-blank, and
c       then a blank is found; returns
c      the start position and the length of the non blank part of the
c       string
      length = 0
      do i = istart,11,-1
         if((chead(i,iline).ne.' ').and.
     &      (chead(i,iline).ne.'(').and.
     &      (chead(i,iline).ne.')').and.
     &      (chead(i,iline).ne.',')) then
             length = length + 1 
          else
             if(length .gt. 0 ) then
                istart = i+1
                return
             end if
          end if
       enddo
       istart = 11 
       return
       end
c______________________________________________________________________
c
      subroutine strpstrt(chead,iline,length,istart)
      character*1 chead(81,50)
      integer length,iline
c       strips off tail of header line (starting from position istart
c       and going backwards until a blank is found; returns
c      the start position and the length of the non blank part of the
c       string
      length = 0
      do 5 i = istart,80,1
         if((chead(i,iline).ne.' ').and.
     &      (chead(i,iline).ne.'(').and.
     &      (chead(i,iline).ne.')').and.
     &      (chead(i,iline).ne.',')) then
             length = length + 1 
          else
             if(length .gt. 0 ) then
                istart = i-1
                return
             end if
          end if
5         continue
       istart = 11 
       return
       end
ccc_____________________________________________________________________
ccc
      subroutine find_ord(nch,lefield,l_x,l_y,l_z,iordinv,chid,ftype)

      integer iordinv(nch)
      logical lefield(nch),l_x(nch),l_z(nch),l_y(nch),l_ok
      character*2 chid(nch)
      character*2 ftype(nch)

      l_ok = .true.
      nh = 2
      do k = 1,nch
         if((.not.lefield(k)).and.(l_x(k))) then
ccc         This is Hx  (If there's more than 1 ... tough!)
            ftype(1) = 'TB'
            iordinv(1) = k
            chid(1) = 'Hx'
            go to 10
         endif
      enddo 
      l_ok = .false.
10    continue
      do k = 1,nch
         if((.not.lefield(k)).and. (l_y(k))) then
ccc         This is Hy  (If there's more than 1 ... tough!)
            ftype(2) = 'TB'
            iordinv(2) = k
            chid(2) = 'Hy'
            go to 20
         endif
      enddo 
      l_ok = .false.
ccc   (both Hx and Hy are required ... currently)
20    continue
      if(.not.l_ok) then
ccc         write(0,*) 'Hx, Hy not found; stopping'
ccc         stop
         nh = 0
      else
         iordinv(3) = 0
         do k = 1,nch
            if((.not.lefield(k)).and. (l_z(k))) then
ccc            This is Hz  (If there's more than 1 ... tough!)
               nh = nh+1
               ftype(3) = 'TB'
               iordinv(3) = k
               chid(3) = 'Hz'
               go to 30
            endif
         enddo 
      endif
30    continue
      do k = 1,nch
         if(.not.lefield(k)) then
            if((k.ne.iordinv(1)).and.(k.ne.iordinv(2)).and.
     &         (k.ne.iordinv(3))) then
ccc         This is yet another H ... just add after the first Hx,Hy,Hz
               nh =   nh+1
               ftype(nh) = 'TB'
               iordinv(nh) = k
               if(l_x(k)) then
                  chid(nh) = 'Hx'
               else if(l_y(k)) then
                  chid(nh) = 'Hy'
               else if(l_z(k)) then
                  chid(nh) = 'Hy'
               endif
            endif
         endif
      enddo 

ccc   now just put the E's in order at the end
      do k = 1,nch
         if(lefield(k)) then
            nh = nh + 1
            iordinv(nh) = k
            if(l_x(k)) then
               chid(nh) = 'Ex'
            else if(l_y(k)) then
               chid(nh) = 'Ey'
            endif
            ftype(nh) = 'TE'
         endif
      enddo 
      return
      end
