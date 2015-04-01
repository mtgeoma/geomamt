      subroutine setup(nsta,ntape,stnames,arrayid,nd,nch,cbandf,
     & stcor,decl,orient,chid,dr,npts,var,isuse,sta,itmax,
     & nt_tot,ih_tot,ie_tot,ns,theta)

      include 'iosize.inc'
      include '../include/four_byte.inc'

      character*6 chid(ntmx)
      character*1 ans
      character*3 sta(*)
      integer nch(nsta),ntape(nsta),npts(ndmax,*),nd(nsta)
     &,isuse(2,ndmax),itemp(2),nt_tot,ih_tot(*),ie_tot(*),
     &   bytes_per_rec(nstamx),lar_id,iclong
      real dr(ndmax,*),orient(2,*),stcor(2,*),decl(*),var(*)
      character*40 stnames(nsta,*),fname
      character*10 arrayid
      character*80 cbandf

      logical YoN
      YoN(ans) = (ans.eq.'Y').or.(ans.eq.'y')

c     read in name of band set up file
      read(1,'(a80)') cbandf
c      write(*,*) cbandf
c     read in info for each station
      i0 = 1
      do 5 i = 1,nsta
         read(1,*) ntape(i),nch(i)
c         write(*,*) i,ntape(i),nch(i)
         read(1,*) (var(k),k=i0,i0+nch(i)-1)
c         write(*,*)  (var(k),k=i0,i0+nch(i)-1)
         i0 = i0+nch(i)
         do 4 j = 1,ntape(i)
            read(1,'(a40)') stnames(i,j)
4           continue
         read(1,'(a3)') sta(i) 
c         write(*,*) sta(i)
5        continue

      read(1,'(a10)') arrayid

ccc   begining of optional lines in array.cfg file
c     read in rotation angle for output vector coordinate system
ccc       angle is degree EAST OF GEOMAGNETIC NORTH of x-axis
ccc     default is geomag coordinates
      theta = 0.0
      read(1,*,end=7) theta
    
c     read in maximum number of iterations for robust iteration
      itmax = 0
      read(1,*,end=7) itmax

c   read in limits on set numbers to use for each decimation level
c             (optional)
c    read in number of decimation levels
      read(1,*,end=7) ndt
      do 6 i = 1,ndt
         read(1,*,end=7) id,itemp(1),itemp(2)
         isuse(1,id) = itemp(1)
         isuse(2,id) = itemp(2)
6        continue

c      write(*,*) 'isuse = ',isuse

7     continue
      write(*,*) 'done reading array.cfg'

      ii=20
      it = 1
      do 10 i=1,nsta
      if(lpack) then
         iorecl(i) = 4*(nch(i)+1)
      else
         iorecl(i) = 4*(2*nch(i)+1)
      end if
      bytes_per_rec(i) = iorecl(i)
      if(l_4byte) iorecl(i) = iorecl(i)/4
      do j=1,ntape(i)
         ii=ii+1
         inunit(i,j)=ii
         cfilein(i,j) = stnames(i,j)
         open(unit=ii,file=cfilein(i,j),status='old',
     &      form='unformatted',access='direct',recl=iorecl(i))
         print*,'unit   ',inunit(i,j),'   file   ',cfilein(i,j),
     &           '   open  ','rl=',iorecl(i)
      enddo

      call rfhead(inunit(i,1),nchi,ndmax,nd(i),nfmaxi
     &  ,npts(1,i),dr(1,i),chid(it),orient(1,it),decl(i),
     &          stcor(1,i),bytes_per_rec(i))
      it = it + nch(i)

       write(*,*) 'nchi,ndmax,nd(i),nfmaxi,npts(1,i)'
       write(*,*) nchi,ndmax,nd(i),nfmaxi,npts(1,i)

         if(nchi.ne.nch(i)) then
            print*,'warning: nch in set up file disagrees with nch in'
     &                      ,' file header'
            print*,'file = ',cfilein(i,j),' nch in set up file = ' 
     &            ,nch(i), '  nch in header = ',nchi
            return
         end if

         if(nd(i).gt.ndmax) then
          print*,'warning: nd(i) .gt. ndmax','  i = ',i,' nd= ',nd(i),
     &              'ndmax =',ndmax
             return
         end if

         if(nfmaxi.gt.nfreqmax) then
          print*,'warning: nfreqmax in file header exceeds max set',
     &               ' in program'
            print*,'file = ', cfilein(i,j),'nfmaxi =',nfmaxi,
     &            '  nfreqmx = ',nfreqmax
         end if
         close(20)
10       continue
       write(*,*) (chid(k),k=1,10)
        
ccc   set up array ih (pointers to H component for each station)
ccc   and also array ie (pointers to E components)  
      nt_tot = 1
      do i = 1,nsta
         ih_tot(i) = nt_tot
         nt_tot = nt_tot + nch(i)
         do k = 1,nch(i)
            if(chid(ih_tot(i)+k-1)(1:1).eq.'E') then
               ie_tot(i) = ih_tot(i)+k-1
               go to 15
            endif
         enddo
         ie_tot(i) = ih_tot(i)+nch(i)
15       continue
      enddo
      ih_tot(nsta+1) = nt_tot
      ntt = nt_tot - 1 + nsta
      ns = nt_tot*(nt_tot-1)/2
      nt_tot = nt_tot-1

      lar_id = iclong(arrayid,10)
c     file for main array output
      open(unit=3,file=arrayid(1:lar_id)//'.M',status='unknown')

      return
      end
c______________________________________________________________________
c
      subroutine bset(nb,ibandlim,idl,period,nd,nbmax,
     1  ierr,samprate,npts,cfile)
c          standard frequency averaging band setup routine
      integer ibandlim(2,nbmax),idl(nbmax),npts(nd)
      real period(nbmax),samprate(nd)
      character*80 cfile
        
      open(unit=59,file=cfile,status='old')
        
      ierr = 0
      read(59,*) nb
      if(nb.gt.nbmax) then
         print*,'error: nb exceeds nbmax; nb = ',nb,'nbmax=',nbmax
         ierr = 1
         return
      end if
           
      do ib=1,nb
         read(59,*) idl(ib),ibandlim(1,ib),ibandlim(2,ib)
         if((idl(ib).gt.nd).or.(idl(ib).lt.0)) then
           print*,'error: decimation levle invalid; idl=',idl(ib)
           ierr = 2
           return
         end if
           
         period(ib)=samprate(idl(ib))*npts(idl(ib))*2./
     &                            (ibandlim(1,ib)+ibandlim(2,ib))
      enddo

      close(59)

      return
      end
c___________________________________________________________
c
      function iclong(croot,n)
      character*1 croot(n)
      do 10 i = n,1,-1
      if(ichar(croot(i)).ne.32) then
         iclong = i
         return
      end if
10    continue
      iclong = 0
      return
      end
