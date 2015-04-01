      program sys_tbl_tf

      include 'fcor.inc'
      integer nch,nfil(nchmx),iftype(nfilmax,nchmx)
      real sampr,sc(nchmx),decl,stcor(2),orient(2,nchmx),cda,cdb
      logical lclkd
      character*1 chid(nchmx),ans
      character*40 cf
      character*80 afparam(nfilmax,nchmx),cdir,cfsp,cfsTF

ccc   get directory path for system parameter (input) and table (output) files
ccc   (same directory for both ... also same directory is used throughout
ccc   the run ... run separately for each directory)
      print*,'directory for sp____ and sTF______'
      read(5,'(a80)') cdir
      do ndir = 80,1,-1
         if(cdir(ndir:ndir).ne.' ') go to 5
      enddo
      ndir = 0
5     continue
      if(ndir.gt.0) then
         cfsp(1:ndir) = cdir(1:ndir)
         cfsTF(1:ndir) = cdir(1:ndir)
         cfsp(ndir+1:ndir+3) = '/sp'
         cfsTF(ndir+1:ndir+4) = '/sTF'
         isp0 = 3
         isTF0 = 4
      else
         cfsp(1:2) = 'sp'
         cfsTF(1:3) = 'sTF'
         isp0 = 2
         isTF0 = 3
      endif

c       return here to start processing another data file
1000  continue
      print*,'data file name (d_****####)'
      read(5,'(a40)') cf
      do ncf = 40,1,-1
         if(cf(ncf:ncf).ne.' ') go to 10
      enddo
      go to 1000
10    continue
      do i = 3,ncf
         cfsp(ndir+isp0+i-2:ndir+isp0+i-2) = cf(i:i)
         cfsTF(ndir+isTF0+i-2:ndir+isTF0+i-2) = cf(i:i)
      enddo
      do i = ndir+ncf+isp0+1,80
         cfsp(i:i) = ' '
      enddo
      do i = ndir+ncf+isTF0+1,80
         cfsTF(i:i) = ' '
      enddo

c  read system parameter file:  number of channels,
c  electrode line lengths, filter parameters, conversion
c  factors, etc. from file spsta###x
      iounit=3
      open(unit=iounit,file=cfsp)
      call getsp1(nch,sampr,sc,nfil,iftype,afparam,
     &   iounit,decl,stcor,orient,chid,cda,cdb,lclkd)

ccc   make table with only system corrections
ccc   and output to file sTF***###
      call systblsu(sampr,nch,nfil,afparam,iftype,sc,cda,cfsTF)

      print*,'Another file?'
      read(5,'(a1)') ans
      if(ans.eq.'y') go to 1000
      stop
      end
