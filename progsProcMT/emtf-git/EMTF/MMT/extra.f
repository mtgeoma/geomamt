      subroutine init_xtr(nsta,nt_tot,ntape,nch,nd,ih_tot,lx,isuse,
     & xx,wrk,iwrk,cwrk,nxu,la,ivunits,arrayid,sta,ixs,ixf,stc_tot,
     & chid,decl_tot,orient_tot,nbt,nsig)

      include 'nstamx.inc'
      include '../include/four_byte.inc'
      integer nsta,nt_tot,nch(nsta),ntape(nsta),ih_tot(*),
     &    iband(2),iwrk(*),isuse(2,*),na(100),ivunits(*),ih(ntmx)
   
      character*3 sta(*),stau(nstamx)
      character*2 cfnum(10)
      character*6 chid(*)
      character*10 arrayid
      real wrk(*),stc_tot(2,*),stcor(2,nstamx),decl_tot(*),decl(nstamx),
     &   orient_tot(2,*),orient(2,ntmx)
      complex cwrk(*),xx(*)
      logical llx,lx(nsta,*),la(nsta,50),ltemp
      integer nchu(nstamx),ixs(*),ixf(*),nsig,lar_id,iclong
      parameter (fmin = .10,nstamin =2)

      data cfnum/'S0','S1','S2','S3','S4','S5','S6','S7','S8','S9'/

ccc   get one frequency from decimation level 1
      id = 1
      iband(1) = 1
      iband(2) = 1
      llx = .false.

c      write(*,*) 'In init_xtr, isuse = ',(isuse(1,k),isuse(2,k),k=1,nd)
ccc   Make array data records
      call mkrec(nd,id,iband,nsta,ntape,nch,nt_tot,isuse,
     &  nf_tot,xx,llx,lx,ixs,ixf,cwrk,iwrk)

ccc   Make a list of all sub-array configurations
ccc   nx is total number of sub-arrays (including main)
ccc   la(.,ix) is logical mask for stations in sub-array ix
ccc   na(ix) is number of data vectors in this sub-array
      nx = 1
      na(1) = 0
      do ista = 1,nsta
         la(ista,1) = .true.
      enddo 

      do i = 1,nf_tot
         do ix = 1,nx
            ltemp = .true.
            do ista = 1,nsta
               ltemp = ltemp.and.( (la(ista,ix).and.lx(ista,i)) 
     &           .or. (.not.la(ista,ix).and. (.not.lx(ista,i)) ) )
            enddo
            if(ltemp) then
               na(ix) = na(ix) + 1
               go to 10
            endif
         enddo
         nx = nx+1
         na(nx) = 1
         do ista = 1,nsta
            la(ista,nx) = lx(ista,i)
         enddo
10       continue             
      enddo

ccc  Eliminate arrays with too small a fraction of the
ccc    total number of data vectors  (or too few stations)
ccc      nxu is number of sub-arrays to analyze
ccc      la(.,ix) gives the corresponding logical masks
      
      nxu = 1
      do ix = 2,nx
         nstatmp = 0 
         do ista = 1,nsta
            if(la(ista,ix)) nstatmp = nstatmp + 1
         enddo 
         if( (float(na(ix))/float(nf_tot).gt. fmin)
     &           .and. (nstatmp .ge. nstamin)) then
            nxu = nxu + 1
            do ista = 1,nsta
               la(ista,nxu) = la(ista,ix)
            enddo
         endif
      enddo

ccc   Open files for output of incoherent noise variances, eigenvalues
ccc   and eignevectors (the old "S2 files") ... for all sub-arrays
ccc   Output files will be called S0, S1, S2, S3 ... S9
ccc   Loop over sub-array index
      do ix = 1,nxu
         ivunits(ix) = 60 + ix
         nt = 0
ccc      make array info arrays for each sub-array
         ih(1) = 1
         nstau = 0
         do ista = 1,nsta
            if(la(ista,ix)) then
               nstau = nstau + 1
               nch1 = ih_tot(ista+1) - ih_tot(ista)
               nt = nt + nch1
               ih(ista+1) = nt+1
               stau(nstau) = sta(ista)
               nchu(nstau) = nch1
               stcor(1,nstau) = stc_tot(1,ista)
               stcor(2,nstau) = stc_tot(2,ista)
               decl(nstau) = decl_tot(ista)
               ioff = ih_tot(ista) - ih(nstau)
               do i = ih(nstau),ih(nstau+1)-1
                  orient(1,i) = orient_tot(1,i+ioff) 
                  orient(2,i) = orient_tot(2,i+ioff) 
               enddo
            endif
         enddo
ccc      Write out as a fixed record length direct access file ...
ccc      with one header record ... should be readable by matlab.
ccc      One record will then be written for each frequency band
ccc      Each record will contain :
ccc         period, nf + ev(nsig) + var(nt) + u(nt,nsig)
         ns = (nt*(nt+1))/2
         irecl =  4*(2+nt+2*ns)
         if(l_4byte) irecl = irecl/4
         lar_id = iclong(arrayid,10)

c         write(*,*) 'ix,arrayid = ',ix,arrayid
ccc      Open file
         open(unit=ivunits(ix),file=arrayid(1:lar_id)//'.'//cfnum(ix),
     &      form='unformatted',access='direct',recl=irecl)
ccc      Write out header record
         write(ivunits(ix),rec=1) irecl,nbt,nt,nstau,nsig,
     &      (nchu(k),ih(k),stcor(1,k),stcor(2,k),decl(k),sta(k),
     &      k=1,nstau),(orient(1,l),orient(2,l),chid(l),l=1,nt)
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine masks(la,ix,lx,luse,lch,nsta_tot,nsta,ih_tot,ih,
     &   nt_tot,nt,nf_tot,nf,stc_tot,stc,decl_tot,decl)

      logical la(nsta_tot,*),lch(nt_tot),luse(nf_tot),lx(nsta_tot,*)
     & ,ltemp
      integer ih(*),ih_tot(*),nt,nt_tot,nf,nf_tot
      real stc_tot(2,*),stc(2,*),decl_tot(*),decl(*)

ccc     using logical masks for stations in la, and logical masks for
ccc      all data points, constructs logical masks for 
ccc     components, sets corresponding to subarr ix.
ccc        Also calculates number of components nt
ccc       and number of frequencies nf.   Call once for each
ccc      sub-array/ band
       
      nsta = 0
      ih(1) = 1
      nt = 0
      do ista = 1,nsta_tot
         do j = ih_tot(ista),ih_tot(ista+1) - 1
            lch(j) = la(ista,ix)
         enddo

         if(la(ista,ix)) then
            nsta = nsta+1
            stc(1,nsta)  = stc_tot(1,ista)
            stc(2,nsta)  = stc_tot(2,ista)
            decl(nsta) = decl_tot(ista)
            nch = ih_tot(ista+1) - ih_tot(ista)
            ih(nsta+1) =  ih(nsta) + nch
            nt = nt + nch
         endif
      enddo
 
      nf = 0
      do i = 1,nf_tot 
         ltemp = .true.
         do ista = 1,nsta_tot
            ltemp = ltemp.and.( (la(ista,ix).and.lx(ista,i))
     &        .or. (.not.la(ista,ix).and. (.not.lx(ista,i)) ) )
         enddo
         if(ltemp) then
            nf = nf + 1
            luse(i) = .true.
         else
            luse(i) = .false.
         endif
      enddo

      return
      end 
c______________________________________________________________________
c
      subroutine mask_ch(la,ix,lch,nsta_tot,nsta,ih_tot,ih,
     &   nt_tot,stc_tot,stc,decl_tot,decl)
      logical la(nsta_tot,*),lch(nt_tot)
      integer ih(*),ih_tot(*),nt,nt_tot
      real stc_tot(2,*),stc(2,*),decl_tot(*),decl(*)

      nsta = 0
      ih(1) = 1
      nt = 0
      do ista = 1,nsta_tot
         do j = ih_tot(ista),ih_tot(ista+1) - 1
            lch(j) = la(ista,ix)
         enddo

         if(la(ista,ix)) then
            nsta = nsta+1
            stc(1,nsta)  = stc_tot(1,ista)
            stc(2,nsta)  = stc_tot(2,ista)
            decl(nsta) = decl_tot(ista)
            nch = ih_tot(ista+1) - ih_tot(ista)
            ih(nsta+1) =  ih(nsta) + nch
            nt = nt + nch
         endif
      enddo
      return
      end
c______________________________________________________________________
c
      subroutine data_sort(x,nf_tot,nt_tot,nsta_tot,ih_tot,lx,la,nxu,
     &   y,nf,nt,npoint,luse,ixs,ixf)

ccc   ==> x ::  data matrix produced by mkrec
ccc   ==> nf_tot,nt_tot,nsta_tot,ih_tot,lx ::  usual variables, arrays
ccc          describing data matrix
ccc   ==> la,lex ::  masks for sub-arrays; number of sub-arrays
ccc   <== y  :: reordered data matrix - all data for each sub-array
ccc            are contiguous (no missing data gaps); row column order
ccc            is switched relative to input matrix x
ccc   <== npoint,nf,nuse :: y(npoint(ix)) contains first element
ccc             of data matrix for sub-array ix; nf(ix),nuse(ix) give
ccc            array dimensions
      complex x(nt_tot,nf_tot),y(*)
      logical la(nsta_tot,nxu),lx(nsta_tot,*),ltemp,luse(nf_tot)
      integer nt(nxu),nf_tot,nf(nxu),npoint(nxu),iuse(200),ih_tot(*),
     &   ixs(*),ixf(*)

ccc   take data (some missing) in array x, sort by sub-arrays,
ccc   put in order appropriate for QR etc

c      write(*,*) 'nxu,nf_tot,nt_tot,nsta_tot, = ',nxu,
c     &    nf_tot,nt_tot,nsta_tot
c      write(*,*) 'ih_tot',(ih_tot(k),k=1,nsta_tot+1)
      ipoint = 0
ccc   outer loop over sub-arrays
c      open(unit=83,file='xx.new')

      do ix = 1,nxu
ccc      figure out which channels to use in sub-array ix
         nt(ix) = 0
         do ista = 1,nsta_tot
            if(la(ista,ix)) then
               do i=ih_tot(ista),ih_tot(ista+1)-1
                  nt(ix) = nt(ix)+1
                  iuse(nt(ix)) = i
               enddo
            endif
         enddo

ccc         write(*,*) 'la',(la(ista,ix),ista=1,nsta_tot)
ccc         write(*,*) 'nt',nt(ix)
ccc       write(*,*) 'iuse',(iuse(k),k=1,nt(ix))

ccc      see which points belong in sub-array ix, count them
         nf(ix) = 0
         ixp = ipoint/nt(ix)
         do i = 1,nf_tot
            ltemp = .true.
            do ista = 1,nsta_tot
               ltemp = ltemp.and.( (la(ista,ix).and.lx(ista,i))
     &         .or. (.not.la(ista,ix).and. (.not.lx(ista,i)) ) )
            enddo
ccc         (i.e., ltemp is true if logical mask for point i agrees with
ccc            logical mask for array ix) 
            if(ltemp) then
               luse(i) = .true. 
               nf(ix) = nf(ix) + 1
               ixs(ixp+nf(ix)) = ixs(i)
               ixf(ixp+nf(ix)) = ixf(i)
            else
               luse(i) = .false.
            endif
         enddo
         
ccc      now move them from array x to array y
         npoint(ix) = ipoint + 1
         print*,'nf(ix)',nf(ix),';  nt(ix)',nt(ix)
         id = 0
         do i = 1,nf_tot
            if(luse(i)) then
ccc            move data  
               id = id + 1
               do j = 1,nt(ix) 
                  ij = ipoint + nf(ix)*(j-1) + id
                  y(ij) = x(iuse(j),i)
               enddo
            endif
         enddo
         ipoint = ipoint + nt(ix)*nf(ix)
      enddo
      return
      end
