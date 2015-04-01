c
c****************************************
c
        subroutine dcimte(ix1,nin,xx,nchp1,nmax
     1   ,ifirst,ifirstd,next,lstart,lstartd,nacc,rmsval)

cnew    ,nd,nwin,ioff,idec,nfc,nfcmax,fc
c    this subroutine does the digital low pass filtering and decimating

c      new version 10 - 16 -87; changes: filtered and decimated data
c          is stored as real (instead of integer) array

c   some of the logic of this routine is fairly subtle; fiddle with it at your own risk
c    the input variables (and the functioning of the routine are better described in the
c   documentation of decset

      include 'decimate.inc'
      integer ix1(nchp1,0:nmax),ifirst(nd),next(nd)
     &   ,nmod(ndmx),ifirstd(nd),nacc(nd)
      logical lstart(nd),lstartd(nd)
      real temp(50),xx(nd,nchp1,0:nmax)

      nch = nchp1-1
      nmod(1) = 1
      do i = 2,nd
          nmod(i) = nmod(i-1) * idec(i)
      enddo

ccc   move nin undecimated data points into real array xx
      n = nin - 1
      do i = 0,n
      j = next(1) + i
      j = mod(j,nmax)
         do k = 1,nchp1
            xx(1,k,j) = ix1(k,i)
         enddo
      enddo

ccc   reset next data pointer for decimation level 1
      next(1) = mod( next(1)+nin , nmax )

ccc   level 1 is now ready to process and to decimate to level 2
      id1 = 1
      id2 = 2

ccc   main decimation loop

15    continue

ccc      find number of points accumulated at level id1
         call ptdist(ifirst(id1),next(id1),nmax,nacc(id1),lstart(id1))

c        decimate to level id2
         call ptdist(ifirstd(id1),next(id1),nmax,naccd,lstartd(id1))
         nn = 0
         j2 = next(id2) - 1
         ii1 = ifirstd(id1)
         ii2 = ii1 + naccd - nfc(id2)
         ii3 = idec(id2)

ccc      filtering loop is executed for j's which correspond to center 
ccc      points for filtered decimation

         do j = ii1,ii2,ii3
            j1 = mod(j,nmax)
            j2 = mod(j2+1,nmax)
            nn = nn+1

ccc         copy record number (center of filter) into record number 
ccc         for next higher level of decimation
            xx(id2,1,j2) = xx(id1,1,j1)
        
            if(nfc(id2) .eq. 1) then
               do k=2,nch+1
                  xx(id2,k,j2) = xx(id1,k,j1)
               enddo
            else

ccc         filter and decimate
               l0 = nfc(id2)
               do k = 1,nch
                  temp(l0) = xx(id1,k+1,j1)
                  do l = 2,l0
                     jl1 = mod(j1+l-1,nmax)
                     jl2 = mod(j1-l+1,nmax)
                     jl2 = mod(jl2+nmax,nmax)
                     l1 = l0 + l - 1
                     l2 = l0 - l + 1
                     temp(l1) = xx(id1,k+1,jl1)
                     temp(l2) = xx(id1,k+1,jl2)
                  enddo
                  call dcfilt(temp,rmsval,nfc,fc,nd,ndmx,id2,xd)
                  xx(id2,k+1,j2) = xd
               enddo
            end if
         enddo            !   end of filtering/decimation loop

ccc     reset next data pointers
        ifirstd(id1) = mod(j,nmax)
        next(id2) = mod( next(id2) + nn,nmax )
        id1 = id1 + 1
        id2 = id2 + 1
        if ( id2 .le. nd) go to 15

ccc     else no more levels;
        call ptdist(ifirst(id1),next(id1),nmax,nacc(id1),
     &                     lstart(id1))
        return
        end
