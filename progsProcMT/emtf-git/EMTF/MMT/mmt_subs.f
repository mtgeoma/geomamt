       subroutine diag(s,n,sdiag)
c       extracts real diagonal from complex hermitian s in symetric
c      storage mode
       complex s(*)
       real sdiag(*)

       do i=1,n
          ii = (i*(i+1))/2
          sdiag(i) = real(s(ii))
       enddo
       return
       end
c______________________________________________________________________
c
      subroutine hnorm(u,nt,ih,ie,nsta)
c     takes a complex vector u of length nt and normalizes
c     all components so RMS of total horizontal magnetic
c     components is 1

      complex u(nt)
      integer ih(*),ie(*)

      hscl = 0.0
      do ista = 1,nsta
         do i = ih(ista),ih(ista+1)-1
            hscl = hscl + abs(u(ih(ista)))**2.
          enddo
      enddo
      hscl = sqrt(hscl/nsta)
      do i = 1,nt
         u(i) = u(i)/hscl
      enddo
      return
      end
c______________________________________________________________________
      subroutine eunits(u,nt,ih,ie,nsta,period)

c       takes vectors and scales E-field components to
c        have same units as H-fields .... the scaling convention
c        is so that a E and H will be the same magnitude for
c          an apparent resistivity of 100 ohm-m
c      
      complex u(nt)
      integer ih(*),ie(*)
      escl = sqrt(period/500.)
      do ista = 1,nsta
         do i = ie(ista),ih(ista+1)-1
            u(i) = u(i) * escl
         enddo
      enddo
      return
      end
c______________________________________________________________________         
c
      subroutine cntxtr(lx,nfreq,nsta,na,la,lall,nxtra)
c        counts up possible combinations of extra sub-arrays
c        identified by integer codes for each record in
c        array lx;  nxtra (output) is number of extra
c        sub-arrays, la gives integer codes for each, na
c       is number of records for the subarray
      integer na(*),la(*),lall,i,j,nxtra,nfreq

      logical lx(nsta,*)
      integer ilx

      nxtra = 0
      do 20 i = 1,nfreq
         ilxt = ilx(lx,nsta,i)
         if(ilxt.ne.lall) then
            do 10 j = 1,nxtra
               if(ilxt.eq.la(j)) then
                  na(j) = na(j)+1
                  go to 20
               end if
10             continue
            nxtra = nxtra+1
            la(nxtra) = ilxt
            na(nxtra) = 1
         end if
20       continue
      return
      end
c______________________________________________________________________
c
      function ilx(lx,nsta,i)
      integer ilx,nsta,i,ib2,k
      logical lx(nsta,*)

      ilx  = 0
      ib2 = 1
      do 5 k = 1,nsta
         if(lx(k,i)) ilx = ilx + ib2
         ib2 = ib2*2
5        continue
      return
      end
c______________________________________________________________________
c
      subroutine mmt_init(machep,maxint,isuse,npts,ndmax,nstamx)
      integer maxint,isuse(2,ndmax),npts(ndmax,nstamx)
      real machep,temp,temp1
ccc........ various initializations

C ******* MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING **********
C       THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC
      machep = 1.
      temp = 1.
      do i = 1,100
         temp = temp/10.
         temp1 = 1. - temp
         if (temp1.eq.1) go to 20
         machep = machep/10.
      enddo
20    continue

c ******** MAXINT IS LARGEST INTEGER (USED FOR MISSING VALUE CODE)
      maxint = 2**30 + (2**30 -1)

ccc  Set up to use all data
      do i = 1,ndmax
         isuse(1,i) = 0
         isuse(2,i) = maxint
         do  j = 1,nstamx
            npts(i,j) = 1
         enddo
      enddo

      return
      end
