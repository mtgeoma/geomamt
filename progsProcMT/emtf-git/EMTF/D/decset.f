c________________________________________________________________
c
        subroutine decset(              !nd,nwin,olap,idec,nfc,
     &   nwmx,  !nfcmx,fc,ioffr,ioffs,missmx,ndmx,nchmx,npw,npwmx,
     &   lfd,idoff,ierr)        !!nch,nfuse,

        include 'iounits.inc'
        include 'decimate.inc'


c       this routine sets up the parameters which control the decimation
c       program;  nd (number of decimation levels (undecimated data is
c       level one) and nfcmx (maximum number of filter coefficients) are
c       set in the calling program.  other parameters are set in this
c       routine:
c               samprate - array of sampling rates; only the
c                       undecimated sampling rate is entered; others
c                       are computed from decimation factors
c               nwin - array of window lengths (no. of points)
c               olap - array of overlap lengths
c               idec - array of decimation factors; idec(1) corresponds
c                       to the undecimated data and should be = 1
c               nfc,fc - arrays of filter lengths and filter coefficients;
c                       filters are assumed symmetric and only half
c                       of the coefficients are given - e.g. the filter
c                        .25, .5, .25 used before decimating to level i
c                        is specified by setting nfc(i) = 2 and
c                        fc(i,.) = (.5,.25) ; note that the filter
c                       coefficients for level one are not used
c               ioffr - array of offsets for decimated records; the
c                       decimation filter to make level i records is
c                       centered on records satisfying
c                       mod(rec#,idec(1)* . . . *idec(i)) = ioffr(i)
c               ioffs - array of offsets for starting record numbers
c                       of sets; sets at level i must start with record
c                       numbers which satisfy
c                       mod(rec#,idec(1)*...*idec(i)*(nwin(i)-olap(i)))
c                                       = ioffs(i)
c               missmx - array defining criteria for filling in missing
c                       data points; missing data points are filled in
c                       using a running median; at level id missmx(1,id)
c                       gives the number of points in the running median
c                       centered at the missing point; missmx(2,id) gives
c                       the maximum number of missing points among these
c                       which are tolerated - if more points are missing
c                       the center point is not filled and the set is not
c                       used

c                nfuse - array defining number of fourier coefficients
c                      output; at deciomation level id, nfuse(id) fcs
c                      will be saved in the output fc file; by proper
c                      choice of the values of nfuse for each decimation
c                      level it is possible to minimize wasted storage
c                      space resulting from excessive overlap of
c                      frequency ranges for adjacent decimation levels


cnew        real fc(ndmx,nfcmx),olap(nd)

cnew        integer nfc(*),nwin(*),ioffr(*),ioffs(*),npwmx(*)
cnew     1   ,missmx(2,*),idec(*),npw(nchmx,*),nfuse(*)

        logical lfd(nch,nd)

c       reads from file cf_decset; first line gives number of dec.
c       levels and sampling rate;  for other lines: two lines
c       per decimation level (second contains filter coefficients,
c       first everything else; see read stmt.)

c       routine also now reads array nfuse from cf_decset file

c   12-87
c     routine now also opens file cf_pwset and reads in instructions
c     for prewhitening of time series; the input file should contain
c     a filter length for each channel/ decimation level pair.
c     if the filter length = -1 a standard first difference operator
c     will be used to pre-whiten the time series; if filter length
c     is = 0 or +1 no filtering will be done; if fiolter length > 1
c     this defines the length of the ar prewhitening filter.
c
c       at this time changes are also being made to allow the set
c     overlap (given by olap) to be non-integer


c       first read in number of decimation levels
     
        open(unit=dec_unit,file=cfdecset,status = 'old', err=1000)
        read(dec_unit,*) nd,idoff

        ierr = 0
        if(nd.gt.ndmx) then

c       nd exceeds ndmx set in main; return error and abort
           ierr = 1
           print*,'error: # of dec. levels exceeds ndmx'
           return
        end if

        print*
        print*,'Decimation Level Setup from: ',
     &  cfdecset(1:iclong(cfdecset,80))
        print*
        print*,'DLEV   NPOINTS   OLAP   DEC-FAC   NFC'

      do i = 1,nd
         read(dec_unit,*) nwin(i),olap(i),idec(i),ioffr(i),ioffs(i),
     &   missmx(1,i),missmx(2,i),nfuse(i),nfc(i)
         write (6,'(i4,i8,f8.2,i8,i8)')i, nwin(i),olap(i),idec(i),
     &      nfuse(i)
         if(nwin(i).gt.nwmx) then
            ierr = -1
            print*,'nwmx ( = ',nwmx,') not set large enough'
            print*,'decimation level = ',i,' nwin = ',nwin(i)
            return
         end if
         if(nfc(i).gt.nfcmx) then
            ierr = -1
            print*,'nfcmx ( = ',nfcmx,') not set large enough'
            print*,'decimation level = ',i,' nfc = ',nfc(i)
            return
         end if
         read(dec_unit,*) (fc(i,j),j=1,nfc(i))
      enddo 
      close(dec_unit)

c     now get prewhitening info
        open(unit = pw_unit, file = cfpwset,status='old', err= 1001)
        print*
        print*,'Prewhitening settings from: ',
     &   cfpwset(1:iclong(cfpwset,80))
        print*
             
        read(pw_unit,*) nd1

        if(nd1 .ne. nd) then
           ierr = -2
           print*,'nd set incorrectly in cf_pwset'
           return
        end if

ccc     March 20, 1998  The charade ceases: pre-whitening parameters
ccc     are now the same for all channels

        do i = 1,nd
           read(pw_unit,*)  npw(1,i)
           do j = 2,nch
              npw(j,i) = npw(1,i)
           enddo
           do j= 1,nch
              lfd(j,i) = npw(j,i) .eq. -1
              npw(j,i) = max(1,abs(npw(j,i)))
              if(lfd(j,i)) npw(j,i) = 2
           enddo
           npwmx(i) = npw(1,i)
        enddo
        close(pw_unit)

        return
 1000   print*,' ERROR OPENING DECSET FILE !'
        stop
 1001   print*,' ERROR OPENING PWSET FILE !'
        stop
        end
