c>>>>>>>>>>>>>>>>>>>>>>>>>>>>    parameters to set   <<<<<<<<<<<<<<<<<<<<<

ccc    ARRAY VERSION
c       nstamx = maximum # of stations
c       nchmx = # maximum of channels of data
      include 'nstamx.h'
c       ntfmax = maximum number of (complex) data points to allow for
c                for a single transfer function band
c       ncbmx = maximum number of data points to allow for in
c               wide coherence sort band
c       nsetmx = maximum number of sets to allow for
      parameter (ntfmax = 100000, ncbmx = ntfmax,nsetmx=ntfmax/5)

c       ntfmax = maximum total number of frequencies per set in fourier
c               coefficient file
c       ndmax = number of decimation levels in fourier coefficient file
c       ntpmax = maximum number of tapes to stack

      parameter (nfreqmax = 64, ndmax = 8,ntpmax = 30)
      parameter (nsmx=ntmx*(ntmx+1)/2,ntmx2 = ntmx*2+1)
      parameter (nhdrec = 20)

c          logical variables which control io functioning

c      lpack = .true. for packed integers
      logical lpack
      parameter (lpack = .true.)

ccc   ALL of this is commented out for array version
ccc     For array processing, llx is a variable declared
c      llx = .true. and nlx = ncbmx
c             to return all array records in band, with logical
c      array lx to mark which stations have good data
c      logical llx
c      integer nlx
c     llx = .false. and nlx = 1 to only return array records for which
c           all data is present
c      parameter (llx = .false.,nlx = 1)
c      parameter (llx = .true.,nlx = ncbmx)

c       lfop = .true. to open and close FC files for each read of a
c        frequency block; used when there are too many stations to allow
c        all units to be open at the same time
      logical lfop
      parameter (lfop = .false.)

c      lfull = .true. for full info return by mkrec : set numbers and
c       frequency band number for each array record; (if .false. just
c       return complex FCs - all that is needed for routine computations).
c        if lfull is true set nfull = ncbmx; else = 1)
      logical lfull
      parameter (lfull = .true.,nfull = ncbmx)

c       ljunk(ista) is true to use records for staton ista which are output in FC file
c       with negative set numbers; use with ljunk(1) = .false., ljunk(2) = .true.
c       to allow good mag field FCs to be used for RR processing, even when
c       corresponding E-fields are crap
      logical ljunk(nstamx)

      integer inunit(nstamx,ntpmax),irecd(ndmax,nfreqmax,nstamx,ntpmax),
     &  iorecl(nstamx),ixs(nfull),ixf(nfull)
      character*80 cfilein(nstamx,ntpmax)
      common /iobolck/inunit,cfilein,iorecl,irecd,ljunk
