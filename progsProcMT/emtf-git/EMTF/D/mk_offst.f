      subroutine mk_offset(time_int,iset_offset,   !ioffs,ioffr,nd,
     &    sampfreq,dt)  !idec,nwin,olap,

      include 'decimate.inc'

ccc   Computes offsets to add to set numbers computed using run start time
ccc   so that output set numbers are given relative to survey start time.
ccc   Also modifies ioffs so that sets for runs at different stations
ccc   are aligned even if run start times differ

      integer time_int,iset_offset(nd)  !,nd,ioffs(nd),nwin(nd)
cnew     &   idec(nd),ioffr(nd)
      real sampfreq,dt  !olap(nd),

ccc   input: time_int       =  time in seconds of run start time,
ccc                             relative to survey start time
ccc          nd             =  number of decimation levels
ccc          idec(nd)       =  decimation factors
ccc          nwind(nd)      =  window lengths
ccc          olap(nd)       =  overlap between adjacent windows
ccc          ioffs(nd)      =  relative starting position for sets
ccc                              at each decimation level
ccc          sampfreq       =  sampling frequency (real)  in hz

ccc   output: ioffs(nd)     =  relative starting position for sets
ccc                              at each decimation level (modified)
ccc         iset_offset(nd) = offset to add to set numbers computed
ccc                            for each decimation level
ccc         dt              = fractional sample time shift, which may be
ccc                           required when fractional sampling rates are used
ccc                            (e.g., 3.125 hz) ; used to correct phase of FCs

      real*8 start_samp,set_num,rem
      integer npts_set,id,k

ccc   compute starting sample number of run start time (double precision real)
      start_samp = dble(time_int)*dble(sampfreq)
c      write(0,*) 'start_samp',start_samp
      nmod = 1
      do id = 1,nd
ccc      calculate npts_set = number of points between set starting points
         npts_set = nwin(id)-olap(id)
c          write(0,*) npts_set
         do k = 1,id
             npts_set = npts_set*idec(k)
         enddo
         ioffs(id) = mod(ioffs(id) + npts_set,npts_set)
c         write(0,*) 'id,npts_set,ioffs(id)',id,npts_set,ioffs(id)
         if(ioffs(id) .gt. 0 ) then
            set_num = 1 + (start_samp-ioffs(id))/dble(npts_set)
         else
            set_num = start_samp/dble(npts_set)
         endif
c         write(0,*) 'set_num',set_num
         iset_offset(id) = nint(set_num+.4999999999999)
         rem =  (iset_offset(id) - set_num)*npts_set
c         write(0,*) 'rem',rem
         ioffs(id) = nint(rem)
         nmod = nmod*idec(id)
         ioffr(id) = mod(ioffs(id),nmod)
c         write(0,*) 'id,ioffr(id),ioffs(id)',id,ioffr(id),ioffs(id)
c         write(0,*) 'iset_offset(id),set_num,npts_set,rem',
c     &           iset_offset(id),set_num,npts_set,rem
         if(id.eq.1) dt = ioffs(1) - rem
c         write(0,*) dt
       enddo
       return
       end
