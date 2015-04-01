C        include file input.inc
C        contains variables used in input routines
C
C        routines including input.inc:
C                dnff    (dnff.f)
C                cininit (inpu_bin.f)
C                rdhd    (inpu_bin.f)
C                rdblk   (inpu_bin.f)



        integer         time_int
        integer         iclk0
        real            sampfreq
        real            dr
        logical         lfmiss
        logical         lfirst
        character*80    cfile
        integer         bytes
        integer         itf

ccc     variables added for ASCII input data file
        logical         l_asc
ccc                          = .true. for ASCII data file
        logical         l_clk_hd
ccc                          = .true. for clock reset info in file header
        logical         l_asc_rec_num
ccc                          = .true. for sample numbers at beginning
ccc                             of each data record



      common /inblk/  dr,iclk0,ioff,nmsmx,lfirst,lfmiss,cfile,time_int,
     &     sampfreq,bytes, itf
ccc     variables added for ASCII input data file
     &      , l_clk_hd,l_asc,l_asc_rec_num
