c    include file IOUNITS.INC
c    
c    define all I/O units here to avoid opening files 
c    under the same unit. Set all units using PARMATER
c    statements.
c
c    This file is included in the folowing routines:
c        dnff
c        cininit        (inpu_bin.f)
c        rdhd        (inpu_bin.f)
c        rdblk        (inpu_bin.f)
c        decset
c        getsp
c        bdrscu
c        outinit        (out_pack.f)
c        resptbl

        integer         br_unit
        integer         pw_unit
        integer         dec_unit
        integer         in_unit
        integer         out_unit
        integer         rsp_unit
        integer         sp_unit
        integer         pth_unit
        integer         clk_unit

        character*80    cfbr            ! bad record file
        character*80    cfdecset        ! decimation level definition file
        character*80    cfpwset         ! pre-whitening definition file
        character*80    cfout           ! FC output file
        character*80    cfsp            ! System parameter file
        character*80    cfrsp
        character*80    cfclk           ! clock reset file (ASCII data files)

        parameter (br_unit = 69)
        parameter (pw_unit = 1)
        parameter (dec_unit = 1)
        parameter (in_unit = 12)
        parameter (out_unit = 15)
        parameter (rsp_unit = 37)
        parameter (sp_unit = 55)
        parameter (pth_unit = 77)
        parameter (clk_unit = 76)

        common /FILES/  cfbr,cfdecset,cfpwset,cfout,cfsp, cfrsp, cfclk
