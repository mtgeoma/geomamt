c      integer MAXS        maximum number of sites
c      integer MAXF        maximum number of frequencies per band
c      integer MAXDAT      maximum number of data frequencies in file
c      integer MAXBND      maximum number of bands
c      integer NREL        maximum number of realizations
c
c      integer MINREL      minimum number of realizations for statistics
c
c      integer LW   # size of real    work array LIWORK in E04UPF
c      integer LIW  # size of integer work array LWORK  in E04UPF
c
c      LIWORK > 3*N + NCLIN + 2*NCNLN + M*(N+3)
c      LWORK > 20*N + M*(N+3)
c
c      N - no. of variables = azim + MAXS(shear) + MAXS(twist) +
c                             MAXS(2 complex impedances*MAXF = 4*MAXF)
c      NCLIN - no. of general linear constraints = 0
c      NCNLN - no. of nonlinear constraints = 0
c      M - no. of subfunctions = 8*MAXS*MAXF
c

      integer MAXS, MAXF, MAXDAT, MAXBND, MAXREL, MINREL

      integer MAXFBD, MAXFS, NVARIABLES, NCLIN, NCNLN, MSUBFUNC, LW, LIW

      parameter ( MAXS   = 25 )
      parameter ( MAXF   = 60 )  
      parameter ( MAXDAT = 60 )  
      parameter ( MAXBND = 18 )
      parameter ( MAXREL =  1 )
      parameter ( MINREL = 11 )

      parameter ( MAXFBD = MAXF*MAXBND )  
      parameter ( MAXFS  = MAXF*MAXS )  

      parameter ( NVARIABLES = 1 + MAXS + MAXS + 4*MAXS*MAXF )
      parameter ( NCLIN = 0 )
      parameter ( NCNLN = 0 )
      parameter ( MSUBFUNC = 8*MAXS*MAXF )
       
      parameter ( LIW = 3*NVARIABLES + NCLIN + 2*NCNLN
     &                   + MSUBFUNC*(NVARIABLES+3) )
      parameter ( LW  = 20*NVARIABLES
     &                   + MSUBFUNC*(NVARIABLES+3) )
