C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE SearchLM(LOGFILE_SCREEN,
     >   REQUIRED_RMS,SD_RMS,CONT_HIGHER_RMS,
     >   MIN_LM,MAX_LM,STARTING_LM,FIX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >   MAX_SEARCH_LM,
     >   PARABOLIC_CORRECTION,DoSmooth,dmode,
     >   NMode,NRes,NPer,NSta,ModTyp,ResTyp,Period,DatRes,DatErr,Cd,
     >   NN,NNT,DatInx,SenInx,StsInx,StaNod,SSIndx,
     >   DFStatus,PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >   dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,PrevTol,
     >   ETOL,MAX_PCG_ITER,FlagMu,
     >   x,fx,rx,RhoX,FFX,SSX,RMS,RMSS,RMSP)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER FlagMu,DoSmooth,MAX_SEARCH_LM,MAX_PCG_ITER,
     >        CONT_HIGHER_RMS,LOGFILE_SCREEN
      REAL*8  REQUIRED_RMS,SD_RMS,MIN_LM,MAX_LM,STARTING_LM,
     >        STEPSIZE_LM,SMOOTH_SZLM,
     >        PrevTol,ETOL,PARABOLIC_CORRECTION 
      INTEGER Nza,Nzb,Nz,Ny,FlagSearch,FIX_LM
      INTEGER dmode,NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER ModTyp(*),ResTyp(NMODMX,NRESMX),NN(NMODMX,NRESMX,2)
      INTEGER DFStatus(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  Dzb(*),Dz(*),Dy(*),Czb(*),Cz(*),Cy(*)
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX),Cd(*)
      REAL*8  BB(NN0MX,LL0MX)
      REAL*8  PriRho(NZ0MX,NY0MX)
      REAL*8  CmHt(MM0MX,LL0MX),CDhat(*),DD(*),Repm(*)
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)
      REAL*8  RhoX(NZ0MX,NY0MX),FFX(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  x,fx,rx,a,fa,b,fb,SSX(NMODMX,NSTAMX)
      REAL*8  RMS(*),RMSS(NMODMX,NRESMX,NSTAMX) 
      REAL*8  RMSP(NMODMX,NRESMX,NPERMX)
 
      REAL*8  maxtolreq

      CALL SimpleSearch(LOGFILE_SCREEN,
     >     REQUIRED_RMS,SD_RMS,
     >     MIN_LM,MAX_LM,STARTING_LM,FIX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >     MAX_SEARCH_LM,
     >     PARABOLIC_CORRECTION,DoSmooth,dmode,
     >     NMode,NRes,NPer,NSta,ModTyp,ResTyp,Period,DatRes,DatErr,
     >     Cd,NN,NNT,DatInx,SenInx,StsInx,StaNod,SSIndx,
     >     DFStatus,PriRho,BB,Repm,CDhat,DD,
     >     CmHtIndx,CmHt,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >     ETOL,MAX_PCG_ITER,PrevTol,
     >     x,fx,rx,RhoX,FFX,SSX,a,fa,b,fb,FlagSearch)


      maxtolreq = REQUIRED_RMS + SD_RMS

1000  CONTINUE

C     CASE 1 : (NEWTOL < PREVTOL)  AND (NEWTOL <= TOLREQ)
C     SOLVE1 : Intercept TOLREQ by using Bisection method
C            : only when FlagSearch = 1

C     FlagSearch = 1, found misfit below TOLREQ
      IF (FlagSearch.EQ.1) THEN
        DoSmooth = 1
        CALL BisectionSearch(LOGFILE_SCREEN,
     >  REQUIRED_RMS,SD_RMS,DoSmooth,dmode,
     >  NMode,NRes,NPer,NSta,ModTyp,ResTyp,Period,DatRes,DatErr,
     >  Cd,NN,NNT,DatInx,SenInx,StsInx,StaNod,SSIndx,
     >  DFStatus,PriRho,BB,Repm,CDhat,DD,
     >  CmHtIndx,CmHt,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >  ETOL,MAX_PCG_ITER,
     >  x,fx,rx,RhoX,FFX,SSX,a,fa,b,fb)
        FlagMu      = 1
        STARTING_LM = x
        GOTO 5000
      ENDIF
C     FlagSearch = 2, fx = TOLREQ
      IF (FlagSearch.EQ.2) THEN
        DoSmooth    = 1
        FlagMu      = 1
        STARTING_LM = x
        GOTO 5000
      ENDIF
C     FlagSearch = 3, fb = TOLREQ
      IF (FlagSearch.EQ.3) THEN
        DoSmooth    = 1
        FlagMu      = 1
        STARTING_LM = x
        GOTO 5000
      ENDIF
C     FlagSearch = 4, fa = TOLREQ
      IF (FlagSearch.EQ.4) THEN
        DoSmooth    = 1
        FlagMu      = 1
        STARTING_LM = x
        GOTO 5000
      ENDIF

C     CASE 2 : (NEWTOL <= PREVTOL) AND (NEWTOL > TOLREQ)
C     SOLVE2 : Going to next iteration
      IF ((fx.LT.PrevTol).AND.(fx.GT.maxtolreq)) THEN
        FlagMu      = 2
        STARTING_LM = x
        GOTO 5000
      ENDIF

C     CASE 3 : (NEWTOL > PREVTOL) AND (NEWTOL > TOLREQ)
C     SOLVE2 : Going to next iteration
      IF ((fx.GT.PrevTol).AND.(fx.GT.maxtolreq)) THEN
        IF (CONT_HIGHER_RMS.EQ.1) THEN
          FlagMu      = 3
          STARTING_LM = x
          GOTO 5000
        ELSE
          IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8000) 
          WRITE(99,8000) 
          STOP
        ENDIF
      ENDIF
 
5000  CONTINUE
      CALL CompMisfitAll(NMode,NRes,NPer,NSta,ResTyp,NNT,NN,DatRes,
     >     DatErr,DatInx,FFX,RMS,RMSS,RMSP)

8000  FORMAT('HIGHER RMS IS FOUND, PROGRAM WILL STOP')

      RETURN
      END ! SearchMu

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CopyRespond(NMode,NRes,NPer,NSta,FFO,FFN)
      INCLUDE 'parameter.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*)
      REAL*8  FFO(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  FFN(NMODMX,NRESMX,NPERMX,NSTAMX)
 
      INTEGER im,ir,ip,is
 
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im)
              FFN(im,ir,ip,is) = FFO(im,ir,ip,is)
            ENDDO ! do is 
          ENDDO ! do ip 
        ENDDO ! do ir 
      ENDDO ! do im 
 
      RETURN
      END ! CopyRespond()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CopyDistortion(NMode,NSta,SSO,SSN)
      INCLUDE 'parameter.h'
      INTEGER NMode,NSta(*)
      REAL*8  SSO(NMODMX,NSTAMX),SSN(NMODMX,NSTAMX)
 
      INTEGER is,im
 
      DO im = 1,NMode
        DO is = 1,NSta(im)
          SSN(im,is) = SSO(im,is)
        ENDDO ! is 
      ENDDO ! im
 
      RETURN
      END ! CopyDistortion

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE Replace(Nzb,Ny,xo,fo,ro,Rold,xn,fn,rn,Rnew)
      INCLUDE 'parameter.h'
      INTEGER Nzb,Ny
      REAL*8  xn,fn,rn,xo,fo,ro
      REAL*8  Rnew(NZ0MX,NY0MX),Rold(NZ0MX,NY0MX)

      xn = xo
      fn = fo
      rn = ro
      CALL CopyMatrixR8(1,Nzb,1,Ny,NZ0MX,NY0MX,Rold,
     >                  1,Nzb,1,Ny,NZ0MX,NY0MX,Rnew)

      RETURN
      END ! Replace

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

C     Search for lm to get the minimum misfit

C     The minimum misfit is found under the condition of ...
C     Flag : 0 found misfit less than required tolerance
C          : 1 number of searching iteration exceeds maximum
C          : 2 f(a), f(x) and f(b) are colinear (less than
C              SIGMISFIT range) shold not continue
C          : 3 Precision reach when  a-x = b-x = SIGMUE
C          : < 0  (ERR0R)

      SUBROUTINE SimpleSearch(LOGFILE_SCREEN,
     >    REQUIRED_RMS,SD_RMS,
     >    MIN_LM,MAX_LM,STARTING_LM,FIX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >    MAX_SEARCH_LM,
     >    PARABOLIC_CORRECTION,DoSmooth,dmode,
     >    NMode,NRes,NPer,NSta,ModTyp,ResTyp,Period,DatRes,DatErr,
     >    Cd,NN,NNT,DatInx,SenInx,StsInx,StaNod,SSIndx,
     >    DFStatus,PriRho,BB,Repm,CDhat,DD,
     >    CmHtIndx,CmHt,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >    ETOL,MAX_PCG_ITER,oldfx,
     >    x,fx,rx,RhoX,FFX,SSX,a,fa,b,fb,FlagSearch)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER DoSmooth,MAX_SEARCH_LM,MAX_PCG_ITER,
     >        LOGFILE_SCREEN,FIX_LM
      REAL*8  REQUIRED_RMS,SD_RMS,MIN_LM,MAX_LM,STARTING_LM,ETOL,
     >        STEPSIZE_LM,SMOOTH_SZLM,PARABOLIC_CORRECTION
      INTEGER Nza,Nzb,Nz,Ny,FlagSearch
      INTEGER dmode,NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER ModTyp(*),ResTyp(NMODMX,NRESMX)
      INTEGER DFStatus(*),NN(NMODMX,NRESMX,2)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  Dzb(*),Dz(*),Dy(*),Czb(*),Cz(*),Cy(*)
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX),Cd(*)
      REAL*8  BB(NN0MX,LL0MX),PriRho(NZ0MX,NY0MX)
      REAL*8  CmHt(MM0MX,LL0MX),CDhat(*),DD(*),Repm(*)
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4),oldfx

      INTEGER schno,FWD_COND,im,is
      REAL*8  mintolreq,maxtolreq,tolreq,tol1
      REAL*8  lma,lmx,lmb,lmu
      REAL*8  a,fa,ra,b,fb,rb,x,fx,rx,u,fu,ru

      REAL*8  RhoA(NZ0MX,NY0MX),FFA(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  RhoB(NZ0MX,NY0MX),FFB(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  RhoX(NZ0MX,NY0MX),FFX(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  RhoU(NZ0MX,NY0MX),FFU(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  SSX(NMODMX,NSTAMX),SSU(NMODMX,NSTAMX)
      REAL*8  SSA(NMODMX,NSTAMX),SSB(NMODMX,NSTAMX)

      INTEGER sdir,minflag,maxflag
      REAL*8  xa,xa2,xb,xb2,fxb,fxa,mm,fch,fchdesired


      FlagSearch = 0
      minflag    = 0
      maxflag    = 0

      tolreq    = REQUIRED_RMS
      mintolreq = tolreq - SD_RMS
      maxtolreq = tolreq + SD_RMS
      IF (DoSmooth.EQ.0) tol1 = STEPSIZE_LM
      IF (DoSmooth.EQ.1) tol1 = SMOOTH_SZLM

C     Initialize SSA and SSB
      DO im = 1, NMode
        DO is = 1, NSta(im)
          SSA(im,is) = D0
          SSB(im,is) = D0
        ENDDO
      ENDDO
      

50    CONTINUE

      a  = STARTING_LM-tol1
      IF (a.LT.MIN_LM) THEN
        IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8300)  MIN_LM
        WRITE(99,8300) MIN_LM
        a = MIN_LM
        minflag = 1
      ENDIF
      b  = STARTING_LM+tol1
      IF (b.GT.MAX_LM) THEN
        IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8310)  MAX_LM
        WRITE(99,8310) MAX_LM
        b = MAX_LM
        maxflag = 1
      ENDIF
      fa = D0
      fb = D0
      ra = D0
      rb = D0
      x  = STARTING_LM
      fx = D0
      rx = D0

C     for fixed lm
      IF (FIX_LM.EQ.1) THEN
        lmx = D10**x
        CALL CompModel(dmode,lmx,Nzb,Ny,NMode,NRes,NPer,NSta,
     >       Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >       PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >       RhoX,rx)
        CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >       Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >       SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >       Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >       RhoX,SSX,FFX,fx)
        IF (FWD_COND.EQ.1) THEN
          STARTING_LM = x + tol1
          GOTO 50
        ENDIF

C       improve misfit  fx < oldfx
        IF ((fx.LT.oldfx).OR.(fx.LE.maxtolreq)) THEN 
          IF (LOGFILE_SCREEN.EQ.1) THEN
             WRITE(6,8210) x,fx,rx
          ENDIF
          WRITE(99,8210) x,fx,rx
          IF (fx.GT.maxtolreq) THEN
            GOTO 2000
          ENDIF
          IF (fx.LE.maxtolreq) THEN
            IF ((fx.LE.maxtolreq).AND.(fx.GE.mintolreq)) THEN
               FlagSearch = 2
               GOTO 2000
            ENDIF

            lma = D10**a
            CALL CompModel(dmode,lma,Nzb,Ny,NMode,NRes,NPer,NSta,
     >         Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >         PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >         RhoA,ra)
            CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >         Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >         SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >         Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >         RhoA,SSA,FFA,fa)
            lmb = D10**b
            CALL CompModel(dmode,lmb,Nzb,Ny,NMode,NRes,NPer,NSta,
     >         Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >         PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >         RhoB,rb)
            CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >         Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >         SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >         Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >         RhoB,SSB,FFB,fb)
            GOTO 95
          ENDIF
        ENDIF

C       not improve misfit, search for LM
        IF (fx.GT.oldfx) THEN
          IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8220)
          WRITE(99,8220)
          lma = D10**a
          CALL CompModel(dmode,lma,Nzb,Ny,NMode,NRes,NPer,NSta,
     >         Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >         PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >         RhoA,ra)
          CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >         Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >         SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >         Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >         RhoA,SSA,FFA,fa)
          lmb = D10**b
          CALL CompModel(dmode,lmb,Nzb,Ny,NMode,NRes,NPer,NSta,
     >         Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >         PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >         RhoB,rb)
          CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >         Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >         SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >         Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >         RhoB,SSB,FFB,fb)
          GOTO 95
        ENDIF
      ENDIF


C     If not fixed lm.

      lmx = D10**x
      lma = D10**a
      lmb = D10**b
      CALL CompModel3(dmode,lma,lmx,lmb,
     >     Nzb,Ny,NMode,NRes,NPer,NSta,Period,Cd,NNT,DatInx,SenInx,
     >     DFStatus,PriRho,BB,Repm,CDhat,DD,
     >     CmHtIndx,CmHt,RhoA,ra,RhoX,rx,RhoB,rb)
      CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >     SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >     Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >     RhoA,SSA,FFA,fa)
           IF (FWD_COND.EQ.1) THEN
             STARTING_LM = x + tol1
             GOTO 50
           ENDIF
      CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >     SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >     Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >     RhoX,SSX,FFX,fx)
           IF (FWD_COND.EQ.1) THEN
             STARTING_LM = x + tol1
             GOTO 50
           ENDIF
      CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >     SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >     Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >     RhoB,SSB,FFB,fb)
           IF (FWD_COND.EQ.1) THEN
             STARTING_LM = x + tol1
             GOTO 50
           ENDIF



95    CONTINUE
      schno = 0
100   CONTINUE
      schno = schno + 1

      IF (LOGFILE_SCREEN.EQ.1) THEN
        IF (schno.EQ.1) WRITE(6,8000)
        WRITE(6,8001) schno,a,x,b
        WRITE(6,8002) ' ',fa,fx,fb
        WRITE(6,8003) ' ',ra,rx,rb
      ENDIF
      IF (schno.EQ.1) WRITE(99,8000)
      WRITE(99,8001) schno,a,x,b
      WRITE(99,8002) ' ',fa,fx,fb
      WRITE(99,8003) ' ',ra,rx,rb

      IF (schno.GT.MAX_SEARCH_LM) THEN
       IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8400) 
       WRITE(99,8400)
       IF (fa.LT.fx) THEN
         CALL Replace(Nzb,Ny,a,fa,ra,RhoA,x,fx,rx,RhoX)
         CALL CopyRespond(NMode,NRes,NPer,NSta,FFA,FFX)
         CALL CopyDistortion(NMode,NSta,SSA,SSX)
         FlagSearch = -1
         GOTO 1000
       ENDIF
       IF (fb.LT.fx) THEN
         CALL Replace(Nzb,Ny,b,fb,rb,RhoB,x,fx,rx,RhoX)
         CALL CopyRespond(NMode,NRes,NPer,NSta,FFB,FFX)
         CALL CopyDistortion(NMode,NSta,SSB,SSX)
         FlagSearch = -1
         GOTO 1000
       ENDIF
      ENDIF

C     tend to go for right side (higher value of lm) for lower model norm
      IF (((fx.LE.maxtolreq).OR.(fa.LE.maxtolreq)).OR.
     >    (fb.LE.maxtolreq)) THEN
       IF ((fx.LE.maxtolreq).AND.(fx.GE.mintolreq)) THEN
         FlagSearch = 2
         GOTO 2000
       ENDIF
       IF ((fb.LE.maxtolreq).AND.(fb.GE.mintolreq)) THEN
         CALL Replace(Nzb,Ny,b,fb,rb,RhoB,x,fx,rx,RhoX)
         CALL CopyRespond(NMode,NRes,NPer,NSta,FFB,FFX)
         CALL CopyDistortion(NMode,NSta,SSB,SSX)
         FlagSearch = 3
         GOTO 2000
       ENDIF
       IF ((fa.LE.maxtolreq).AND.(fa.GE.mintolreq)) THEN
         CALL Replace(Nzb,Ny,a,fa,ra,RhoA,x,fx,rx,RhoX)
         CALL CopyRespond(NMode,NRes,NPer,NSta,FFA,FFX)
         CALL CopyDistortion(NMode,NSta,SSA,SSX)
         FlagSearch = 4
         GOTO 2000
       ENDIF

       FlagSearch = 1
       GOTO 1000
      ENDIF


C     a high value in the middle : need new ranges
      IF ((fx.GT.fa).AND.(fx.GT.fb)) THEN
        sdir = 0
        IF (fa.LT.fb) THEN
          GOTO 500
        ENDIF
        IF (fb.LT.fa) THEN
          GOTO 600
        ENDIF
      ENDIF

      IF ((fx.GT.fa).AND.(fx.LT.fb)) THEN
C       minimum is on the left side
        GOTO 500
      ENDIF

      IF ((fx.GT.fb).AND.(fx.LT.fa)) THEN
C       minimum is on the right side
        GOTO 600
      ENDIF

      IF ((fx.LT.fa).AND.(fx.LT.fb)) THEN
C       minimum is in between [a,b]
        GOTO 700
      ENDIF

      GOTO 1000

500   CONTINUE
      sdir = 2
      CALL Replace(Nzb,Ny,x,fx,rx,RhoX,b,fb,rb,RhoB)
      CALL CopyRespond(NMode,NRes,NPer,NSta,FFX,FFB)
      CALL CopyDistortion(NMode,NSta,SSX,SSB)
      CALL Replace(Nzb,Ny,a,fa,ra,RhoA,x,fx,rx,RhoX)
      CALL CopyRespond(NMode,NRes,NPer,NSta,FFA,FFX)
      CALL CopyDistortion(NMode,NSta,SSA,SSX)
      a  = a - tol1

      IF (a.LT.MIN_LM) THEN
        IF (minflag.EQ.1)  GOTO 1000
        IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8300)  MIN_LM
        WRITE(99,8300) MIN_LM
        a = MIN_LM
        minflag = 1
      ENDIF
      fa = D0
      ra = D0
      lma= D10**a 
      CALL CompModel(dmode,lma,Nzb,Ny,NMode,NRes,NPer,NSta,
     >     Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >     PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >     RhoA,ra)
      CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >     SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >     Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >     RhoA,SSA,FFA,fa)
      GOTO 100


600   CONTINUE
      sdir = 3
      CALL Replace(Nzb,Ny,x,fx,rx,RhoX,a,fa,ra,RhoA)
      CALL CopyRespond(NMode,NRes,NPer,NSta,FFX,FFA)
      CALL CopyDistortion(NMode,NSta,SSX,SSA)
      CALL Replace(Nzb,Ny,b,fb,rb,RhoB,x,fx,rx,RhoX)
      CALL CopyRespond(NMode,NRes,NPer,NSta,FFB,FFX)
      CALL CopyDistortion(NMode,NSta,SSB,SSX)

      b  = b + tol1
      IF (b.GT.MAX_LM) THEN
        IF (maxflag.EQ.1) GOTO 1000
        IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8310)  MAX_LM
        WRITE(99,8310) MAX_LM
        b = MAX_LM
        maxflag = 1
      ENDIF
      fb = D0
      rb = D0
      lmb= D10**b
      CALL CompModel(dmode,lmb,Nzb,Ny,NMode,NRes,NPer,NSta,
     >     Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >     PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >     RhoB,rb)
      CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >     SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >     Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >     RhoB,SSB,FFB,fb)
      GOTO 100

700   CONTINUE
      sdir = 1

C     calculate the minimum of a parabola through these three points
C     parabola equation  (x-h)^2 = 4p(y-k)
      xa  = x-a
      xa2 = (x-a)**D2
      xb  = x-b
      xb2 = (x-b)**D2
      fxa = fx-fa
      fxb = fx-fb
      u   = x - 0.5*(xa2*fxb-xb2*fxa)/(xa*fxb-xb*fxa) 
      mm  = ((a-u)**D2)/((b-u)**D2)
      fu  = (mm*fb - fa)/(mm-D1) 

      fchdesired = PARABOLIC_CORRECTION
      fch        = (DABS(fu-fx)/fx)*D100
      IF ((fch.GT.fchdesired).AND.(fu.LT.fx)) THEN

       lmu = D10**u
       CALL CompModel(dmode,lmu,Nzb,Ny,NMode,NRes,NPer,NSta,
     >      Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >      PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >      RhoU,ru)
       CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >      Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >      SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >      Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >      RhoU,SSU,FFU,fu)

        IF (fu.GT.fx) THEN
c         IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8150)
c         WRITE(99,8150)
        ELSE
          IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8110) fu
          WRITE(99,8110) fu

          CALL Replace(Nzb,Ny,u,fu,ru,RhoU,x,fx,rx,RhoX)
          CALL CopyRespond(NMode,NRes,NPer,NSta,FFU,FFX)
          CALL CopyDistortion(NMode,NSta,SSU,SSX)
          IF ((fx.LE.maxtolreq).AND.(fx.GE.mintolreq)) THEN
            FlagSearch = 2
            GOTO 2000
          ENDIF
          IF (fx.LE.mintolreq) THEN
            FlagSearch = 1
            GOTO 1000
          ENDIF
        ENDIF
      ENDIF
    
      IF ((fx.GT.oldfx).AND.(fu.LT.oldfx)) THEN
        lmu = D10**u
        CALL CompModel(dmode,lmu,Nzb,Ny,NMode,NRes,NPer,NSta,
     >       Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >       PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >       RhoU,ru)
        CALL CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >       Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >       SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >       Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >       RhoU,SSU,FFU,fu)

        IF (fu.GT.fx) THEN
c         IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8150)
c         WRITE(99,8150)
        ELSE
          IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8110)
          WRITE(99,8110)
          CALL Replace(Nzb,Ny,u,fu,ru,RhoU,x,fx,rx,RhoX)
          CALL CopyRespond(NMode,NRes,NPer,NSta,FFU,FFX)
          CALL CopyDistortion(NMode,NSta,SSU,SSX)
          IF ((fx.LE.maxtolreq).AND.(fx.GE.mintolreq)) THEN
            FlagSearch = 2
            GOTO 2000
          ENDIF
          IF (fx.LE.mintolreq) THEN
            FlagSearch = 1
            GOTO 1000
          ENDIF
        ENDIF
      ENDIF

      IF ((fx.GT.oldfx).AND.(fu.GT.oldfx)) THEN
c       IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8150)
c       WRITE(99,8150)
      ENDIF

1000  CONTINUE

2000  CONTINUE

8000  FORMAT('  # |                 LEFT     MIDDLE    RIGHT')
8001  FORMAT(i3,' |   LOG10    ',3f10.4)
8002  FORMAT(a3,' |    RMS     ',3f10.4)
8003  FORMAT(a3,' | MODEL NORM ',3f10.4)

8110  FORMAT('PARABOLIC INTERPOLATION REDUCES RMS TO ',f10.4)
8150  FORMAT('PARABOLIC INTERPOLATION CAN NOT HELP REDUCING RMS')


8210  FORMAT('LOG10 = ',f10.4,' RMS = ',f10.4,' MODEL NORM = ',f10.4)
8220  FORMAT('RMS IS NOT IMPROVED, SEARCH FOR LGM IS ENFORCED')

8300  FORMAT('    | NO VALUE BEYOND ',f10.4,' (LEFT) IS USED')
8310  FORMAT('    | NO VALUE BEYOND ',f10.4,' (RIGHT) IS USED')

8400  FORMAT('STOP BRACKETING FOR MINIMUM, ',
     >       'MAXIMUM NUMBER OF SEARCHING LGM IS REACHED')
        
      RETURN
      END ! SimpleSearch

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE BisectionSearch(LOGFILE_SCREEN,
     >    REQUIRED_RMS,SD_RMS,DoSmooth,dmode,
     >    NMode,NRes,NPer,NSta,ModTyp,ResTyp,Period,DatRes,DatErr,
     >    Cd,NN,NNT,DatInx,SenInx,StsInx,StaNod,SSIndx,
     >    DFStatus,PriRho,BB,Repm,CDhat,DD,
     >    CmHtIndx,CmHt,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >    ETOL,MAX_PCG_ITER,
     >    x,fx,rx,RhoX,FFX,SSX,a,fa,b,fb)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER DoSmooth,MAX_PCG_ITER,LOGFILE_SCREEN
      REAL*8  REQUIRED_RMS,SD_RMS,ETOL
      INTEGER Nza,Nzb,Nz,Ny
      INTEGER dmode,NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER ModTyp(*),ResTyp(NMODMX,NRESMX)
      INTEGER DFStatus(*),NN(NMODMX,NRESMX,2)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  Dzb(*),Dz(*),Dy(*),Czb(*),Cz(*),Cy(*)
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX),Cd(*)
      REAL*8  BB(NN0MX,LL0MX),PriRho(NZ0MX,NY0MX)
      REAL*8  CmHt(MM0MX,LL0MX),CDhat(*),DD(*),Repm(*)
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)

      REAL*8  mintolreq,maxtolreq,tolreq,lmxit
      REAL*8  a,fa,b,fb,x,fx,rx,u,fu,xit,fxit,rxit

      REAL*8  RhoX(NZ0MX,NY0MX),FFX(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  RhoXit(NZ0MX,NY0MX),FFXit(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  SSX(NMODMX,NSTAMX),SSXit(NMODMX,NSTAMX)

      INTEGER FWD_COND
      REAL*8  xa,xa2,xb,xb2,fxb,fxa,mm,pfocus

      tolreq = REQUIRED_RMS
      mintolreq = tolreq - SD_RMS
      maxtolreq = tolreq + SD_RMS

C     Initialize SSXit 
      DO im = 1, NMode
        DO is = 1, NSta(im)
          SSXit(im,is) = D0
        ENDDO
      ENDDO


C     calculate the minimum of a parabola through these three points
C     parabola equation  (x-h)^2 = 4p(y-k)
C     pfocus = p, h = u,  k = fu
      xa  = x-a
      xa2 = (x-a)**D2
      xb  = x-b
      xb2 = (x-b)**D2
      fxa = fx-fa
      fxb = fx-fb
      u   = x - 0.5*(xa2*fxb-xb2*fxa)/(xa*fxb-xb*fxa) 
      mm  = ((a-u)**D2)/((b-u)**D2)
      fu  = (mm*fb - fa)/(mm-D1) 
      pfocus = (x-u)**D2/(D4*(fx-fu))


      IF ((fx.GT.mintolreq).AND.(fx.LT.maxtolreq)) THEN
        GOTO 2000
      ENDIF

      IF ((fx.LT.mintolreq).OR.(fx.GT.maxtolreq)) THEN
1000    CONTINUE
c       calculate lm that give us the misfit at tolreq     
        xit = DSQRT(D4*pfocus*(tolreq-fu)) + u


        lmxit = D10**xit
        CALL CompModel(dmode,lmxit,Nzb,Ny,NMode,NRes,NPer,NSta,
     >       Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >       PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >       RhoXit,rxit)
        CALL CompSubRespAll(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >       Period,StaNod,SSIndx,DatRes,DatInx,SenInx,StsInx,
     >       dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >       ETOL,MAX_PCG_ITER,FWD_COND,
     >       RhoXit,SSXit,FFXit)
        CALL CompMisfit(NMode,NRes,NPer,NSta,NNT,
     >       DatRes,DatErr,DatInx,FFXit,fxit)

        IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8102)  fxit,xit
        WRITE(99,8102) fxit,xit

        IF ((fxit.GT.mintolreq).AND.(fxit.LT.maxtolreq)) THEN
          CALL Replace(Nzb,Ny,xit,fxit,rxit,RhoXit,x,fx,rx,RhoX)
          CALL CopyRespond(NMode,NRes,NPer,NSta,FFXit,FFX)
          CALL CopyDistortion(NMode,NSta,SSXit,SSX)
          GOTO 2000
        ENDIF
c       IF (fxit.GT.maxtolreq) THEN
c       ENDIF
c       IF (fxit.LT.mintolreq) THEN
c       ENDIF

        IF (xit.LT.a) THEN
          b  = x
          fb = fx
          x  = a
          fx = fa
          a  = xit
          fa = fxit
          GOTO 1500
        ENDIF
        IF ((xit.GT.a).AND.(xit.LT.x)) THEN 
          b  = x
          fb = fx
          x  = xit
          fx = fxit
          GOTO 1500
        ENDIF
        IF ((xit.GT.x).AND.(xit.LT.b)) THEN 
          a  = x
          fa = fx
          x  = xit
          fx = fxit
          GOTO 1500
        ENDIF
        IF (xit.GT.b) THEN
          a  = x
          fa = fx
          x  = b
          fx = fb
          b  = xit
          fb = fxit
          GOTO 1500
        ENDIF

1500    CONTINUE
        xa  = x-a
        xa2 = (x-a)**D2
        xb  = x-b
        xb2 = (x-b)**D2
        fxa = fx-fa
        fxb = fx-fb
        u   = x - 0.5*(xa2*fxb-xb2*fxa)/(xa*fxb-xb*fxa) 
        mm  = ((a-u)**D2)/((b-u)**D2)
        fu  = (mm*fb - fa)/(mm-D1) 
        pfocus = (x-u)**D2/(D4*(fx-fu))
        GOTO 1000
      ENDIF

2000  CONTINUE

8102  FORMAT('PARABOLIC INTERP. ESTIMATES  RMS = ',f10.4,
     >       ' FROM X = ',f10.4)


      RETURN
      END !  BisectionSearch

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompModel(dmode,lm,Nzb,Ny,NMode,NRes,NPer,NSta,
     >           Period,Cd,NNT,DatInx,SenInx,DFStatus,
     >           PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >           RhoOut,MNorm)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER dmode,Nzb,Ny,NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER DFStatus(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  lm,Period(NMODMX,NPERMX),Cd(*)
      REAL*8  BB(NN0MX,LL0MX),MNorm
      REAL*8  PriRho(NZ0MX,NY0MX),RhoOut(NZ0MX,NY0MX)
      REAL*8  CmHt(MM0MX,LL0MX),CDhat(*),DD(*),Repm(*)

      INTEGER ii,jj,iy,iz
      REAL*8  xx(NN0MX)
      REAL*8  newrho

C     data subspace inversion mode
      IF (dmode.EQ.0) THEN
        CALL CompModelSubSpace(lm,NNT,BB,Repm,CDhat,xx,MNorm)
      ENDIF

C     data space inversion mode
      IF (dmode.EQ.1) THEN
        CALL CompModelFullSpace(lm,NNT,Repm,CDhat,Cd,xx,MNorm)
      ENDIF
      CALL ConstantMatrixR8(RhoOut,NZ0MX,NY0MX,Nzb,Ny,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj = jj+1
          newrho = D0
          DO ii = 1,NNT(3)
            newrho = newrho + CmHt(jj,ii)*xx(ii)
          ENDDO
          IF (DFStatus(jj).EQ.0) THEN
            newrho = DLOG10(PriRho(iz,iy))
          ELSE
            newrho = newrho + DLOG10(PriRho(iz,iy))
          ENDIF
          RhoOut(iz,iy) = D10**newrho
        ENDDO
      ENDDO ! iy

      RETURN
      END ! CompModel


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompModel3(dmode,lma,lmx,lmb,Nzb,Ny,
     >           NMode,NRes,NPer,NSta,Period,Cd,NNT,DatInx,SenInx,
     >           DFStatus,PriRho,BB,Repm,CDhat,
     >           DD,CmHtIndx,CmHt,
     >           RhoA,MNormA,RhoX,MNormX,RhoB,MNormB)

      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER dmode,Nzb,Ny,NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER DFStatus(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  Period(NMODMX,NPERMX),Cd(*)
      REAL*8  BB(NN0MX,LL0MX),PriRho(NZ0MX,NY0MX)
      REAL*8  CmHt(MM0MX,LL0MX),CDhat(*),DD(*),Repm(*)
      REAL*8  lma,lmx,lmb,MNormA,MNormX,MNormB
      REAL*8  RhoA(NZ0MX,NY0MX),RhoX(NZ0MX,NY0MX),RhoB(NZ0MX,NY0MX)

      INTEGER ii,jj,iy,iz,imd
      REAL*8  lm,MNorm
      REAL*8  newrhox,newrhoa,newrhob
      REAL*8  xa(NN0MX),xx(NN0MX),xb(NN0MX),x1(NN0MX)

      imd = 1
100   CONTINUE
      IF (imd.EQ.1) lm = lmx
      IF (imd.EQ.2) lm = lma
      IF (imd.EQ.3) lm = lmb
      IF (imd.EQ.4) GOTO 105

C     data subspace inversion mode
      IF (dmode.EQ.0) THEN
        CALL CompModelSubSpace(lm,NNT,BB,Repm,CDhat,x1,MNorm)
      ENDIF
C     data space inversion mode
      IF (dmode.EQ.1) THEN
        CALL CompModelFullSpace(lm,NNT,Repm,CDhat,Cd,x1,MNorm)
      ENDIF

      IF (imd.EQ.1) THEN
        DO ii = 1,NNT(3)
          xx(ii) = x1(ii)
        ENDDO
        MNormX = MNorm
      ENDIF
      IF (imd.EQ.2) THEN
        DO ii = 1,NNT(3)
          xa(ii) = x1(ii)
        ENDDO
        MNormA = MNorm
      ENDIF
      IF (imd.EQ.3) THEN
        DO ii = 1,NNT(3)
          xb(ii) = x1(ii)
        ENDDO
        MNormB = MNorm
      ENDIF

      imd = imd+1
      GOTO 100
105   CONTINUE

      CALL ConstantMatrixR8(RhoX,NZ0MX,NY0MX,Nzb,Ny,D0)
      CALL ConstantMatrixR8(RhoA,NZ0MX,NY0MX,Nzb,Ny,D0)
      CALL ConstantMatrixR8(RhoB,NZ0MX,NY0MX,Nzb,Ny,D0)

      jj = 0
      DO iy = 1,Ny
        DO iz = 1,Nzb
          jj = jj+1
          newrhox = D0
          newrhoa = D0
          newrhob = D0
          DO ii = 1,NNT(3)
            newrhox = newrhox + CmHt(jj,ii)*xx(ii)
            newrhoa = newrhoa + CmHt(jj,ii)*xa(ii)
            newrhob = newrhob + CmHt(jj,ii)*xb(ii)
          ENDDO

          IF (DFStatus(jj).EQ.0) THEN
            newrhox = DLOG10(PriRho(iz,iy))
            newrhoa = DLOG10(PriRho(iz,iy))
            newrhob = DLOG10(PriRho(iz,iy))
          ELSE
            newrhox = newrhox + DLOG10(PriRho(iz,iy))
            newrhoa = newrhoa + DLOG10(PriRho(iz,iy))
            newrhob = newrhob + DLOG10(PriRho(iz,iy))
          ENDIF
          RhoX(iz,iy) = D10**newrhox
          RhoA(iz,iy) = D10**newrhoa
          RhoB(iz,iy) = D10**newrhob
        ENDDO
      ENDDO ! iy


      RETURN
      END ! CompModel3

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompModelSubSpace(lm,NNT,BB,Repm,CDhat,x1,MNorm)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER NNT(*)
      REAL*8  lm,BB(NN0MX,LL0MX),CDhat(*),Repm(*),x1(*),MNorm

      INTEGER ii,jj,kk,info
      REAL*8  RR(LLHMX),b1(LL0MX),sumx

C     data subspace inversion mode

      DO ii = 1,(NNT(3)+1)*NNT(3)/2
        RR(ii) = Repm(ii)
      ENDDO
      DO ii = 1,NNT(3)
        jj = ii
        kk = ii+(jj-1)*jj/2
        RR(kk) = lm + Repm(kk)
        b1(ii) = CDhat(ii)
      ENDDO
C     Using Cholesky decomposition to solve 
      CALL DPPSV('U',NNT(3),1,RR,b1,LL0MX,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR WHILE COMPUTING MODEL !!!'
        WRITE(6,*) '!!! Please, your input files and restart   !!!'
        STOP
      ENDIF
C     compute model norm
      CALL CompModelNorm(NNT,b1,Repm,MNorm)
C     Compute r'*b1
      DO ii = 1,NNT(3)
        sumx = D0
        DO jj = 1,ii
          sumx = sumx + BB(jj,ii)*b1(jj)
        ENDDO
        x1(ii) = sumx
      ENDDO

      RETURN
      END ! CompModelSubSpace

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompModelFullSpace(lm,NNT,Repm,CDhat,Cd,x1,MNorm)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER NNT(*)
      REAL*8  lm,CDhat(*),Repm(*),x1(*),MNorm,Cd(*)

      INTEGER ii,jj,kk,info
      REAL*8  RR(LLHMX)

C     data space inversion mode

      DO ii = 1,(NNT(3)+1)*NNT(3)/2
        RR(ii) = Repm(ii)
      ENDDO
      DO ii = 1,NNT(3)
        jj = ii
        kk = ii+(jj-1)*jj/2
        RR(kk) = lm*Cd(ii) + Repm(kk)
        x1(ii) = CDhat(ii)
      ENDDO
C     Using Cholesky decomposition to solve 
      CALL DPPSV('U',NNT(3),1,RR,x1,NN0MX,info)
      IF (info.NE.0) THEN
        WRITE(6,*) '!!! ATTENTION, ERROR WHILE COMPUTING MODEL !!!'
        WRITE(6,*) '!!! Please, your input files and restart   !!!'
        STOP
      ENDIF
C     compute roughness
      CALL CompModelNorm(NNT,x1,Repm,MNorm)

      RETURN
      END ! CompModelFullSpace
 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompSubRespAll(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >           Period,StaNod,SSIndx,DatRes,DatInx,SenInx,StsInx,
     >           dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >           ETOL,MAX_PCG_ITER,FWD_COND,
     >           RhoIn,SSOut,FOUT)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER Nza,Nzb,Nz,Ny,MAX_PCG_ITER,FWD_COND
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      REAL*8  Dzb(*),Dz(*),Dy(*),Czb(*),Cz(*),Cy(*),ETOL
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)
      REAL*8  FOUT(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  RhoIn(NZ0MX,NY0MX)
      REAL*8  SSOut(NMODMX,NSTAMX)

      COMPLEX*16 HXI(MMIMX),HXB(MMBMX)
      COMPLEX*16 EXI(MMIMX),EXB(MMBMX)

      INTEGER im,ir,ip,is,ip2,np1,np2,np3,np4
      REAL*8  stall(NMODMX,NPERMX,NSTAMX)
      REAL*8  per,App(NSTAMX),Phs(NSTAMX)
      REAL*8  AA(MMIMX,4),BD(MMBMX)
      REAL*8  stmul(NPERMX),stmed
      COMPLEX*16 Zyx(NSTAMX),Hxs(NSTAMX),Eys(NSTAMX),Tipper(NSTAMX)
      COMPLEX*16 Zxy(NSTAMX),Exs(NSTAMX),Hys(NSTAMX),Hzs(NSTAMX)

      np1 = NMODMX
      np2 = NRESMX
      np3 = NPERMX
      np4 = NSTAMX

      DO im = 1,NMode
        IF (ModTyp(im).EQ.1)
     >     CALL SetupA_TM(Nzb,Ny,Dzb,Dy,Czb,Cy,RhoIn,AA,BD)
        IF ((ModTyp(im).EQ.2).OR.(ModTyp(im).EQ.3))
     >     CALL SetupA_TE(Nza,Nz,Ny,Dz,Dy,Cz,Cy,RhoIn,AA,BD)
    
        DO ip = 1,NPer(im)
          per = Period(im,ip)

          IF (ModTyp(im).EQ.1) THEN
            CALL Fwd2DTM_PCG(per,Nzb,Ny,Dzb,Dy,RhoIn,AA,BD,
     >           im,NSta,StaNod,ETOL,MAX_PCG_ITER,FWD_COND,
     >           HXI,HXB,App,Phs,Zyx,Hxs,Eys)
          ENDIF
          IF ((ModTyp(im).EQ.2).OR.(ModTyp(im).EQ.3)) THEN
            CALL Fwd2DTE_PCG(per,Nza,Nz,Ny,Dz,Dy,RhoIn,AA,BD,
     >           im,NSta,StaNod,ETOL,MAX_PCG_ITER,FWD_COND,
     >           EXI,EXB,App,Phs,Tipper,Zxy,Exs,Hys,Hzs)
          ENDIF 

          IF (FWD_COND.EQ.1) GOTO 1000

          DO is = 1,NSta(im) 
            DO ir = 1,NRes(im)
             IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN 
               FOUT(im,ir,ip,is) = App(is)
               stall(im,ip,is)   = DatRes(im,ir,ip,is) - App(is)
             ENDIF ! 
             IF ((ResTyp(im,ir).EQ.2).OR.(ResTyp(im,ir).EQ.4)) THEN 
               FOUT(im,ir,ip,is) = Phs(is)
             ENDIF ! 
             IF (ResTyp(im,ir).EQ.5) THEN 
               FOUT(im,ir,ip,is) = DREAL(Tipper(is))
             ENDIF ! 
             IF (ResTyp(im,ir).EQ.6) THEN 
               FOUT(im,ir,ip,is) = IMAG(Tipper(is))
             ENDIF ! 
            ENDDO ! ir
          ENDDO ! is
        ENDDO ! ip

C       static shift
        DO ir = 1,Nres(im)
          IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN 
            DO is = 1,NSta(im)
              IF (SSIndx(im,is).EQ.1) THEN
                ip2 = 0
                DO ip = 1,NPer(im)
                  IF (StsInx(im,ip,is).EQ.1) THEN
                    ip2 = ip2 + 1
                    stmul(ip2) = stall(im,ip,is)
                  ENDIF
                ENDDO
                CALL Median(stmul,ip2,stmed)
                SSOut(im,is) = stmed
              ENDIF
            ENDDO
          ENDIF
        ENDDO ! ir
      ENDDO ! im

      DO im = 1,NMode
        DO ir = 1,NRes(im)
          IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN 
            DO ip = 1,NPer(im)
              DO is = 1,NSta(im)
               FOUT(im,ir,ip,is) = FOUT(im,ir,ip,is)+SSOut(im,is)
              ENDDO ! is
            ENDDO ! ip
          ENDIF
        ENDDO ! ir
      ENDDO ! im

1000  CONTINUE

      RETURN
      END ! CompSubRespAll


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompMisfitAll(NMode,NRes,NPer,NSta,ResTyp,
     >           NNT,NN,DatRes,DatErr,DatInx,FF,
     >           RMS,RMSS,RMSP)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER NN(NMODMX,NRESMX,2),ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  FF(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  RMS(*),RMSS(NMODMX,NRESMX,NSTAMX)
      REAL*8  RMSP(NMODMX,NRESMX,NPERMX)
 
      INTEGER im,ir,ip,is,np,ns,imr
      REAL*8  msfit(NMODMX,NRESMX),rserr(NMODMX,NRESMX)
      REAL*8  misfit(NMODMX,NRESMX,NPERMX,NSTAMX),overallmisfit
      REAL*8  mss,msp,rss,rsp,mf

      overallmisfit = D0
      DO im = 1,NMode
       DO ir = 1,NRes(im)
        msfit(im,ir) = D0
        rserr(im,ir) = D0
        DO ip = 1,NPer(im)
         DO is = 1,NSta(im)
          IF (DatInx(im,ir,ip,is).EQ.1) THEN
           mf = DatRes(im,ir,ip,is) - FF(im,ir,ip,is)
           misfit(im,ir,ip,is) = mf/DatErr(im,ir,ip,is)
           msfit(im,ir)  = msfit(im,ir)+misfit(im,ir,ip,is)**D2
          ENDIF
         ENDDO ! is
        ENDDO ! ip
        overallmisfit = msfit(im,ir)+overallmisfit
       ENDDO ! ir
      ENDDO ! im

      RMS(1)    = DSQRT(overallmisfit/NNT(2))
      imr = 1
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          imr = imr+1
          RMS(imr)    = DSQRT(msfit(im,ir)/NN(im,ir,1))
 
          DO is = 1,NSta(im)
            np  = 0
            mss = D0
            rss = D0
            DO ip = 1,NPer(im)
              IF (DatInx(im,ir,ip,is).EQ.1) THEN
                mss = mss + misfit(im,ir,ip,is)**D2
                np  = np+1
              ENDIF
            ENDDO ! ip
            RMSS(im,ir,is)    = DSQRT(mss/np)
          ENDDO ! is
 
          DO ip = 1,NPer(im)
            msp = D0
            rsp = D0
            ns  = 0

            DO is = 1,NSta(im)
              IF (DatInx(im,ir,ip,is).EQ.1) THEN
                msp = msp + misfit(im,ir,ip,is)**D2
                ns  = ns+1
              ENDIF
            ENDDO ! is
            RMSP(im,ir,ip)    = DSQRT(msp/ns)
          ENDDO ! ip
 
        ENDDO ! ir
      ENDDO ! im

      RETURN
      END ! CompMisfitAll()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

      SUBROUTINE CompMisfit(NMode,NRes,NPer,NSta,NNT,
     >           DatRes,DatErr,DatInx,FF,Misfit)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  FF(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  Misfit

      INTEGER im,ir,ip,is
      REAL*8  overallmisfit,mf
 
      overallmisfit = D0
      DO im = 1,NMode
       DO ir = 1,NRes(im)
        DO ip = 1,NPer(im)
         DO is = 1,NSta(im)
          IF (DatInx(im,ir,ip,is).EQ.1) THEN
           mf = DatRes(im,ir,ip,is) - FF(im,ir,ip,is)
           mf = mf/DatErr(im,ir,ip,is)
           overallmisfit = overallmisfit+mf**D2
          ENDIF
         ENDDO ! is
        ENDDO ! ip
       ENDDO ! ir
      ENDDO ! im

      Misfit    = DSQRT(overallmisfit/NNT(2))

      RETURN
      END ! CompMisfit()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C


      SUBROUTINE CompRespMisfit(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >           Period,StaNod,SSIndx,NNT,DatRes,DatErr,DatInx,
     >           SenInx,StsInx,dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,
     >           Czb,Cz,Cy,ETOL,MAX_PCG_ITER,FWD_COND,DoSmooth,
     >           RhoC,SSC,FFC,fc)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER DoSmooth,Nza,Nzb,Nz,Ny,MAX_PCG_ITER,FWD_COND
      INTEGER NMode,NRes(*),NPer(*),NSta(*),NNT(*)
      INTEGER ModTyp(*),ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      REAL*8  Dzb(*),Dz(*),Dy(*),Czb(*),Cz(*),Cy(*),ETOL
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)
      REAL*8  RhoC(NZ0MX,NY0MX),FFC(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  SSC(NMODMX,NSTAMX),fc


      CALL CompSubRespAll(NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     Period,StaNod,SSIndx,DatRes,DatInx,SenInx,StsInx,
     >     dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,
     >     ETOL,MAX_PCG_ITER,FWD_COND,
     >     RhoC,SSC,FFC)
      CALL CompMisfit(NMode,NRes,NPer,NSta,NNT,
     >     DatRes,DatErr,DatInx,FFC,fc)

      RETURN
      END ! CompRespMisfit

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Model Norm in model space =  || (m-mo)^T*Cm^{-1}*(m-m0) || 
C                in data  space =  || b^{T}*R*b || 

      SUBROUTINE CompModelNorm(NNT,b1,Repm,MNorm)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
      INTEGER NNT(*)
      REAL*8  b1(*),Repm(*),MNorm

      INTEGER ii,jj,kk
      REAL*8  roh(LL0MX),sumx

      REAL*8  DDOT

      DO ii = 1,NNT(3)
        sumx = D0
        DO jj = 1,NNT(3)
          IF (jj.GE.ii) THEN
            kk = ii + (jj-1)*jj/2
          ELSE
            kk = jj + (ii-1)*ii/2
          ENDIF
          sumx = sumx + REPM(kk)*b1(jj)
        ENDDO
        roh(ii) = sumx
      ENDDO
      MNorm = DDOT(NNT(3),b1,1,roh,1)

      RETURN
      END ! CompModelNorm()

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
