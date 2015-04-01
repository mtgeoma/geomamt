
C     ----------------------------- REBOCC ----------------------------------
C
C     2D MT AND MV INVERSION : REDUCED BASIS OCCAM (DATA SUBSPACE) METHOD
C     REBOCC inversion Version 1.0; 
C     Released date :  Jan 20, 1999
C
C     Author Weerachai SIRIPUNVARAPORN
C            College of Oceanic and Atmospheric Science
C            Oregon State University
C            Corvallis, OR 97331
C            wsiripun@oce.orst.edu
C     Last update January 15, 1999.
C
C     When use this program, please refer to 
C    
C     Siripunvaraporn, W. and G. Egbert, 1999,
C     REBOCC: An Efficient Data Subspace Inversion for MT data,
C     Geophysics, submitted.
C
C
C     These programs are provided with the understanding that they will not
C     be redistributed to others. 
C     Please forward requests for copies to  
C             wsiripun@oce.orst.edu or egbert@oce.orst.edu
C     Any changes or modifications of the programs are prohibited.
C
C     If you use these programs, plese send me an email so I can put your name
C     on the list for future releases, bug fixed or any modifications.
C     Any questions or comments about the programs, 
C     please direct it to    wsiripun@oce.orst.edu
C
C     Copyright (c) 1999 by Weerachai Siripunvaraporn and Gary Egbert
C     ------------------------------------------------------------------------

      program rebocc

      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

C
C     MAIN VARIABLES
C

      CHARACTER*80 CStart,CData(NMODMX),CSens(NMODMX),CDist(NMODMX)
      CHARACTER*80 CModel,CPriorModel,CModelControl,COut
      CHARACTER*80 CTitle(NMODMX)

      INTEGER NMode,Ny0,Ny,Nzb,Nza,Nz
      REAL*8  Dy0(NY0MX) ,Dy(NY0MX)
      REAL*8  Dz(NZ0MX)  ,Dzb(NZ0MX) ,Dza(NZ0MX)
      REAL*8  Cy(NY1MX)  ,Cz(NZ1MX)  ,Czb(NZ1MX)
      REAL*8  YDis(NY1MX),ZDis(NY1MX)
      REAL*8  IniRho(NZ0MX,NY0MX)
      REAL*8  PriRho(NZ0MX,NY0MX)

      INTEGER NNT(3),NN(NMODMX,NRESMX,2)
      INTEGER ModTyp(NMODMX),NRes(NMODMX),NPer(NMODMX),NSta(NMODMX)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      INTEGER StaNod(NMODMX,NSTAMX)
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  StaPos(NMODMX,NSTAMX),StaLoc(NMODMX,NSTAMX)
      REAL*8  SSPara(NMODMX,NSTAMX)
      REAL*8  ErrFlr(NMODMX,NRESMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  SknDepth(NMODMX,NPERMX)
      REAL*8  Cd(NN0MX)

      INTEGER ModGrd(NZ0MX,NY0MX),DFStatus(MM0MX)

      REAL*8  CmNorm(MM0MX)
      REAL*8  HDiff(2,NZ0MX,NY0MX),VDiff(2,NZ0MX,NY0MX)
      REAL*8  HGamma(NZ0MX,NY0MX),VGamma(NZ0MX)
      REAL*8  CurrRho(NZ0MX,NY0MX)
      REAL*8  PrevRho(NZ0MX,NY0MX)

      REAL*8  dATM(MM0MX,3),dATE(MM0MX,4)

      REAL*8  BB(NN0MX,LL0MX),Tau(NN0MX),Work(NN0MX)

      INTEGER CmHtIndx(NMODMX,NRESMX,NPERMX)
      REAL*8  CmHt(MM0MX,LL0MX)
      REAL*8  FF(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  TOL(ITERMX),ITOL,PrevTol
      REAL*8  CDhat(NN0MX),CDhatNorm,QCDhatNorm,MinM
      REAL*8  DD(LL0MX)
      REAL*8  Repm(LLHMX)

      REAL*8  RMS(NMODMX*NRESMX+1),RMSS(NMODMX,NRESMX,NSTAMX)
      REAL*8  RMSP(NMODMX,NRESMX,NPERMX)
     
      INTEGER Nit,Info,dmode,FlagMu
      INTEGER length,imodel,DoSmooth,Sit
      REAL*8  mtol,XMin,TolMin,MNormMin,MNormPrev,MNormChange
      REAL*8  RhoChange
      CHARACTER*3 num
      CHARACTER*5 mmd
      CHARACTER*80 cout2,cout3,ctext

      INTEGER INITIAL_ITER_NO,LOGFILE_SCREEN,DTIME_STEP,FWD_Only
      INTEGER MAX_SEARCH_LM,MAX_PCG_ITER,CONT_HIGHER_RMS,
     >        CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,FIX_LM
      REAL*8  MIN_MODEL_CHANGE,MIN_MNORM_CHANGE
      REAL*8  ETOL,SD_RMS,REQUIRED_RMS,PARABOLIC_CORRECTION
      REAL*8  DLENGTH_HOR,DLENGTH_VER,BackGround_Rho
      REAL*8  MIN_LM,MAX_LM,STARTING_LM,STEPSIZE_LM,SMOOTH_SZLM
      INTEGER MAX_ITERATION,MAX_SMOOTHING

      REAL*8  StartPos(NMODMX)

      INTEGER jj
      REAL*8  DDOT

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     BEGIN ACCEPTING INPUTS AND PARAMETERS 
C

      CALL SetDefaultValue(SD_RMS,ETOL,MAX_PCG_ITER,
     >     MIN_MODEL_CHANGE,MIN_MNORM_CHANGE,PARABOLIC_CORRECTION,
     >     INITIAL_ITER_NO,LOGFILE_SCREEN,CONT_HIGHER_RMS,
     >     CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,FWD_Only,
     >     DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,
     >     MIN_LM,MAX_LM,STARTING_LM,FIX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >     MAX_SEARCH_LM)
      imodel = 30

      WRITE(6,9000)
      WRITE(6,*) ' PLEASE ENTER THE STARTUP FILE NAME : '
      READ(5,'(a80)')  CStart


      CALL ReadStartFile(CStart,
     >     NMode,CData,CSens,StartPos,CDist,
     >     COut,INITIAL_ITER_NO,LOGFILE_SCREEN,
     >     CModel,CPriorModel,BackGround_Rho,CModelControl,FWD_Only,
     >     REQUIRED_RMS,MAX_ITERATION,MAX_SMOOTHING,
     >     CONT_HIGHER_RMS,SD_RMS,
     >     CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,
     >     MIN_LM,MAX_LM,STARTING_LM,FIX_LM,
     >     STEPSIZE_LM,SMOOTH_SZLM,MAX_SEARCH_LM,ETOL,MAX_PCG_ITER,
     >     MIN_MODEL_CHANGE,MIN_MNORM_CHANGE,
     >     PARABOLIC_CORRECTION)


      CALL Lenb(COut,length)
      cout3 = COut(1:length)//'.log'
      OPEN(UNIT=99,FILE=cout3,STATUS='unknown')
      IF (LOGFILE_SCREEN.EQ.1) THEN 
        WRITE(6,9000)
        WRITE(6,9998)
        WRITE(6,9000)
      ENDIF
      WRITE(99,9000)
      WRITE(99,9998)
      WRITE(99,9000)
 


C     Read Starting Model
      CALL ReadModel(CModel,
     >     Ny0,Nzb,Nza,Nz,Dy0,Dzb,Dza,Dz,IniRho,YDis,ZDis,Cy,Czb,Cz,
     >     BackGround_Rho)


C     Read Prior Model 
      IF (CPriorModel.EQ.'default') THEN
        CALL ConstantMatrixR8(PriRho,NZ0MX,NY0MX,Nzb,Ny0,D1)
      ELSE
        CALL ReadPriorModel(CPriorModel,Ny0,Nzb,Nza,Dy0,Dzb,Dza,
     >                      PriRho)


        CALL Lenb(COut,length)
        cout2 = COut(1:length)//'.pri'
        OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')
        CALL WriteModel(imodel,-1,D0,D0,
     >     NMode,Nza,Nzb,Ny0,Dza,Dzb,Dy0,PriRho)
        CLOSE(imodel)
      ENDIF

C     Read Model Control 
      CALL ReadModelControl(CModelControl,Ny0,Nzb,
     >   DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,ModGrd,DFStatus)
      DO jj = 1,Ny0*Nzb
        IF (DFStatus(jj).EQ.0) THEN
          IF (CPriorModel.EQ.'default') THEN
            WRITE(6,7010)
            WRITE(6,*) '!!! Program need Prior Model file !!!'
            WRITE(6,7050)
            STOP
          ENDIF
          GOTO 100
        ENDIF
      ENDDO
100   CONTINUE



      CALL ReadData(LOGFILE_SCREEN,CData,NMode,
     >     CTitle,ModTyp,NRes,ResTyp,NPer,NSta,Period,
     >     StaLoc,ErrFlr,DatRes,DatErr,DatInx,Cd,COut)


      CALL ReadSenInc(CSens,NMode,NRes,NPer,NSta,DatInx,
     >     SenInx,NNT,NN,COut,ModTyp)
      CALL SetSknDepth(BackGround_Rho,NMode,NPer,Period,SknDepth)


      CALL ReadStaticInc(CDist,NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >     DatInx,SSIndx,SSPara,StsInx)


      IF (NNT(3).GT.LL0MX) THEN
        WRITE(6,7010)
        WRITE(6,*) '!!! LL0MX is smaller than actual need of ',NNT(3)
        WRITE(6,7100)
        STOP
      ENDIF
      IF (NNT(2).GT.NN0MX) THEN
        WRITE(6,7010)
        WRITE(6,*) '!!! NN0MX is smaller than actual need of ',NNT(2)
        WRITE(6,7100)
        STOP
      ENDIF

      CALL LocateModelPosition(NMode,NPer,NSta,Period,StartPos,StaLoc,
     >     Nzb,Nz,Ny0,Dy0, Ny,Dy,Cy,YDis,IniRho,PriRho,StaPos,StaNod,
     >     ModGrd,DFStatus)
      CALL CopyMatrixR8(1,Nzb,1,Ny,NZ0MX,NY0MX,IniRho,
     >     1,Nzb,1,Ny,NZ0MX,NY0MX,CurrRho)

      CALL WriteResponse(-1,COut,'a',
     >     NMode,NRes,NPer,NSta,DatRes,ModTyp,ResTyp,Period,StaPos,D1)

      CALL WriteErrorResponse(-2,COut,NMode,NRes,NPer,NSta,DatErr,
     >     ModTyp,ResTyp,Period,StaPos,ErrFlr)

      IF (Ny.GT.NY0MX) THEN
        WRITE(6,7010)
        WRITE(6,*) '!!! NY0MX is smaller than actual need of ',Ny
        WRITE(6,7100)
        STOP
      ENDIF
      IF (Nz.GT.NZ0MX) THEN
        WRITE(6,7010)
        WRITE(6,*) '!!! NZ0MX is smaller than actual need of ',Nz
        WRITE(6,7100)
        STOP
      ENDIF




      IF (LOGFILE_SCREEN.EQ.1) THEN
         CALL WriteIntro(6,NMode,ModTyp,NRes,NPer,NSta,
     >        ResTyp,NNT,NN,Ny,Nzb,
     >        REQUIRED_RMS,SD_RMS,MAX_ITERATION,MAX_SMOOTHING,
     >        MIN_LM,STARTING_LM,MAX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >        CONT_HIGHER_RMS,CONT_NOTFOUND_RMS,
     >        CONT_HIGHER_MNORM,FIX_LM,FWD_Only,BackGround_Rho,
     >        DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,
     >        CPriorModel,DFStatus)
      ENDIF
      CALL WriteIntro(99,NMode,ModTyp,NRes,NPer,NSta,
     >        ResTyp,NNT,NN,Ny,Nzb,
     >        REQUIRED_RMS,SD_RMS,MAX_ITERATION,MAX_SMOOTHING,
     >        MIN_LM,STARTING_LM,MAX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >        CONT_HIGHER_RMS,CONT_NOTFOUND_RMS,
     >        CONT_HIGHER_MNORM,FIX_LM,FWD_Only,BackGround_Rho,
     >        DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,
     >        CPriorModel,DFStatus)

      CALL SetSmooth(NMode,NSta,Nzb,Ny,ZDis,YDis,StaNod,StaPos,
     >   DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,HGamma,VGamma)

      CALL Setup1DCM(Nzb,Ny,Dzb,Dy,Czb,Cy,ModGrd,HGamma,VGamma,
     >   HDiff,VDiff)

      CALL SetupNormCM(Nzb,Ny,Dzb,Dy,HDiff,VDiff,DTIME_STEP,
     >   CmNorm)

      CALL SetupDiffOper(NMode,ModTyp,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Cz,Cy,
     >   dATM,dATE)

C     compute interpolation matrix
      IF (NNT(3).LT.NNT(2)) THEN
        IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8005)  
        CALL CompInterPol(NMode,NRes,NPer,Nsta,NNT,DatInx,SenInx,
     >     Period,SknDepth,StaPos,DatErr,BB,Tau,Work)
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     Start iteration 
C


      DoSmooth = 0
      Sit = 0
      Nit = 0
1000  CONTINUE
      Nit = Nit + 1
      IF (DoSmooth.EQ.1) Sit = Sit + 1



      IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8007)  
C     Form Sensitivity Matrix 
      CALL  SubSens2d(Nit,NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >      Period,StaNod,NNT,DFStatus,SSPara,DatInx,SenInx,
     >      dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,CurrRho,
     >      ETOL,MAX_PCG_ITER,FWD_Only,
     >      CmHtIndx,CmHt,FF)


      CALL CompMisfitAll(NMode,NRes,NPer,NSta,ResTyp,NNT,NN,DatRes,
     >     DatErr,DatInx,FF,RMS,RMSS,RMSP)

      mtol = RMS(1)

      IF (Nit.EQ.1) THEN
        ITOL = mtol
        PrevTol = ITOL
        MNormPrev = D0

        WRITE(num,'(i3)') INITIAL_ITER_NO
        IF (INITIAL_ITER_NO.LE.9) num(1:2) = '00'
        IF ((INITIAL_ITER_NO.GE.10).AND.(INITIAL_ITER_NO.LE.99)) 
     >     num(1:1) = '0'
        CALL Lenb(COut,length)
        cout2 = COut(1:length)//'_model.'//num
        OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')
        CALL WriteModel(imodel,Nit-1,mtol,D0,
     >       NMode,Nza,Nzb,Ny,Dza,Dzb,Dy,CurrRho)

        CALL WriteResponse(Nit-1,COut,num,
     >     NMode,NRes,NPer,NSta,FF,ModTyp,ResTyp,Period,StaPos,RMS)
        CALL WriteRMS(Nit-1,COut,num, 
     >     NMode,NRes,NPer,NSta,ModTyp,ResTyp,RMS,RMSS,RMSP)
        CALL WriteDistortion(Nit-1,COut,num,
     >     NMode,NRes,NPer,NSta,ModTyp,ResTyp,SSIndx,SSPara,StsInx) 



        IF (FWD_Only.EQ.1) THEN
          IF (LOGFILE_SCREEN.EQ.1) THEN
            WRITE(6,9000)
            WRITE(6,9040) ITOL
            WRITE(6,9000)
          ENDIF
          WRITE(99,9000)
          WRITE(99,9040) ITOL
          WRITE(99,9000)
          STOP
        ENDIF
      ELSE
        TOL(Nit-1) = mtol
        PrevTol    = mtol
        MNormPrev  = MNormMin
      ENDIF

      CALL CopyMatrixR8(1,Nzb,1,Ny,NZ0MX,NY0MX,CurrRho,
     >                  1,Nzb,1,Ny,NZ0MX,NY0MX,PrevRho)

      CALL CompDHAT(Nzb,Ny,NMode,NRes,NPer,NSta,NNT,
     >     Period,SknDepth,StaPos,DatRes,DatErr,FF,DatInx,SenInx,
     >     CurrRho,PriRho,DFStatus,CmHtIndx,CmHt,
     >     CDhat)

      CALL CompSubRepm(LOGFILE_SCREEN,
     >     Nzb,Ny,Dzb,Dy,NMode,NRes,NPer,NSta,Period,NNT,
     >     DatInx,SenInx,DFStatus,DTIME_STEP,HDiff,VDiff,
     >     CmNorm,CmHtIndx,CmHt,DD,Repm)

      IF (LOGFILE_SCREEN.EQ.1) WRITE(6,8020)  
      IF (NNT(3).LT.NNT(2)) THEN
        CDhatNorm = DDOT(NNT(2),CDhat,1,CDhat,1)
        dmode = 0 
        CALL MulRepm(NNT,BB,Repm)
        CALL DORMQR('L','T',NNT(2),1,NNT(3),BB,NN0MX,Tau,
     >       CDhat,NN0MX,Work,NN0MX,Info)
        QCDhatNorm = DDOT(NNT(3),CDhat,1,CDhat,1)
        MinM = DSQRT((CDhatNorm-QCDhatNorm)/NNT(2))
      ENDIF
      IF (NNT(3).EQ.NNT(2)) THEN
        dmode = 1
      ENDIF

      IF (Nit.EQ.1) THEN
        IF (LOGFILE_SCREEN.EQ.1) THEN
          WRITE(6,9000)
          WRITE(6,9040) ITOL
          WRITE(6,9100)  MinM  
          WRITE(6,9000)
        ENDIF
        WRITE(99,9000)
        WRITE(99,9040) ITOL
        WRITE(99,9100) MinM
        WRITE(99,9000)
      ELSE
        IF (LOGFILE_SCREEN.EQ.1) THEN
           WRITE(6,9000)
           WRITE(6,9100)  MinM  
           WRITE(6,9000)
        ENDIF
        WRITE(99,9000)
        WRITE(99,9100) MinM
        WRITE(99,9000)
      ENDIF

      CALL SearchLM(LOGFILE_SCREEN,
     >     REQUIRED_RMS,SD_RMS,CONT_HIGHER_RMS,
     >     MIN_LM,MAX_LM,STARTING_LM,FIX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >     MAX_SEARCH_LM,
     >     PARABOLIC_CORRECTION,DoSmooth,dmode,
     >     NMode,NRes,NPer,NSta,ModTyp,ResTyp,Period,DatRes,DatErr,Cd,
     >     NN,NNT,DatInx,SenInx,StsInx,StaNod,SSIndx,
     >     DFStatus,PriRho,BB,Repm,CDhat,DD,CmHtIndx,CmHt,
     >     dATM,dATE,Nza,Nzb,Nz,Ny,Dzb,Dz,Dy,Czb,Cz,Cy,PrevTol,
     >     ETOL,MAX_PCG_ITER,FlagMu,
     >     XMin,TolMin,MNormMin,CurrRho,FF,SSPara,
     >     RMS,RMSS,RMSP)


      WRITE(num,'(i3)') Nit+INITIAL_ITER_NO
      IF (Nit+INITIAL_ITER_NO.LE.9) num(1:2) = '00'
      IF ((Nit+INITIAL_ITER_NO.GE.10).AND.
     >    (Nit+INITIAL_ITER_NO.LE.99)) num(1:1) = '0'

      CALL Lenb(COut,length)
      cout2 = COut(1:length)//'_model.'//num
      OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')
      CALL WriteModel(imodel,Nit,TolMin,MNormMin,NMode,
     >     Nza,Nzb,Ny,Dza,Dzb,Dy,CurrRho)

      CALL WriteResponse(Nit,COut,num,
     >     NMode,NRes,NPer,NSta,FF,ModTyp,ResTyp,Period,StaPos,RMS)
      CALL WriteRMS(Nit,COut,num, 
     >     NMode,NRes,NPer,NSta,ModTyp,ResTyp,RMS,RMSS,RMSP)
      CALL WriteDistortion(Nit,COut,num,
     >     NMode,NRes,NPer,NSta,ModTyp,ResTyp,SSIndx,SSPara,StsInx) 

      CALL ModelChange(Nzb,Ny,PrevRho,CurrRho,RhoChange,
     >     MNormPrev,MNormMin,MNormChange)

      IF (Nit.GT.0) THEN
        IF (LOGFILE_SCREEN.EQ.1) THEN
          WRITE(6,9000)
          WRITE(6,9010) Nit,TolMin
          IF (DoSmooth.EQ.1) WRITE(6,9020) Sit,MNormMin
          IF (DoSmooth.EQ.0) WRITE(6,9035) MNormMin
          WRITE(6,9200) DABS(RhoChange)
          IF (DoSmooth.EQ.1) WRITE(6,9210) MNormChange
          WRITE(6,9000)
        ENDIF
        WRITE(99,9000)
        WRITE(99,9010) Nit,TolMin
        IF (DoSmooth.EQ.1) WRITE(99,9020) Sit,MNormMin
        IF (DoSmooth.EQ.0) WRITE(99,9035) MNormMin
        WRITE(99,9200) DABS(RhoChange)
        IF (DoSmooth.EQ.1) WRITE(99,9210) MNormChange
        WRITE(99,9000)
      ENDIF

      IF (DoSmooth.EQ.0) THEN
        IF (Nit.GE.MAX_ITERATION) GOTO 6000
        IF (FlagMu.EQ.2)  GOTO 1000
        IF (FlagMu.EQ.3)  THEN
           IF (CONT_HIGHER_RMS.EQ.0) THEN
             ctext = 
     >       '!!! Higher RMS than previous iteration is found !!!' 
             IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,1)
             CALL Diaglog(99,ctext,1)
             STOP
           ELSE
             GOTO 1000
           ENDIF
        ENDIF
        IF (DABS(RhoChange).LE.MIN_MODEL_CHANGE) THEN
          WRITE(mmd,'(f5.2)') MIN_MODEL_CHANGE
          ctext = 
     >    '!!! Relative model change is less than '//mmd//' % !!!'
          IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,1)
          CALL Diaglog(99,ctext,1)
          STOP
        ENDIF
      ELSE
        IF (Sit.NE.0) THEN
          IF (Sit.GE.MAX_SMOOTHING) GOTO 6000
          IF (FlagMu.NE.1)  THEN
            IF (CONT_NOTFOUND_RMS.EQ.1) THEN
              ctext = 
     > '!!! Desired RMS can not be achieved, go for next iteration !!!'
              IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,0)
              CALL Diaglog(99,ctext,0)
              GOTO 1000
            ELSE
              ctext = 
     > '!!! Desired RMS can not be achieved, choose to quit'
              IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,1)
              CALL Diaglog(99,ctext,1)
              STOP
            ENDIF
          ELSE 
            IF (FlagMu.EQ.1) THEN
              IF ((MNormChange.LE.MIN_MNORM_CHANGE).AND.
     >            (MNormChange.GT.D0)) THEN 
                 WRITE(mmd,'(f5.2)') MIN_MNORM_CHANGE
                 ctext = 
     >   '!!! Relative model norm change is less than '//mmd//' % !!!'
                 IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,1)
                 CALL Diaglog(99,ctext,1)
                 STOP
              ELSE
                IF (MNormChange.LT.D0) THEN
                  IF (CONT_HIGHER_MNORM.EQ.0) THEN
                   ctext = 
     > '!!! Higher Model Norm is found, choose to quit !!!'
                   IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,1)
                   CALL Diaglog(99,ctext,1)
                   STOP
                  ELSE
                    ctext = 
     > '!!! Higher Model Norm is found, go for next iteration !!!'
                    IF (LOGFILE_SCREEN.EQ.1) CALL Diaglog(6,ctext,0)
                    CALL Diaglog(99,ctext,0)
                    GOTO 1000
                  ENDIF
                ENDIF
                GOTO 1000
              ENDIF
            ENDIF
          ENDIF 
        ELSE
          GOTO 5000
        ENDIF
      ENDIF

5000  CONTINUE
      IF (FlagMu.EQ.1) THEN
        GOTO 1000
      ENDIF
6000  CONTINUE
      IF (LOGFILE_SCREEN.EQ.1) THEN
        WRITE(6,9010) Nit,TolMin
        WRITE(6,9020) Sit,MNormMin
      ENDIF
      WRITE(99,9010) Nit,TolMin
      WRITE(99,9020) Sit,MNormMin
      CLOSE(99)

7000  CONTINUE
7010  FORMAT('!!! ATTENTION, ERROR IN INPUT FILE !!!')
7050  FORMAT('!!! Please, have file ready and restart !!!')
7100  FORMAT('!!! Please, correct and restart !!!')

8005  FORMAT('*** forming interpolation matrix ***')
8007  FORMAT('*** forming sensitivity matrix ***')
8020  FORMAT('*** transforming cross-product matrix ***')

 
9000  FORMAT('------------------------------------------------------')
9010  FORMAT('ITERATION # ',i5,' RMS                  = ',f12.4)
9020  FORMAT('SMOOTHING # ',i5,' MODEL NORM           = ',f12.4)
9030  FORMAT('ITERATION # ',i5,' RMS                  = ',f12.4)
9035  FORMAT('                  MODEL NORM           = ',f12.4)

9040  FORMAT('STARTING RMS                           = ',f12.4)

9100  FORMAT('APPROX MIN ACHIEVABLE RMS OF NEXT ITER = ',f12.4)

9200  FORMAT('        RELATIVE CHANGE OF MODEL (%)   = ',f12.4)
9210  FORMAT('     RELATIVE CHANGE OR MODEL NORM (%) = ',f12.4)

9998  FORMAT('----- REBOCC INVERSION FOR MT DATA : VERSION 1.0 -----')

C     End iteration
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      END
C
C     END OF 2D MT AND MV INVERSION PROGRAM
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE Diaglog(iu,ctext,istop)
      INTEGER iu,istop
      CHARACTER*80 ctext

      WRITE(iu,9000)
      WRITE(iu,*) ctext
      IF (istop.EQ.1) THEN
        WRITE(iu,9810)
      ENDIF
      WRITE(iu,9000)
  
9000  FORMAT('------------------------------------------------------')
9810  FORMAT('!!!!!!!!!! REBOCC inversion is finished !!!!!!!!!')

      END ! Diaglog

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
