CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C     SET CONSTANT PARAMETERS                                        C
C                                                                    C
      SUBROUTINE SetDefaultValue(
     >    SD_RMS,ETOL,MAX_PCG_ITER,
     >    MIN_MODEL_CHANGE,MIN_MNORM_CHANGE,PARABOLIC_CORRECTION,
     >    INITIAL_ITERNO,LOGFILE_SCREEN,CONT_HIGHER_RMS,
     >    CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,FWD_Only,
     >    DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,
     >    MIN_LM,MAX_LM,STARTING_LM,FIX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >    MAX_SEARCH_LM)
   
      INCLUDE "constant.h"

      INTEGER MAX_PCG_ITER,MAX_SEARCH_LM,INITIAL_ITERNO,FWD_Only,
     >        LOGFILE_SCREEN,CONT_HIGHER_RMS,CONT_NOTFOUND_RMS,
     >        DTIME_STEP,CONT_HIGHER_MNORM,FIX_LM
      REAL*8  SD_RMS,ETOL,MIN_MODEL_CHANGE,MIN_MNORM_CHANGE,
     >        DLENGTH_HOR,DLENGTH_VER,PARABOLIC_CORRECTION,
     >        MIN_LM,MAX_LM,STARTING_LM,STEPSIZE_LM,SMOOTH_SZLM

      REAL*8  pi4

      D0   =   0.0
      D1   =   1.0
      D2   =   2.0
      D3   =   3.0
      D4   =   4.0
      D5   =   5.0
      D6   =   6.0
      D7   =   7.0
      D8   =   8.0
      D9   =   9.0
      D10  =  10.0
      D90  =  90.0
      D180 = 180.0
      D100 = 100.0

      pi4      = 4.0E-07
      PI       = D2*ASIN(D1)
      CondAir  = 1.E-10
      Mue      = pi4*PI

C     Max no of iter. for searching for lagrange multipler
      MIN_LM        = -D5
      MAX_LM        = D10
      STARTING_LM   = D2
      STEPSIZE_LM   = D1/D2
      SMOOTH_SZLM   = D1/D10
C     default is to use LM to search for min misfit in each iteration
      FIX_LM        = 0 
      MAX_SEARCH_LM = 10

C     standard variation of RMS requires
      SD_RMS        = 0.05

C     PCG : forward modeling
      ETOL          = 1.E-08
      MAX_PCG_ITER  = 500

C     Inversion Progress
      MIN_MODEL_CHANGE     = D1
      MIN_MNORM_CHANGE = D1

      INITIAL_ITERNO   = 0
      LOGFILE_SCREEN  = 1
C     inversion will continue if higher rms found on next iteration
C               will stop if set to 0
      CONT_HIGHER_RMS = 1

C     will continue for next iteration, if smoothing not found RMS
      CONT_NOTFOUND_RMS     = 0
      CONT_HIGHER_MNORM     = 0

C     Forward Modeling only (No inversion), default is NO or 0.
      FWD_Only = 0

C     Diffusion Model Covariance
      DTIME_STEP = 1
      DLENGTH_HOR  = D0
      DLENGTH_VER  = D0

C     if parabolic interpolation help found minimum misfit
      PARABOLIC_CORRECTION = D5


      RETURN
      END ! SUBROUTINE SetConst()
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ReadStartFile(CStart,
     >        NMode,CData,CSens,StartPos,CDist,
     >        COut,INITIAL_ITERNO,LOGFILE_SCREEN,
     >        CModel,CPriorModel,BackGround_Rho,CModelControl,FWD_Only,
     >        REQUIRED_RMS,MAX_ITERATION,MAX_SMOOTHING,
     >        CONT_HIGHER_RMS,SD_RMS,
     >        CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,
     >        MIN_LM,MAX_LM,STARTING_LM,FIX_LM,
     >        STEPSIZE_LM,SMOOTH_SZLM,
     >        MAX_SEARCH_LM,ETOL,MAX_PCG_ITER,
     >        MIN_MODEL_CHANGE,MIN_MNORM_CHANGE,
     >        PARABOLIC_CORRECTION)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER NMode,INITIAL_ITERNO,LOGFILE_SCREEN,
     >        MAX_ITERATION,MAX_SMOOTHING,CONT_HIGHER_RMS,
     >        MAX_SEARCH_LM,MAX_PCG_ITER,FWD_Only,
     >        CONT_HIGHER_MNORM,CONT_NOTFOUND_RMS,FIX_LM
      REAL*8  StartPos(NMODMX),BackGround_Rho,
     >        REQUIRED_RMS,SD_RMS,MIN_LM,MAX_LM,STARTING_LM,
     >        STEPSIZE_LM,SMOOTH_SZLM,
     >        ETOL,MIN_MODEL_CHANGE,MIN_MNORM_CHANGE,
     >        PARABOLIC_CORRECTION
      CHARACTER*80 CStart,CData(NMODMX),CSens(NMODMX),CDist(NMODMX),
     >             CModel,CPriorModel,CModelControl,COut

      INTEGER get_no_mode,get_data_file,get_sens_file,get_starting_pos,
     >        get_starting_model,get_required_rms,
     >        get_max_iteration,get_smoothing,get_output_file
      INTEGER get_distort_file 

      INTEGER im,bwrd,ewrd,chkerr,ic
      CHARACTER*80 ctmp,ccom
      CHARACTER*3  answer

      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd

      chkerr = 0

C     required parameters
      get_no_mode        = 0
      get_data_file      = 0
      get_sens_file      = 0
      get_starting_pos   = 0
      get_starting_model = 0
      get_required_rms   = 0
      get_max_iteration  = 0
      get_smoothing      = 0
      get_output_file    = 0
C     not required parameters
      get_distort_file   = 0

      CModelControl = 'default'
      CPriorModel   = 'default'
      BackGround_Rho = D0

      OPEN(UNIT=10,FILE=CStart,STATUS='OLD',ERR=9000)



      im = 0
      ic = 0
100   CONTINUE
      READ(10,'(a80)',END=6000) ctmp
      ic = ic + 1


c     blank line
      IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
         ic = 0
         GOTO 100
      ENDIF
c     comment line
      IF (FindStr(ctmp(1:10),'#').GT.0) THEN
         ic = 0
         GOTO 100
      ENDIF

      IF (FindStr(ctmp,'number_of_mode').GT.0) THEN
        IF (ic.NE.1) THEN
          WRITE(6,8000)
          WRITE(6,8100) 'NUMBER_OF_MODE'
          WRITE(6,8010)
          STOP
        ENDIF
        get_no_mode = 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=200) NMode
        IF (NMode.LE.0) CALL IncorrectInput('NUMBER_OF_MODE      ')
        IF (NMode.GT.NMODMX) THEN
          WRITE(6,8000)
          WRITE(6,8110) NMode,NMODMX
          WRITE(6,8010)
          STOP
        ENDIF
        IF (NMode.GT.3) THEN
          WRITE(6,8000) 
          WRITE(6,8120)
          WRITE(6,8010)
          STOP
        ENDIF
        GOTO 100
200     CALL IncorrectInput('NUMBER_OF_MODE      ')
      ENDIF

      IF (FindStr(ctmp,'data_file').GT.0) THEN
        IF (get_no_mode.EQ.0) THEN
          WRITE(6,8000) 
          WRITE(6,8130) 'NUMBER_OF_MODE','DATA_FILE'
          WRITE(6,8010)
          STOP
        ELSE
          im = im+1
          IF (im.GT.NMode) THEN
            WRITE(6,8000) 
            WRITE(6,8150) 'DATA_FILE','NUMBER_OF_MODE'
            WRITE(6,8010)
            STOP
          ENDIF
          get_data_file = get_data_file + 1
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          CData(im) = ctmp(bwrd:ewrd)
          OPEN(UNIT=11,FILE=CData(im),STATUS='OLD',ERR=300)
          CLOSE(11)
          CDist(im) = 'default'
        ENDIF
        GOTO 100
300     CONTINUE
        WRITE(6,8000) 
        WRITE(6,8300) CData(im)
        WRITE(6,8010)
        STOP
      ENDIF

      IF (FindStr(ctmp,'sens_inclusion').GT.0) THEN
        IF ((get_data_file.EQ.0).OR.
     >      (get_data_file-get_sens_file.EQ.0)) THEN
          WRITE(6,8000)
          WRITE(6,8130) 'DATA_FILE','SENS_INCLUSION'
          WRITE(6,8010)
          STOP
        ENDIF
        IF (get_data_file-get_sens_file.GT.1) THEN
          WRITE(6,8000)
          WRITE(6,8140) get_sens_file+1,'SENS_INCLUSION',get_data_file
          WRITE(6,8010)
          STOP
        ENDIF
        IF (get_data_file-get_sens_file.EQ.1) THEN
          get_sens_file = get_sens_file + 1
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          ccom = ctmp(bwrd:ewrd)
          CSens(im) = ctmp(bwrd:ewrd)
          IF ((ccom(1:6).EQ.'stripe').OR.(ccom(1:7).EQ.'checker')) THEN
          ELSE
            OPEN(UNIT=11,FILE=CSens(im),STATUS='OLD',ERR=400)
            CLOSE(11)
          ENDIF
        ENDIF
        GOTO 100
400     CONTINUE
        WRITE(6,8000) 
        WRITE(6,8300) CSens(im)
        WRITE(6,8010)
        STOP
      ENDIF

      IF (FindStr(ctmp,'left_offset').GT.0) THEN
        IF ((get_data_file.EQ.0).OR.
     >      (get_data_file-get_starting_pos.EQ.0)) THEN
          WRITE(6,8000)
          WRITE(6,8130) 'DATA_FILE','LEFT_OFFSET'
          WRITE(6,8010)
          STOP
        ENDIF
        IF (get_data_file-get_starting_pos.GT.1) THEN
          WRITE(6,8000)
          WRITE(6,8140) get_starting_pos+1,'LEFT_OFFSET',
     >                  get_data_file
          WRITE(6,8010)
          STOP
        ENDIF
        IF (get_data_file-get_starting_pos.EQ.1) THEN
          get_starting_pos = get_starting_pos + 1
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          READ(ctmp(bwrd:ewrd),*,ERR=500) StartPos(im)
          GOTO 100
500       CALL IncorrectInput('LEFT_OFFSET         ')
        ENDIF
      ENDIF


      IF (FindStr(ctmp,'distort_file').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        get_distort_file = get_distort_file + 1
        IF (get_distort_file.GT.get_data_file) THEN
          WRITE(6,8000)
          WRITE(6,8130) 'DATA_FILE','DISTORT_FILE'
          WRITE(6,8010)
          STOP
        ENDIF
        IF (FindStr(ccom,'default').GT.0) THEN
          CDist(im) = 'default'
          GOTO 100
        ENDIF
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) THEN
          CDist(im) = 'default'
          GOTO 100
        ENDIF
        CDist(im) = ccom
        OPEN(UNIT=11,FILE=CDist(im),STATUS='OLD',ERR=600)
        CLOSE(11)
        GOTO 100
600     CONTINUE
        WRITE(6,8000)
        WRITE(6,8300) CDist(im)
        WRITE(6,8010)
        STOP
      ENDIF

      IF (FindStr(ctmp,'output_file').GT.0) THEN
        get_output_file = 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        COut = ctmp(bwrd:ewrd)
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) THEN
           CALL IncorrectInput('OUTPUT_FILE         ')
        ENDIF
        GOTO 100
      ENDIF

      IF (FindStr(ctmp,'initial_iterno').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        IF (FindStr(ccom,'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=700) INITIAL_ITERNO
        GOTO 100
700     CALL IncorrectInput('initial_iterno      ')
      ENDIF

      IF (FindStr(ctmp,'logfile_screen').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        IF (FindStr(ccom,'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        IF (FindStr(ccom,'yes').GT.0) THEN 
           LOGFILE_SCREEN = 1
           GOTO 100
        ENDIF
        IF (FindStr(ccom,'no').GT.0) THEN
           LOGFILE_SCREEN = 0
           GOTO 100
        ENDIF
800     CALL IncorrectInput('logfile_screen      ')
      ENDIF

      IF (FindStr(ctmp,'starting_model').GT.0) THEN
        get_starting_model = 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        CModel = ctmp(bwrd:ewrd)
        OPEN(UNIT=11,FILE=CModel,STATUS='OLD',ERR=900)
        CLOSE(11)
        GOTO 100
900     CONTINUE
        WRITE(6,8000)
        WRITE(6,8300) CModel
        WRITE(6,8010)
        STOP
      ENDIF

      IF (FindStr(ctmp,'prior_model').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) THEN
          CPriorModel = 'default'
          GOTO 100
        ENDIF
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) THEN
          CPriorModel = 'default'
          GOTO 100
        ENDIF
        CPriorModel = ctmp(bwrd:ewrd)
        OPEN(UNIT=11,FILE=CPriorModel,STATUS='OLD',ERR=1000)
        CLOSE(11)
        GOTO 100
1000    CONTINUE
        WRITE(6,8000)
        WRITE(6,8300) CPriorModel
        WRITE(6,8010)
        STOP
      ENDIF

      IF (FindStr(ctmp,'background_rho').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) THEN
          BackGround_Rho = D0
          GOTO 100
        ENDIF
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) THEN
          BackGround_Rho = D0
          GOTO 100
        ENDIF
        READ(ctmp(bwrd:ewrd),*,ERR=1100) BackGround_Rho
        GOTO 100
1100    CALL IncorrectInput('background_rho      ')
      ENDIF

      IF (FindStr(ctmp,'model_control').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) THEN
          CModelControl = 'default'
          GOTO 100
        ENDIF
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) THEN
          CModelControl = 'default'
          GOTO 100
        ENDIF
        CModelControl = ctmp(bwrd:ewrd)
        OPEN(UNIT=11,FILE=CModelControl,STATUS='OLD',ERR=1200)
        CLOSE(11)
        GOTO 100
1200    CONTINUE
        WRITE(6,8000)
        WRITE(6,8300) CModelControl
        WRITE(6,8010)
        STOP
      ENDIF

      IF (FindStr(ctmp,'fwd_only').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) THEN
          FWD_Only = 0
          GOTO 100
        ENDIF
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) THEN
          FWD_Only = 0
          GOTO 100
        ENDIF
        READ(ctmp(bwrd:ewrd),*,ERR=1300) answer
        IF (FindStr(answer,'no').GT.0)  THEN 
          FWD_Only = 0
          GOTO 100
        ENDIF
        IF (FindStr(answer,'yes').GT.0) THEN
          FWD_Only = 1
          GOTO 100
        ENDIF
1300    CALL IncorrectInput('fwd_only            ')
        GOTO 100
      ENDIF


      IF (FindStr(ctmp,'desired_rms').GT.0) THEN
        get_required_rms = 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=1400) REQUIRED_RMS
        IF (REQUIRED_RMS.LE.0) GOTO 1400
        GOTO 100
1400    CALL IncorrectInput('DESIRED_RMS         ')
      ENDIF

      IF (FindStr(ctmp,'max_iteration').GT.0) THEN
        get_max_iteration = 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=1500) MAX_ITERATION
        IF (MAX_ITERATION.LE.0) GOTO 1500
        GOTO 100
1500    CALL IncorrectInput('MAX_ITERATION       ')
      ENDIF

      IF (FindStr(ctmp,'max_smoothing').GT.0) THEN
        get_smoothing = 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=1600) MAX_SMOOTHING
        IF (MAX_SMOOTHING.LE.0) GOTO 1600
        GOTO 100
1600    CALL IncorrectInput('MAX_SMOOTHING       ')
      ENDIF

      IF (FindStr(ctmp,'cont_high_rms').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=1700) answer
        IF (FindStr(answer,'no').GT.0) THEN 
          CONT_HIGHER_RMS = 0
          GOTO 100
        ENDIF
        IF (FindStr(answer,'yes').GT.0) THEN
          CONT_HIGHER_RMS = 1
          GOTO 100
        ENDIF
1700    CALL IncorrectInput('cont_high_rms       ')
      ENDIF

      IF (FindStr(ctmp,'sd_rms').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=1800) SD_RMS
        IF (SD_RMS.LE.0) GOTO 1800
        GOTO 100
1800    CALL IncorrectInput('sd_rms              ')
      ENDIF

      IF (FindStr(ctmp,'cont_notfound').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=1900) answer
        IF (FindStr(answer,'no').GT.0)  THEN 
           CONT_NOTFOUND_RMS = 0
           GOTO 100
        ENDIF 
        IF (FindStr(answer,'yes').GT.0) THEN
           CONT_NOTFOUND_RMS = 1
           GOTO 100
        ENDIF 
1900    CALL IncorrectInput('cont_notfound       ')
      ENDIF

      IF (FindStr(ctmp,'cont_highnorm').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2000) answer
        IF (FindStr(answer,'no').GT.0)  THEN
           CONT_HIGHER_MNORM = 0
           GOTO 100
        ENDIF
        IF (FindStr(answer,'yes').GT.0) THEN
           CONT_HIGHER_MNORM = 1
           GOTO 100
        ENDIF
2000    CALL IncorrectInput('cont_highnorm       ')
      ENDIF

c     IF (FindStr(ctmp,'min_lgm').GT.0) THEN
c       bwrd = BegWrd(ctmp,2)
c       ewrd = EndWrd(ctmp,2)
c       IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
c       IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
c       READ(ctmp(bwrd:ewrd),*,ERR=2100) MIN_LM
c       GOTO 100
c2100    CALL IncorrectInput('min_lgm              ')
c     ENDIF

c     IF (FindStr(ctmp,'max_lgm').GT.0) THEN
c       bwrd = BegWrd(ctmp,2)
c       ewrd = EndWrd(ctmp,2)
c       IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
c       IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
c       READ(ctmp(bwrd:ewrd),*,ERR=2200) MAX_LM
c       GOTO 100
c2200    CALL IncorrectInput('max_lgm              ')
c      ENDIF

      IF (FindStr(ctmp,'starting_lgm').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2300) STARTING_LM
        GOTO 100
2300    CALL IncorrectInput('starting_lgm         ')
      ENDIF

      IF (FindStr(ctmp,'search_lgm').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2400) answer
        IF (FindStr(answer,'no').GT.0)  THEN 
           FIX_LM = 1
           GOTO 100
        ENDIF
        IF (FindStr(answer,'yes').GT.0) THEN
           FIX_LM = 0
           GOTO 100
        ENDIF
2400    CALL IncorrectInput('search_lgm          ')
      ENDIF

      IF (FindStr(ctmp,'stepsize_lgm').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2500) STEPSIZE_LM
        IF (STEPSIZE_LM.LE.0) GOTO 2500
        GOTO 100
2500    CALL IncorrectInput('stepsize_lgm        ')
      ENDIF

      IF (FindStr(ctmp,'smooth_szlgm').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2600) SMOOTH_SZLM
        IF (SMOOTH_SZLM.LE.0) GOTO 2600
        GOTO 100
2600    CALL IncorrectInput('smooth_szlgm         ')
      ENDIF

      IF (FindStr(ctmp,'max_search_lgm').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2700) MAX_SEARCH_LM
        IF (MAX_SEARCH_LM.LE.0) GOTO 2700
        GOTO 100
2700    CALL IncorrectInput('max_search_lgm       ')
      ENDIF

      IF (FindStr(ctmp,'etol').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2800) ETOL
        IF (ETOL.LE.0) GOTO 2800
        GOTO 100
2800    CALL IncorrectInput('etol                ')
      ENDIF

      IF (FindStr(ctmp,'max_pcg_iter').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=2900) MAX_PCG_ITER
        IF (MAX_PCG_ITER.LE.0) GOTO 2900
        GOTO 100
2900    CALL IncorrectInput('max_pcg_iter        ')
      ENDIF

      IF (FindStr(ctmp,'model_change').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=3000) MIN_MODEL_CHANGE
        IF (MIN_MODEL_CHANGE.LE.0) GOTO 3000
        GOTO 100
3000    CALL IncorrectInput('model_change        ')
      ENDIF

      IF (FindStr(ctmp,'mnorm_change').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=3100) MIN_MNORM_CHANGE
        IF (MIN_MNORM_CHANGE.LE.0) GOTO 3100
        GOTO 100
3100    CALL IncorrectInput('mnorm_change        ')
      ENDIF

      IF (FindStr(ctmp,'parabolic_cor').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        IF (FindStr(ctmp(bwrd:ewrd),'default').GT.0) GOTO 100
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80)) GOTO 100
        READ(ctmp(bwrd:ewrd),*,ERR=3200) PARABOLIC_CORRECTION
        IF (PARABOLIC_CORRECTION.LE.0) GOTO 3200
        GOTO 100
3200    CALL IncorrectInput('parabolic_cor       ')
      ENDIF

      WRITE(6,8000)
      WRITE(6,8400) ctmp
      WRITE(6,8010)
      STOP


5000  FORMAT(a80)

6000  CONTINUE

C     required parameters
      IF (get_no_mode.EQ.0) THEN
        WRITE(6,8000)
        WRITE(6,8100) 'NUMBER_OF_MODE'
        chkerr = chkerr + 1
      ENDIF
      IF (get_data_file.NE.NMode) THEN
        WRITE(6,8000)
        WRITE(6,8500) 'DATA_FILE',NMode
        chkerr = chkerr + 1
      ENDIF
      IF (get_sens_file.NE.NMode) THEN
        WRITE(6,8000)
        WRITE(6,8500) 'SENS_INCLUSION',NMode
        chkerr = chkerr + 1
      ENDIF
      IF (get_starting_pos.NE.NMode) THEN
        WRITE(6,8000)
        WRITE(6,8500) 'LEFT_OFFSET',NMode
        chkerr = chkerr + 1
      ENDIF
      IF (get_output_file.EQ.0) THEN
        WRITE(6,8000)
        WRITE(6,8550) 'OUTPUT_FILE'
        chkerr = chkerr + 1
      ENDIF
      IF (get_starting_model.EQ.0) THEN
        WRITE(6,8000)
        WRITE(6,8550) 'STARTING_MODEL'
        chkerr = chkerr + 1
      ENDIF
      IF (get_required_rms.EQ.0) THEN
        WRITE(6,8000)
        WRITE(6,8550) 'REQUIRED_RMS'
        chkerr = chkerr + 1
      ENDIF
      IF (get_max_iteration.EQ.0) THEN
        WRITE(6,8000)
        WRITE(6,8550) 'MAX_ITERATION'
        chkerr = chkerr + 1
      ENDIF
      IF (get_smoothing.EQ.0) THEN
        WRITE(6,8000)
        WRITE(6,8550) 'MAX_SMOOTHING'
        chkerr = chkerr + 1
      ENDIF
      IF (chkerr.GT.0) THEN
        WRITE(6,8010) 
        STOP
      ELSE
        GOTO 7000
      ENDIF


7000  CONTINUE
      CLOSE(10)

      GOTO 9999

8000  FORMAT('!!!!!! ATTENTION, ERROR IN THE STARTUP FILE !!!!!!')
8010  FORMAT('!!!!!!      please, correct and restart     !!!!!!')

8100  FORMAT('!!! Program requires ',a15,' as the first command. ')
8110  FORMAT('!!! Number of Mode = ',i4,' exceeds NMODMX = ',i4,'.')
8120  FORMAT('!!! For REBOCC V1, 3 mode inversion is maximum.')
8130  FORMAT('!!! Please, enter ',a20,' prior to ',a20,' .')
8140  FORMAT('!!! There are only ',i4,a20,' while need ',i4,' .') 
8150  FORMAT('!!! There are more ',a20,' than specify in ',a20)

8300  FORMAT('!!! File name <',a40,'> cannot be found. ')
8400  FORMAT('!!! Keyword ',a15,' is not recognized in startup file.')

8500  FORMAT('!!! There is not enough ',a15,' for ',i1,' mode(s).')
8550  FORMAT('!!! There is no ',a15,' in the startup file.')

9000  CONTINUE
      WRITE(6,8000)
      WRITE(6,8300) CStart
      WRITE(6,8010)
      STOP

9999  CONTINUE

      RETURN
      END ! SUBROUTINE ReadStartFile()
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE IncorrectInput(CString)
      CHARACTER*20 CString

      WRITE(6,8000) 
      WRITE(6,8985) CString
      WRITE(6,8999)
      STOP

8000  FORMAT('!!!!!! ATTENTION, ERROR IN INPUT FILE !!!!!!')
8985  FORMAT('!!!!!!  Incorrecly input at ',a20,' !!!!!!')
8999  FORMAT('!!!!!!  Please, correct and restart   !!!!!!')
 
      END ! INCORRECT_INPUT
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ReadModel(CModel,
     >          Ny0,Nzb,Nza,Nz,Dy0,Dzb,Dza,Dz,Rho,YDis,ZDis,Cy,Czb,Cz,
     >          BGRho)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

C         IO VARIABLES   
      CHARACTER*80 CModel,ctitle
      INTEGER Ny0,Nzb,Nza,Nz
      REAL*8  Dy0(*),Dzb(*),Dza(*),Dz(*),Rho(NZ0MX,NY0MX),
     >        YDis(*),ZDis(*),Cy(*),Czb(*),Cz(*),BGRho
      
C         LOCAL VARIABLES 
      CHARACTER*80 ctmp,ccom

      INTEGER iy,iz,modelmode
      REAL*8  avg_rho

      INTEGER ic,maxcommand,bwrd,ewrd
      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd

      maxcommand = 3
      ic  = 0
      Ny0 = 0
      Nzb = 0
      Nza = 0

C     by default, Nza = 10
      Nza = 10
      Dza(1)  = 10.
      Dza(2)  = 30.
      Dza(3)  = 100.
      Dza(4)  = 300.
      Dza(5)  = 1000.
      Dza(6)  = 3000.
      Dza(7)  = 10000.
      Dza(8)  = 30000.
      Dza(9)  = 100000.
      Dza(10) = 300000.

      OPEN(UNIT=10,FILE=CModel,STATUS='OLD')

10    CONTINUE
      IF (ic.EQ.maxcommand) THEN
        CLOSE(10)
        GOTO 900
      ENDIF
      READ(10,'(a80)',END=20) ctmp

      GOTO 30
20    IF (ic.LT.maxcommand) THEN
        WRITE(6,1000)
        WRITE(6,*) '!!! End of file reached before completed !!!'
        WRITE(6,1100)
        STOP
      ENDIF

30    CONTINUE

c     blank line
      IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
         GOTO 10
      ENDIF
c     comment line
      IF (FindStr(ctmp(1:10),'#').GT.0) THEN
         GOTO 10
      ENDIF
      ic = ic+1

      IF (FindStr(ctmp,'title').GT.0) THEN
        maxcommand = maxcommand + 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*) ctitle
        GOTO 10
      ENDIF

      IF (FindStr(ctmp,'ny').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=40) Ny0
        GOTO 45
40      CONTINUE
        CALL IncorrectInput('       NY           ')
45      CONTINUE
        IF (Ny0.LT.1)     CALL DimenWrong(0,'NY     ',1)
        IF (Ny0.GT.NY0MX) CALL DimenWrong(1,'NY     ',NY0MX)

        READ(10,*,ERR=50) (Dy0(iy),iy=1,Ny0)
        GOTO 55
50      CONTINUE
        CALL IncorrectInput('       DY           ')
55      CONTINUE
        GOTO 10
      ENDIF

      IF (FindStr(ctmp(1:3),'nzb').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=60) Nzb
        GOTO 65
60      CONTINUE
        CALL IncorrectInput('       NZB          ')
65      CONTINUE
        IF (Nzb.LT.1)     CALL DimenWrong(0,'NZB    ',1)
        IF (Nzb.GT.NZ0MX) CALL DimenWrong(1,'NZB    ',NZ0MX)

        READ(10,*,ERR=70) (Dzb(iz),iz=1,Nzb)
        GOTO 75
70      CONTINUE
        CALL IncorrectInput('       DZB          ')
75      CONTINUE
        GOTO 10
      ENDIF

      IF (FindStr(ctmp(1:3),'nza').GT.0) THEN
        maxcommand = maxcommand + 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=80) Nza
        GOTO 85
80      CONTINUE
        CALL IncorrectInput('       nza          ')
85      CONTINUE
        IF (Nza.LT.1)     CALL DimenWrong(0,'nza    ',1)
        IF (Nza.GT.NZ0MX) CALL DimenWrong(1,'nza    ',NZ0MX)

        READ(10,*,ERR=90) (Dza(iz),iz=1,Nza)
        GOTO 95
90      CONTINUE
        CALL IncorrectInput('       dza          ')
95      CONTINUE
        GOTO 10
      ENDIF

      IF (Nza+Nzb.LT.1)     CALL DimenWrong(0,'NZA+NZB',1)
      IF (Nza+Nzb.GT.NZ0MX) CALL DimenWrong(1,'NZA+NZB',NZ0MX)

      IF (FindStr(ctmp,'resistivity_model').GT.0) THEN
        IF (Ny0.EQ.0) THEN
          WRITE(6,1000)
          WRITE(6,*) 
     >    '!!! Ny must be specified before <RESISTIVITY_MODEL> !!!'
          WRITE(6,1100)
          STOP
        ENDIF
        IF (Nzb.EQ.0) THEN
          WRITE(6,1000)
          WRITE(6,*) 
     >    '!!! NZB must be specified before <RESISTIVITY_MODEL> !!!'
          WRITE(6,1100)
          STOP
        ENDIF

        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)

        modelmode = -1
C       by default, resistivity model is in general form
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80))  modelmode = 0
        IF (FindStr(ccom,'index').GT.0)     modelmode = 1
        IF (FindStr(ccom,'layer').GT.0)     modelmode = 2
        IF (FindStr(ccom,'halfspace').GT.0) modelmode = 3

        IF (modelmode.EQ.-1) THEN
          WRITE(6,1000)
          WRITE(6,*) 
     >    '!!! Keyword is not recognized after RESISTIVITY_MODEL!!!'
          WRITE(6,1100)
          STOP
        ELSE
          CALL ReadRho(10,modelmode,Nzb,Ny0,Rho)
        ENDIF
        GOTO 10
      ENDIF

      WRITE(6,1000)
      WRITE(6,1200) ctmp
      WRITE(6,1100)
      STOP

900   CONTINUE

      Nz = Nza+Nzb
      DO iz = 1,Nza
        Dz(Nza-iz+1) = Dza(iz)
      ENDDO
      DO iz = 1,Nzb
        Dz(Nza+iz)   = Dzb(iz)
      ENDDO

C     COMPUTE CUMULATIVE DISTANCE FROM THE LEFT SIDE AND SURFACE
      CALL CumulativeDistance(Ny0,Dy0,YDis)
      CALL CumulativeDistance(Nzb,Dzb,ZDis)

      CALL DistanceBetweenBlocks(Ny0,Dy0,Cy)
      CALL DistanceBetweenBlocks(Nzb,Dzb,Czb)
      CALL DistanceBetweenBlocks(Nz,Dz,Cz)

C     Estimate Background Resistivity
C     By default, using an average value of the first layer of the 
C     resistivity starting model 

      IF (BGRho.EQ.D0) THEN
        avg_rho = D0
        DO iy = 1,Ny0
          avg_rho =  avg_rho + Rho(1,iy)
        ENDDO 
        BGRho = avg_rho/Ny0
      ENDIF



1000  FORMAT('!!!!!!  ATTENTION, ERROR IN STARTING MODEL !!!!!')
1100  FORMAT('!!!!!!     Please, correct and restart     !!!!!')

1200  FORMAT('!!! Keyword ',a20,' is not recognized !!!')

      RETURN
      END ! SUBROUTINE ReadModel()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
      SUBROUTINE ReadRho(io,icase,nz0,ny0,rh)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER nz0,ny0,io
      REAL*8  rh(NZ0MX,NY0MX)
 
      INTEGER iz,iy,ir,nyt,len,irhomax,irho(NZ0MX,NY0MX),nrho,icase,
     >        bwrd,ewrd,iyt
      REAL*8  rhoval(NZ0MX)
      CHARACTER*256 cindex

      INTEGER  BegWrd,EndWrd
      EXTERNAL BegWrd,EndWrd

c     output model
      IF (icase.EQ.0) THEN
        DO iz = 1,nz0
         READ(io,*,ERR=100) (rh(iz,iy),iy=1,ny0)
        ENDDO
        GOTO 200
100     WRITE(6,1000)
        WRITE(6,1030)
        WRITE(6,1100)
        STOP
200     CONTINUE
      ENDIF

c     index model
      IF (icase.EQ.1) THEN
        DO iz = 1,nz0
          nyt = 0
300       CONTINUE
          READ(io,'(a256)') cindex
          bwrd = BegWrd(cindex,1)
          ewrd = EndWrd(cindex,1)
          len  = ewrd - bwrd + 1
          iyt = 0
          DO iy = bwrd,ewrd
            iyt = iyt + 1
            READ(cindex(iy:iy),'(i1)',ERR=400) irho(iz,iyt+nyt)
          ENDDO
          nyt = len + nyt
          IF (nyt.LT.ny0) GOTO 300
        ENDDO

        irhomax = 0
        DO iz = 1,nz0
         DO iy = 1,ny0
          irhomax = MAX(irho(iz,iy),irhomax)
          IF (irho(iz,iy).LE.0) THEN
           WRITE(6,1000)
           WRITE(6,1040)
           WRITE(6,1100)
           STOP
          ENDIF
         ENDDO ! iy
        ENDDO ! iz
        nrho = irhomax
        READ(io,*,ERR=400) (rhoval(ir),ir=1,nrho)

        CALL AssignValIndx(1,nz0,1,ny0,NZ0MX,NY0MX,irho,rhoval,rh)
        GOTO 450
400     WRITE(6,1000)
        WRITE(6,1030)
        WRITE(6,1100)
        STOP
450     CONTINUE
      ENDIF

c     1-d layer model
      IF (icase.EQ.2) THEN
        READ(io,*,ERR=500) (rhoval(iz),iz=1,nz0)
        DO iz = 1,nz0
         DO iy = 1,ny0
           rh(iz,iy) = rhoval(iz)
         ENDDO ! iy
        ENDDO ! iz
        GOTO 550
500     WRITE(6,1000)
        WRITE(6,1030)
        WRITE(6,1100)
        STOP
550     CONTINUE
      ENDIF

c     half space initial model
      IF (icase.EQ.3) THEN
         READ(io,*,ERR=600) rhoval(1)
         IF (rhoval(1).LT.D0) GOTO 600
         GOTO 610
600      CONTINUE
         WRITE(6,1000)
         WRITE(6,1030)
         WRITE(6,1100)
         STOP
610      CONTINUE
         CALL ConstantMatrixR8(rh,NZ0MX,NY0MX,nz0,ny0,rhoval(1))
      ENDIF

      DO iz = 1,nz0
        DO iy = 1,ny0
          IF (rh(iz,iy).LE.D0) THEN
            WRITE(6,1000)
            WRITE(6,1050)
            WRITE(6,1100)
            STOP
          ENDIF
        ENDDO ! iy
      ENDDO ! iz

1000  FORMAT('!!!!!!  ATTENTION, ERROR IN READING MODEL !!!!!')
1030  FORMAT('!!!      Error reading resistivity value ')
1040  FORMAT('!!! Resistivity indices must be positive integers.')
1050  FORMAT('!!! Resistivity values must be positive values.')
1100  FORMAT('!!!!!!     Please, correct and restart     !!!!!')

      RETURN
      END ! ReadRho

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C


      SUBROUTINE DimenWrong(imore,ctext,ivalue)
      INTEGER imore,ivalue
      CHARACTER*7 ctext

      IF (imore.EQ.0) THEN
        WRITE(6,1000)
        WRITE(6,1010) ctext,ivalue
        WRITE(6,1100)
        STOP
      ENDIF

      IF (imore.EQ.1) THEN
        WRITE(6,1000)
        WRITE(6,1020) ctext,ivalue
        WRITE(6,1100)
        STOP
      ENDIF

1000  FORMAT('!!!!!!  ATTENTION, ERROR IN STARTING MODEL !!!!!')
1010  FORMAT('!!! ',a3,' should have value more than ',i3)
1020  FORMAT('!!! ',a3,' should have value less than ',i3)
1100  FORMAT('!!!!!!     Please, correct and restart     !!!!!')

      RETURN
      END ! DimenWrong

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ReadPriorModel(CPriorModel,
     >           Ny0,Nzb,Nza,Dy0,Dzb,Dza,PRho)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

C         IO VARIABLES   
      CHARACTER*80 CPriorModel
      INTEGER Ny0,Nzb,Nza
      REAL*8  PRho(NZ0MX,NY0MX)
      REAL*8  Dy0(*),Dzb(*),Dza(*)
      
C         LOCAL VARIABLES 
      CHARACTER*80 ctmp,ctitle,ccom
      INTEGER iy,iz,modelmode
      INTEGER My0,Mzb
      REAL*8  Ey0(NY0MX),Ezb(NZ0MX)

      INTEGER ic,maxcommand,bwrd,ewrd
      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd
  
      My0 = 0
      Mzb = 0
      maxcommand = 1
      ic = 0

      OPEN(UNIT=10,FILE=CPriorModel,STATUS='OLD')

10    CONTINUE
      IF (ic.EQ.maxcommand) THEN
        CLOSE(10)
        GOTO 900
      ENDIF
      READ(10,'(a80)',END=20) ctmp

      GOTO 30
20    IF (ic.LT.maxcommand) THEN
        WRITE(6,1000)
        WRITE(6,*) '!!! End of file reached before completed !!!'
        WRITE(6,1100)
        STOP
      ENDIF

30    CONTINUE

c     blank line
      IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
         GOTO 10
      ENDIF
c     comment line
      IF (FindStr(ctmp(1:10),'#').GT.0) THEN
         GOTO 10
      ENDIF
      ic = ic+1

      IF (FindStr(ctmp,'title').GT.0) THEN
        maxcommand = maxcommand + 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*) ctitle
        GOTO 10
      ENDIF

      IF (FindStr(ctmp,'ny').GT.0) THEN
        maxcommand = maxcommand + 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=40) My0
        GOTO 45
40      CONTINUE
        CALL IncorrectInput('       NY           ')
45      CONTINUE
        IF (My0.LT.1)     CALL DimenWrong(0,'NY     ',1)
        IF (My0.GT.NY0MX) CALL DimenWrong(1,'NY     ',NY0MX)
        IF (My0.NE.Ny0) THEN
          WRITE(6,1000)
          WRITE(6,1005) 'NY'
          WRITE(6,1100)
          STOP
        ENDIF

        READ(10,*,ERR=50) (Ey0(iy),iy=1,My0)
        DO iy = 1,My0
          IF (Ey0(iy).NE.Dy0(iy)) THEN
            WRITE(6,1000)
            WRITE(6,1005) 'DY '
            WRITE(6,1100)
            STOP
          ENDIF
        ENDDO
        GOTO 55
50      CONTINUE
        CALL IncorrectInput('       DY           ')
55      CONTINUE
        GOTO 10
      ENDIF

      IF (FindStr(ctmp(1:3),'nzb').GT.0) THEN
        maxcommand = maxcommand + 1
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        READ(ctmp(bwrd:ewrd),*,ERR=60) Mzb
        GOTO 65
60      CONTINUE
        CALL IncorrectInput('       NZB          ')
65      CONTINUE
        IF (Mzb.LT.1)     CALL DimenWrong(0,'NZB    ',1)
        IF (Mzb.GT.NZ0MX) CALL DimenWrong(1,'NZB    ',NZ0MX)
        IF (Mzb.NE.Nzb) THEN
          WRITE(6,1000)
          WRITE(6,1005) 'NZB'
          WRITE(6,1100)
          STOP
        ENDIF 

        READ(10,*,ERR=70) (Ezb(iz),iz=1,Mzb)
        DO iz = 1,Mzb
          IF (Ezb(iz).NE.Dzb(iz)) THEN
            WRITE(6,1000)
            WRITE(6,1005) 'DZB'
            WRITE(6,1100)
            STOP
          ENDIF
        ENDDO
        GOTO 75
70      CONTINUE
        CALL IncorrectInput('       DZB          ')
75      CONTINUE
        GOTO 10
      ENDIF

      IF (FindStr(ctmp,'resistivity_model').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)

        modelmode = -1
C       by default, resistivity model is in general form
        IF ((bwrd.EQ.80).AND.(ewrd.EQ.80))  modelmode = 0
        IF (FindStr(ccom,'index').GT.0)     modelmode = 1
        IF (FindStr(ccom,'layer').GT.0)     modelmode = 2
        IF (FindStr(ccom,'halfspace').GT.0) modelmode = 3

        IF (modelmode.EQ.-1) THEN
          WRITE(6,1000)
          WRITE(6,*) 
     >    '!!! Keyword is not recognized after RESISTIVITY_MODEL!!!'
          WRITE(6,1100)
          STOP
        ELSE
          CALL ReadRho(10,modelmode,Nzb,Ny0,PRho)
        ENDIF

        GOTO 10
      ENDIF

      WRITE(6,1000)
      WRITE(6,1200) ctmp
      WRITE(6,1100)
      STOP

900   CONTINUE
      CLOSE(10)

1000  FORMAT('!!!!!!  ATTENTION, ERROR IN PRIOR MODEL !!!!!')
1005  FORMAT('!!! ',a3,' of Prior Model must be the same as ',
     >       ' that of the Starting Model.')
1100  FORMAT('!!!!!!     Please, correct and restart     !!!!!')

1200  FORMAT('!!! Keyword ',a20,' is not recognized !!!')

      RETURN
      END ! SUBROUTINE ReadPriorModel()
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ReadData(LOGFILE_SCREEN,
     >           CData,NMode,CTitle,ModTyp,NRes,ResTyp,
     >           NPer,NSta,Period,StaLoc,ErrFlr,DatRes,DatErr,DatInx,
     >           Cd,COut)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

C            IO VARIABLES 
      CHARACTER*80 CData(NMODMX),CTitle(NMODMX),COut
      INTEGER NMode,LOGFILE_SCREEN

      INTEGER ModTyp(*),NRes(*),NPer(*),NSta(*)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  Period(NMODMX,NPERMX)
      REAL*8  ErrFlr(NMODMX,NRESMX)
      REAL*8  DatRes(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  DatErr(NMODMX,NRESMX,NPERMX,NSTAMX)
      REAL*8  Cd(*)

C           LOCAL VARIABLES
      INTEGER im,ir,ip,is,idx,bwrd,ewrd,length,
     >        ic,irr,useperiod,
     >        maxcommand,idr1,idr2,ier1,ier2,iir1,iir2,ires1,ires2
      CHARACTER*80 ctmp,ccom,cout2
      CHARACTER*30 ctext
      CHARACTER*256 cindex
      CHARACTER*1   crr
      REAL*8  StaLoc(NMODMX,NSTAMX)
      REAL*8  incfac,freq,minerr,eflr
      REAL*8  frequency(NPERMX)
      INTEGER iunit 

      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd

      maxcommand = 5 + 3*NRESMX

C     Read Data
      DO im = 1,NMode
        idr1 = 0
        idr2 = 0
        ier1 = 0
        ier2 = 0
        iir1 = 0
        iir2 = 0
        ires1 = 0
        ires2 = 0

        NRes(im) = 0
        NPer(im) = 0
        NSta(im) = 0

        OPEN(UNIT=10,FILE=CData(im),STATUS='OLD')


        ic = 0
10      CONTINUE
        IF (ic.EQ.maxcommand) THEN
          CLOSE(10)
          GOTO 500
        ENDIF
        READ(10,5000,END=20) ctmp

        GOTO 30
20      IF (ic.LT.maxcommand) THEN
          WRITE(6,2000)
          WRITE(6,*) '!!! End of file reached before completed !!!'
          WRITE(6,2010)
          STOP
        ENDIF

30      CONTINUE

c       blank line
        IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
           GOTO 10
        ENDIF
c       comment line
        IF (FindStr(ctmp(1:10),'#').GT.0) THEN
           GOTO 10
        ENDIF
        ic = ic+1
       

        IF (FindStr(ctmp,'title').GT.0) THEN
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          READ(ctmp(bwrd:ewrd),*) CTitle(im)
          GOTO 10
        ENDIF

        IF (FindStr(ctmp,'mode_type').GT.0) THEN
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          ccom = ctmp(bwrd:ewrd)

          ModTyp(im) = 0
          IF (FindStr(ccom,'tm').GT.0) THEN
            ModTyp(im) = 1
          ENDIF
          IF (FindStr(ccom,'te').GT.0) THEN
            ModTyp(im) = 2
          ENDIF
          IF (FindStr(ccom,'tp').GT.0) THEN
            ModTyp(im) = 3
          ENDIF

          IF (ModTyp(im).EQ.0) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!!! Mode type (TM/TE/TP) must be selected !!!'
            WRITE(6,2010)
            STOP
          ENDIF
          GOTO 10
        ENDIF

        IF (FindStr(ctmp,'number_of_response').GT.0) THEN
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          READ(ctmp(bwrd:ewrd),*,ERR=50) NRes(im)
          GOTO 60
50        CONTINUE
          CALL IncorrectInput('NUMBER_OR_RESPONSE  ')
60        CONTINUE
          IF (NRes(im).LE.0) THEN
            WRITE(6,2000)
            WRITE(6,*) 
     >      '!!! NUMBER_OR_RESPONSE should be positive value !!!'
            WRITE(6,2010)
            STOP
          ENDIF
          IF (NRes(im).GT.NRESMX) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!! NRESMX is smaller than actual need of ',
     >                 NRes(im),' !!!'
            WRITE(6,2010)
            STOP
          ENDIF
          maxcommand = 5 + 3*NRes(im)
          GOTO 10
        ENDIF

        IF ((FindStr(ctmp,'number_of_frequency').GT.0)
     >     .OR.(FindStr(ctmp,'number_of_period').GT.0)) THEN
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          READ(ctmp(bwrd:ewrd),*,ERR=100) NPer(im)
          IF (FindStr(ctmp,'period').GT.0)    useperiod = 1
          IF (FindStr(ctmp,'frequency').GT.0) useperiod = 0
          GOTO 110
100       CONTINUE
          IF (FindStr(ctmp,'period').GT.0)    useperiod = 1
          IF (FindStr(ctmp,'frequency').GT.0) useperiod = 0
          IF (useperiod.EQ.1) THEN
            CALL IncorrectInput('NUMBER_OF_PERIOD    ')
          ELSE
            CALL IncorrectInput('NUMBER_OF_FREQUENCY ')
          ENDIF
          STOP
110       CONTINUE

          IF (NPer(im).LE.0) THEN
            WRITE(6,2000)
            WRITE(6,*) 
     >   '!!! NUMBER_OR_FREQUENCY/PERIOD should be positive value !!!'
            WRITE(6,2010)
            STOP
          ENDIF
          IF (NPer(im).GT.NPERMX) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!! NPERMX is smaller than actual need of ',
     >                 NPer(im),' !!!'
            WRITE(6,2010)
            STOP
          ENDIF

          READ(10,*,ERR=150) (frequency(ip),ip=1,NPer(im))
          GOTO 160
150       CONTINUE
          CALL IncorrectInput('<period(s)>         ')
160       CONTINUE
          IF (useperiod.EQ.1) THEN
            DO ip = 1,NPer(im)
              Period(im,ip) = frequency(ip)
            ENDDO ! ip
          ELSE
            DO ip = 1,NPer(im)
              Period(im,ip) = D1/frequency(ip)
            ENDDO ! ip
          ENDIF
          GOTO 10
        ENDIF

        IF (FindStr(ctmp,'number_of_station').GT.0) THEN
          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          READ(ctmp(bwrd:ewrd),*,ERR=200) NSta(im)
          GOTO 210
200       CONTINUE
          CALL IncorrectInput('NUMBER_OF_STATION   ')
210       CONTINUE

          IF (NSta(im).LE.0) THEN
            WRITE(6,2000)
            WRITE(6,*) 
     >   '!!! NUMBER_OR_STATION should be positive value !!!'
            WRITE(6,2010)
            STOP
          ENDIF
          IF (NSta(im).GT.NSTAMX) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!! NSTAMX is smaller than actual need of ',
     >                 NSta(im),' !!!'
            WRITE(6,2010)
            STOP
          ENDIF

          READ(10,*,ERR=250) (StaLoc(im,is),is=1,NSta(im))
          GOTO 260
250       CONTINUE
          CALL IncorrectInput('<station(s)        >')
260       CONTINUE
          GOTO 10
        ENDIF

        IF (FindStr(ctmp,'data_response_no_').GT.0) THEN
          IF (NRes(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_RESPONSE must be given ',
     >                  'prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          IF (NPer(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_PERIOD/NUMBER_OF_FREQUENCY',
     >       ' must be given prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          IF (NSta(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_STATION ',
     >       ' must be given prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          bwrd = BegWrd(ctmp,1)
          ewrd = EndWrd(ctmp,1)
          IF (FindStr(ctmp(ewrd:ewrd),'*').GT.0) THEN
            READ(ctmp(ewrd-1:ewrd-1),'(i1)') irr
          ELSE
            READ(ctmp(ewrd:ewrd),'(i1)') irr
          ENDIF
          IF ((irr.LT.1).OR.(irr.GT.NRes(im))) GOTO 400
          IF ((idr1.EQ.1).AND.(irr.EQ.1)) GOTO 400
          IF ((idr2.EQ.1).AND.(irr.EQ.2)) GOTO 400
          IF (irr.EQ.1) idr1 = 1
          IF (irr.EQ.2) idr2 = 1

          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          ccom = ctmp(bwrd:ewrd)

          ResTyp(im,irr) = 0
          IF ((ModTyp(im).EQ.1).OR.(ModTyp(im).EQ.2)) THEN
            IF (FindStr(ccom,'app').GT.0) THEN
              ResTyp(im,irr) = 1
              IF (FindStr(ccom,'log').GT.0) THEN
                ResTyp(im,irr) = 3
              ENDIF
              IF (irr.EQ.1) ires1 = 1
              IF (irr.EQ.2) ires2 = 1
            ENDIF
            IF (FindStr(ccom,'phs').GT.0) THEN
              IF (FindStr(ccom,'deg').GT.0) THEN
                ResTyp(im,irr) = 2
              ENDIF
              IF (FindStr(ccom,'rad').GT.0) THEN
                ResTyp(im,irr) = 4
              ENDIF
              IF (irr.EQ.1) ires1 = 2
              IF (irr.EQ.2) ires2 = 2
            ENDIF
          ENDIF
          IF (ModTyp(im).EQ.3) THEN
            IF (FindStr(ccom,'rel').GT.0) THEN
               ResTyp(im,irr) = 5
               IF (irr.EQ.1) ires1 = 5
               IF (irr.EQ.2) ires2 = 5
            ENDIF
            IF (FindStr(ccom,'img').GT.0) THEN
               ResTyp(im,irr) = 6
               IF (irr.EQ.1) ires1 = 6
               IF (irr.EQ.2) ires2 = 6
            ENDIF
          ENDIF

          IF ((idr1.EQ.1).AND.(idr2.EQ.1)) THEN
            IF (ires1.EQ.ires2) THEN
              WRITE(6,2000)
              WRITE(6,*) '!!!!! Same response type are used !!!!!'
              WRITE(6,2010)
              STOP
            ENDIF
          ENDIF

          IF (ResTyp(im,irr).EQ.0) THEN
            WRITE(6,2000)
            IF ((ModTyp(im).EQ.1).OR.(ModTyp(im).EQ.2)) THEN
             WRITE(6,*) '!!! Response type (app/phsdeg/applog/phsrad)'
             WRITE(6,*) '!!! must be selected for TM or TE mode!!!'
            ENDIF
            IF (ModTyp(im).EQ.3) THEN
             WRITE(6,*) '!!! Response type ',
     >       ' (rel/img) must be selected for TP mode!!!'
            ENDIF
            WRITE(6,2010)
            STOP
          ENDIF

          IF (FindStr(ctmp(19:19),'*').GT.0) THEN
            DO ip = 1,NPer(im)
              READ(10,*,ERR=400) (DatRes(im,irr,ip,is),is=1,NSta(im))
            ENDDO ! ip 
          ELSE
            DO ip = 1,NPer(im)
              READ(10,*,ERR=400) 
     >            freq,(DatRes(im,irr,ip,is),is=1,NSta(im))
            ENDDO ! ip 
          ENDIF
          GOTO 410
400       CONTINUE
          WRITE(crr,'(i1)') irr
          ctext = 'DATA_RESPONSE_NO_'//crr//'  '
          CALL IncorrectInput(ctext)
          GOTO 410
410       CONTINUE
          GOTO 10
        ENDIF

        IF (FindStr(ctmp,'error_response_no_').GT.0) THEN
          IF (NRes(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_RESPONSE must be given ',
     >                  'prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          IF (NPer(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_PERIOD/NUMBER_OF_FREQUENCY',
     >       ' must be given prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          IF (NSta(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_STATION ',
     >       ' must be given prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          bwrd = BegWrd(ctmp,1)
          ewrd = EndWrd(ctmp,1)
          IF (FindStr(ctmp(ewrd:ewrd),'*').GT.0) THEN
            READ(ctmp(ewrd-1:ewrd-1),'(i1)') irr
          ELSE
            READ(ctmp(ewrd:ewrd),'(i1)') irr
          ENDIF

          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          READ(ctmp(bwrd:ewrd),*,ERR=450) ErrFlr(im,irr)
          IF ((irr.LT.1).OR.(irr.GT.NRes(im)))  GOTO 450
          IF ((ier1.EQ.1).AND.(irr.EQ.1)) GOTO 450
          IF ((ier2.EQ.1).AND.(irr.EQ.2)) GOTO 450
          IF (irr.EQ.1) ier1 = 1
          IF (irr.EQ.2) ier2 = 1

          IF (FindStr(ctmp(20:20),'*').GT.0) THEN
            DO ip = 1,NPer(im)
              READ(10,*,ERR=450) (DatErr(im,irr,ip,is),is=1,NSta(im))
            ENDDO ! ip 
          ELSE
            DO ip = 1,NPer(im)
              READ(10,*,ERR=450) 
     >          freq,(DatErr(im,irr,ip,is),is=1,NSta(im))
            ENDDO ! ip 
          ENDIF
          GOTO 460
450       CONTINUE
          WRITE(crr,'(i1)') irr
          ctext = 'ERROR_RESPONSE_NO_'//crr//'  '
          CALL IncorrectInput(ctext)
460       CONTINUE
          GOTO 10
        ENDIF

        IF (FindStr(ctmp,'data_inclusion_no_').GT.0) THEN
          IF (NRes(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_RESPONSE must be given ',
     >                  'prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          IF (NPer(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_PERIOD/NUMBER_OF_FREQUENCY',
     >       ' must be given prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF
          IF (NSta(im).EQ.0) THEN
             WRITE(6,2000)
             WRITE(6,*) '!!! NUMBER_OF_STATION ',
     >       ' must be given prior to DATA_RESPONSE_NO'
             WRITE(6,2010)
             STOP
          ENDIF

          bwrd = BegWrd(ctmp,1)
          ewrd = EndWrd(ctmp,1)
          READ(ctmp(ewrd:ewrd),'(i1)') irr
          WRITE(crr,'(i1)') irr
          ctext = 'DATA_INCLUSION_NO_'//crr//'  '
          IF ((irr.LT.1).OR.(irr.GT.NRes(im)))  GOTO 470
          IF ((iir1.EQ.1).AND.(irr.EQ.1)) GOTO 470
          IF ((iir2.EQ.1).AND.(irr.EQ.2)) GOTO 470
          IF (irr.EQ.1) iir1 = 1
          IF (irr.EQ.2) iir2 = 1

          bwrd = BegWrd(ctmp,2)
          ewrd = EndWrd(ctmp,2)
          ccom = ctmp(bwrd:ewrd)

          CALL ReadIndex(10,ccom,im,irr,NPer,NSta,DatInx,ctext)
          GOTO 480
470       CONTINUE
          CALL IncorrectInput(ctext)
480       CONTINUE
          GOTO 10
        ENDIF

        WRITE(6,2000)
        WRITE(6,2020) ctmp
        WRITE(6,2010)
        STOP

500     CONTINUE
      ENDDO


      iunit = 61

      CALL Lenb(COut,length)
      DO im = 1,NMode
        IF (ModTyp(im).EQ.1) cout2 = COut(1:length)//'.dix_tm'
        IF (ModTyp(im).EQ.2) cout2 = COut(1:length)//'.dix_te'
        IF (ModTyp(im).EQ.3) cout2 = COut(1:length)//'.dix_tp'
        OPEN(UNIT=iunit,FILE=cout2,STATUS='unknown')
        DO ir = 1,NRes(im)
          WRITE(crr,'(i1)') ir
          ctext = 'DATA_INCLUSION_NO_'//crr//' index '
          WRITE(iunit,*) ctext
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im) 
              WRITE(cindex(is:is),'(i1)') DatInx(im,ir,ip,is)
            ENDDO 
            WRITE(iunit,*) cindex(1:NSta(im))
          ENDDO ! ip 
        ENDDO ! ir
        CLOSE(iunit)
      ENDDO


C     Adjust Data

C     Increase the error from inclusion matrix 
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im) 
              IF (DatInx(im,ir,ip,is).GE.1) THEN
                incfac = D2**(DatInx(im,ir,ip,is)-1)
                DatErr(im,ir,ip,is) = incfac*DatErr(im,ir,ip,is)
                DatInx(im,ir,ip,is) = 1
              ELSE
                DatErr(im,ir,ip,is) = 999.99
              ENDIF
            ENDDO ! is 
          ENDDO ! ip
        ENDDO ! ir
      ENDDO ! im

C     Set minimum errors from error floor
C     Compare minimum errors with current errors
C     For app and phase of TM, TE, ErrFlr is relative error floor.
C     For real and imag of tipper, ErrFlr is the absolute error floor.
      DO im = 1,NMode
       DO ir = 1,NRes(im)

        IF (ResTyp(im,ir).EQ.1) THEN
          DO ip = 1,NPer(im)
           DO is = 1,NSta(im) 
            IF (DatInx(im,ir,ip,is).EQ.1) THEN
              minerr = ErrFlr(im,ir)*DatRes(im,ir,ip,is)
              DatErr(im,ir,ip,is) = MAX(DatErr(im,ir,ip,is),minerr)
            ENDIF
           ENDDO ! is 
          ENDDO ! ip
        ENDIF ! ResTyp = 1

        IF (ResTyp(im,ir).EQ.3) THEN
          eflr  = ErrFlr(im,ir)
          minerr = DLOG10(D1+eflr)
c         IF (LOGFILE_SCREEN.EQ.1) WRITE(6,1000) minerr
c         WRITE(99,1000) minerr
          DO ip = 1,NPer(im)
           DO is = 1,NSta(im) 
            IF (DatInx(im,ir,ip,is).EQ.1) THEN
              DatErr(im,ir,ip,is) = MAX(DatErr(im,ir,ip,is),minerr)
            ENDIF
           ENDDO ! is 
          ENDDO ! ip
        ENDIF

        IF (ResTyp(im,ir).EQ.2) THEN
          eflr = ErrFlr(im,ir)
          minerr = eflr*D90/PI
c         IF (LOGFILE_SCREEN.EQ.1) WRITE(6,1010) minerr
c         WRITE(99,1010) minerr
          DO ip = 1,NPer(im)
           DO is = 1,NSta(im) 
            IF (DatInx(im,ir,ip,is).EQ.1) THEN
              DatErr(im,ir,ip,is) = MAX(DatErr(im,ir,ip,is),minerr)
            ENDIF
           ENDDO ! is 
          ENDDO ! ip
        ENDIF

        IF (ResTyp(im,ir).EQ.4) THEN
          eflr = ErrFlr(im,ir)
          minerr = (eflr*D90/PI)*(PI/D180)
c         IF (LOGFILE_SCREEN.EQ.1) WRITE(6,1010) (minerr*(D180/PI))
c         WRITE(99,1010) (minerr*(D180/PI))
          DO ip = 1,NPer(im)
           DO is = 1,NSta(im) 
            IF (DatInx(im,ir,ip,is).EQ.1) THEN
              DatErr(im,ir,ip,is) = MAX(DatErr(im,ir,ip,is),minerr)
            ENDIF
           ENDDO ! is 
          ENDDO ! ip
        ENDIF

        IF ((ResTyp(im,ir).EQ.5).OR.(ResTyp(im,ir).EQ.6)) THEN
          minerr = ErrFlr(im,ir)
          DO ip = 1,NPer(im)
           DO is = 1,NSta(im) 
            IF (DatInx(im,ir,ip,is).EQ.1) THEN
              DatErr(im,ir,ip,is) = MAX(DatErr(im,ir,ip,is),minerr)
            ENDIF
           ENDDO ! is 
          ENDDO ! ip
        ENDIF
       ENDDO ! ir
      ENDDO ! im


C     Convert app. resitivity to log10(app) and phs (rad) to degree
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          IF (ResTyp(im,ir).EQ.1) THEN
            DO ip = 1,NPer(im)
              DO is = 1,NSta(im) 
c               IF (DatInx(im,ir,ip,is).EQ.1) THEN
                  DatErr(im,ir,ip,is) = 
     >            DatErr(im,ir,ip,is)/(DatRes(im,ir,ip,is)*DLOG(D10))
                  DatRes(im,ir,ip,is) = DLOG10(DatRes(im,ir,ip,is))
c               ENDIF
              ENDDO ! is 
            ENDDO ! ip
          ENDIF

          IF (ResTyp(im,ir).EQ.2) THEN
            DO ip = 1,NPer(im)
              DO is = 1,NSta(im) 
c               IF (DatInx(im,ir,ip,is).EQ.1) THEN
                  DatRes(im,ir,ip,is) = 
     >            -D1*DatRes(im,ir,ip,is)*(PI/D180)
                  DatErr(im,ir,ip,is) = DatErr(im,ir,ip,is)*(PI/D180)
c               ENDIF
              ENDDO ! is 
            ENDDO ! ip
          ENDIF

          IF (ResTyp(im,ir).EQ.3) THEN
c           no conversion
          ENDIF
          IF (ResTyp(im,ir).EQ.4) THEN
c           no conversion
          ENDIF
          IF (ResTyp(im,ir).EQ.5) THEN
c           no conversion
          ENDIF
          IF (ResTyp(im,ir).EQ.6) THEN
c           no conversion
          ENDIF

        ENDDO ! ir
      ENDDO ! im

C     Calculate Cd matrix
      idx = 0
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im) 
              IF (DatInx(im,ir,ip,is).EQ.1) THEN
                idx = idx + 1
                Cd(idx) =  DatErr(im,ir,ip,is)**D2
              ENDIF
            ENDDO ! is 
          ENDDO ! ip
        ENDDO ! ir
      ENDDO ! im


1000  FORMAT('*** Absolute Error Floor of ',f5.2,' is assigned for ',
     >       'the log10 app resistivity ***')
1010  FORMAT('*** Absolute Error Floor of ',f5.2,' deg is assigned ',
     >       ' for the phase ***')

2000  FORMAT('!!!!! ATTENTION, ERROR IN DATA FILE !!!!!')
2010  FORMAT('!!!!!  Please, correct and restart  !!!!!') 
2020  FORMAT('!!! Keyword ',a20,' is not recognized !!!')

5000  FORMAT(a80)


      RETURN
      END ! SUBROUTINE ReadData()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE ReadIndex(io,ccom,im,ir,NPer,NSta,Inx,ctext)
      INCLUDE 'parameter.h'

      INTEGER io,im,ir,NPer(*),NSta(*)
      INTEGER Inx(NMODMX,NRESMX,NPERMX,NSTAMX)
      CHARACTER*80 ccom
      CHARACTER*20 ctext

      INTEGER ip,is,len,nst,bwrd,ewrd,ist
      CHARACTER*256 cindex 
      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd

      IF (FindStr(ccom,'all').GT.0) THEN
        DO ip = 1,NPer(im)
          DO is = 1,NSta(im)
            Inx(im,ir,ip,is) = 1
          ENDDO ! is
        ENDDO ! ip
        GOTO 560
      ENDIF
      IF (FindStr(ccom,'index').GT.0) THEN
        DO ip = 1,NPer(im)
          nst = 0
500       CONTINUE
          READ(io,'(a256)') cindex
          bwrd = BegWrd(cindex,1)
          ewrd = EndWrd(cindex,1)
          len  = ewrd - bwrd + 1
          ist = 0
          DO is = bwrd,ewrd
           ist = ist + 1
           READ(cindex(is:is),'(i1)',ERR=550) Inx(im,ir,ip,ist+nst)
          ENDDO
          nst = len + nst
          IF (nst.LT.NSta(im)) GOTO 500
        ENDDO ! ip 
        GOTO 560
      ENDIF

550   CONTINUE
      CALL IncorrectInput(ctext)
560   CONTINUE

      END ! ReadIndex

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C     ModGrd = 1:  normal
C     ModGrd = 0:  freez the block
C            = 2:  no diffusion to the right
C            = 3:  no diffusion to the left
C            = 4:  no diffusion to the top
C            = 5:  no diffusion to the bottom
C            = 6:  no diffusion to the right and top
C            = 7:  no diffusion to the left  and top
C            = 8:  no diffusion to the right and bottom
C            = 9:  no diffusion to the left  and bottom

      SUBROUTINE ReadModelControl(CModelControl,Ny0,Nzb,
     >           DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,
     >           ModGrd,DFStatus)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      CHARACTER*80 CModelControl
      INTEGER Ny0,Nzb,ModGrd(NZ0MX,NY0MX),DTIME_STEP,DFStatus(*)
      REAL*8  DLENGTH_HOR,DLENGTH_VER

      CHARACTER*80 ctmp,ccom
      CHARACTER*256 cindex
      INTEGER iy,iz,nyt,len,ic,ewrd,bwrd,maxcommand,iyt
      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd


      DTIME_STEP = 1
      DLENGTH_HOR  = D0
      DLENGTH_VER  = D0
      CALL ConstantMatrixI4(ModGrd,NZ0MX,NY0MX,Nzb,Ny0,1)

      IF (CModelControl.EQ.'default') THEN 
C       Diffusion Model Covariance
        CALL ConstantMatrixI4(ModGrd,NZ0MX,NY0MX,Nzb,Ny0,1)
        GOTO 6000
      ELSE  
        OPEN(UNIT=10,FILE=CModelControl,STATUS='OLD')
        GOTO 99
      ENDIF

99    CONTINUE
      maxcommand = 4
      ic = 0
100   CONTINUE
      IF (ic.EQ.maxcommand) THEN
        GOTO 6000
      ENDIF
      READ(10,'(a80)',END=6000) ctmp
      ic = ic + 1

c     blank line
      IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
         ic = 0
         GOTO 100
      ENDIF
c     comment line
      IF (FindStr(ctmp(1:10),'#').GT.0) THEN
         ic = 0
         GOTO 100
      ENDIF

      IF (FindStr(ctmp,'time_step').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        IF (FindStr(ccom,'default').GT.0) THEN
          GOTO 100
        ELSE
          READ(ccom,*,ERR=200) DTIME_STEP
          IF (DTIME_STEP.LE.0) GOTO 200
          GOTO 100
200       CALL IncorrectInput('time_step           ')
        ENDIF
      ENDIF

      IF (FindStr(ctmp,'hor_length').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        IF (FindStr(ccom,'default').GT.0) THEN
          GOTO 100
        ELSE
          READ(ccom,*,ERR=300) DLENGTH_HOR
          IF (DLENGTH_HOR.LE.0) GOTO 300
          GOTO 100
300       CALL IncorrectInput('hor_length           ')
        ENDIF
      ENDIF

      IF (FindStr(ctmp,'ver_length').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        IF (FindStr(ccom,'default').GT.0) THEN
          GOTO 100
        ELSE
          READ(ccom,*,ERR=400) DLENGTH_VER
          IF (DLENGTH_VER.LE.0) GOTO 400
          GOTO 100
400       CALL IncorrectInput('ver_length           ')
        ENDIF
      ENDIF


      IF (FindStr(ctmp,'model_control').GT.0) THEN
        bwrd = BegWrd(ctmp,2)
        ewrd = EndWrd(ctmp,2)
        ccom = ctmp(bwrd:ewrd)
        IF (FindStr(ccom,'default').GT.0) THEN
          CALL ConstantMatrixI4(ModGrd,NZ0MX,NY0MX,Nzb,Ny0,1)
          GOTO 100
        ENDIF
        IF (FindStr(ccom,'index').GT.0) THEN
          DO iz = 1,Nzb
            nyt = 0
700         CONTINUE
            READ(10,'(a256)') cindex
            bwrd = BegWrd(cindex,1)
            ewrd = EndWrd(cindex,1)
            len  = ewrd - bwrd + 1
            iyt = 0
            DO iy = bwrd,ewrd
              iyt = iyt + 1
              READ(cindex(iy:iy),'(i1)',ERR=800) ModGrd(iz,iyt+nyt)
            ENDDO
            nyt = len + nyt
            IF (nyt.LT.NY0) GOTO 700
          ENDDO
          GOTO 810
800       CONTINUE
          CALL IncorrectInput('model_control       ')
810       CONTINUE
          GOTO 100
        ENDIF

        WRITE(6,8000)
        WRITE(6,8020) ctmp
        WRITE(6,8010)
        STOP
      ENDIF

      WRITE(6,8000)
      WRITE(6,8020) ctmp
      WRITE(6,8010)
      STOP

6000  CONTINUE
      CLOSE(10)

      CALL CheckStatus(Nzb,Ny0,ModGrd,DFStatus)

8000  FORMAT('!!!!!! ATTENTION, ERROR IN THE MODEL CONTROL FILE !!!!!')
8010  FORMAT('!!!!!!      please, correct and restart     !!!!!!')
8020  FORMAT('!!! Keyword ',a15,' is not recognized !!!!!')

      RETURN
      END ! SUBROUTINE ReadModelGrid()
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE AssignValIndx(ix1,ix2,jx1,jx2,np1,np2,midx,val,mout)
      INTEGER ix1,ix2,jx1,jx2,np1,np2
      INTEGER midx(np1,np2)
      REAL*8  mout(np1,np2),val(*)

      INTEGER i,j

      DO i = ix1,ix2
        DO j = jx1,jx2
          mout(i,j) = val(midx(i,j))
        ENDDO
      ENDDO

      RETURN
      END ! SUBROUTINE  AssignValIndx

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE CheckStatus(Mz,My,MGrd,DFStatus)
      INCLUDE 'parameter.h'

      INTEGER My,Mz
      INTEGER MGrd(NZ0MX,NY0MX),DFStatus(*)

      INTEGER iy,iz,jj

C     Diffusion status of each block 
C     DFStatus = 1; if that block allow to be diffused
C              = 0; if that block is not allowed to be diffused.

      jj = 0
      DO iy = 1,My
        DO iz = 1,Mz
          jj = jj + 1
          IF (MGrd(iz,iy).EQ.0) THEN
            DFStatus(jj) = 0
          ELSE 
            DFStatus(jj) = 1
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END ! SUBROUTINE  CheckStatus

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

C     Sensitivity inclusion matrix=> SenInx(im,ir,ip,is)

      SUBROUTINE ReadSenInc(CSenInc,NMode,NRes,NPer,NSta,DatInx,
     >           SenInx,NNT,NN,COut,ModTyp)
      INCLUDE 'parameter.h'

      CHARACTER*80 CSenInc(NMODMX),COut
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*)
      INTEGER SenInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER NNT(*),NN(NMODMX,NRESMX,2)

      INTEGER im,ir,ip,is,nsum,pp,ss,bwrd,ewrd
      INTEGER stripe(NPERMX),endstp,ipd,isw,iflag,length,
     >        irr,iir1,iir2,ic,iunit
      CHARACTER*1  crr
      CHARACTER*30 ctext
      CHARACTER*80 ctmp,csn,cout2,ccom
      CHARACTER*256 cindex

      INTEGER  SumIndx,FindStr,BegWrd,EndWrd
      EXTERNAL SumIndx,FindStr,BegWrd,EndWrd

      DO im = 1,NMode
        csn = CSenInc(im)
        IF ((csn(1:6).EQ.'stripe').OR.(csn(1:7).EQ.'checker')) THEN
C         rule :  first and last periods are alway selected for good interpolation.
C              :  start from long period (more smoother) to short period
          IF (csn(1:6).EQ.'stripe') THEN
            READ(csn(10:11),'(i2)') pp
  
            IF (pp.LE.1) THEN
              WRITE(6,2000) 
              WRITE(6,*) '!!!   This is a SUB-SPACE version,',
     >                   ' p should be more than 1 !!!'
              WRITE(6,*) '!!!   Want full basis version, Please ',
     >                   ' use DASOCC inversion !!!'
              WRITE(6,2010)
              STOP
            ENDIF
            IF (pp.GE.10) THEN
              WRITE(6,*) 
     >   '????           Warning: Few periods used to compute   ????'
              WRITE(6,*)
     >   '???? representer may make the inversion not converging ????'
              WRITE(6,*) 
     >   '????           Two periods (at least) per decade to   ????'
              WRITE(6,*)
     >   '???? form representer is recommended for best results  ????'                     
            ENDIF
            IF (pp.GT.NPer(im)) THEN
              WRITE(6,2000) 
              WRITE(6,*) '!!! p value cannot exceed NUMBER OF PERIOD'
              WRITE(6,2010)
              STOP
            ENDIF

            CALL ConstantVectorI4(stripe,NPer(im),0)

            IF ((NPer(im)-(NPer(im)/pp)*pp).EQ.1) THEN
              endstp = NPer(im)-1
            ELSE
              endstp = ((NPer(im)/pp) - 1)*pp
            ENDIF

            DO ip = NPER(im),NPer(im)-endstp,-pp
              stripe(ip) = 1
            ENDDO
            ipd = (NPer(im)-endstp)/2 +  mod(Nper(im)-endstp,2) 
            stripe(ipd) = 1
            stripe(1) = 1

            DO ir = 1,NRes(im)
              DO ip = 1,NPER(im)
                IF (stripe(ip).EQ.1) THEN
                  DO is = 1,NSta(im)
                    SenInx(im,ir,ip,is) = 1
                  ENDDO
                ELSE
                  DO is = 1,NSta(im)
                    SenInx(im,ir,ip,is) = 0
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

          ENDIF
          IF (csn(1:7).EQ.'checker') THEN
            isw = FindStr(csn(9:20),':')
            READ(csn(11:7+isw),*) pp
            READ(csn(8+isw+3:8+isw+4),*) ss

            IF (pp.LE.1) THEN
              WRITE(6,2000) 
              WRITE(6,*) '!!!   This is a SUB-SPACE version,',
     >                   ' p should be more than 1 !!!'
              WRITE(6,*) '!!!   Want full basis version, Please ',
     >                   ' use DASOCC inversion !!!'
              WRITE(6,2010)
              STOP
            ENDIF
            IF (ss.LE.0) THEN
              WRITE(6,2000) 
              WRITE(6,*) '!!!  s must be positive number !!!' 
              WRITE(6,2010)
              STOP
            ENDIF
            IF (pp.GE.10) THEN
              WRITE(6,*) 
     >   '????           Warning: Few periods used to compute   ????'
              WRITE(6,*)
     >   '???? representer may make the inversion not converging ????'
              WRITE(6,*) 
     >   '????           Two periods (at least) per decade to   ????'
              WRITE(6,*)
     >   '???? form representer is recommended for best results  ????'                     
            ENDIF
            IF (ss.GE.4) THEN
              WRITE(6,*) 
     >   '????           Warning: Few stations used to compute   ????'
              WRITE(6,*)
     >   '???? representer may make the inversion not converging ????'
            ENDIF
            IF (pp.GT.NPer(im)) THEN
              WRITE(6,2000) 
              WRITE(6,*) '!!! p value cannot exceed NUMBER OF PERIOD'
              WRITE(6,2010)
              STOP
            ENDIF
            IF (ss.GT.NSta(im)) THEN
              WRITE(6,2000) 
              WRITE(6,*) '!!! s value cannot exceed NUMBER OF STATION'
              WRITE(6,2010)
              STOP
            ENDIF

            CALL ConstantVectorI4(stripe,NPer(im),0)

            IF ((NPer(im)-(NPer(im)/pp)*pp).EQ.1) THEN
              endstp = NPer(im)-1
            ELSE
              endstp = ((NPer(im)/pp) - 1)*pp
            ENDIF

            DO ip = NPER(im),NPer(im)-endstp,-pp
              stripe(ip) = 1
            ENDDO
            ipd = (NPer(im)-endstp)/2 +  mod(Nper(im)-endstp,2) 
            stripe(ipd) = 1
            stripe(1) = 1

            DO ir = 1,NRes(im)
              iflag = 0
              DO ip = 1,NPER(im)
                IF (stripe(ip).EQ.1) THEN
                  DO is = 1,NSta(im)
                    SenInx(im,ir,ip,is) = 0
                  ENDDO
                  IF (iflag.EQ.0) THEN
                    DO is = 1,NSta(im),ss
                      SenInx(im,ir,ip,is) = 1
                    ENDDO
                    iflag = 1
                  ELSE
                    DO is = 1+ss/2,NSta(im),ss
                      SenInx(im,ir,ip,is) = 1
                    ENDDO
                    iflag = 0
                  ENDIF
                ELSE
                  DO is = 1,NSta(im)
                    SenInx(im,ir,ip,is) = 0
                  ENDDO
                ENDIF
              ENDDO
            ENDDO

          ENDIF
        ELSE
          OPEN(UNIT=10,FILE=CSenInc(im),STATUS='OLD')
          iir1 = 0
          iir2 = 0
          ic   = 0

100       CONTINUE

          IF (ic.EQ.NRes(im)) THEN
            CLOSE(10)
            GOTO 500
          ENDIF

          READ(10,'(a80)',END=200) ctmp
          GOTO 300

200       IF (ic.LT.NRes(im)) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!! End of file reached before completed !!!'
            WRITE(6,2010)
            STOP
          ENDIF

300       CONTINUE

c         blank line
           IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
             GOTO 100
          ENDIF
c         comment line
          IF (FindStr(ctmp(1:10),'#').GT.0) THEN
             GOTO 100
          ENDIF
          ic = ic+1
       
          IF (FindStr(ctmp,'sens_inclusion_no_').GT.0) THEN
            bwrd = BegWrd(ctmp,1)
            ewrd = EndWrd(ctmp,1)
            READ(ctmp(ewrd:ewrd),'(i1)') irr
            WRITE(crr,'(i1)') irr
            ctext = 'SENS_INCLUSION_NO_'//crr//'  '
            IF ((irr.LT.1).OR.(irr.GT.NRes(im)))  GOTO 470
            IF ((iir1.EQ.1).AND.(irr.EQ.1)) GOTO 470
            IF ((iir2.EQ.1).AND.(irr.EQ.2)) GOTO 470
            IF (irr.EQ.1) iir1 = 1
            IF (irr.EQ.2) iir2 = 1

            bwrd = BegWrd(ctmp,2)
            ewrd = EndWrd(ctmp,2)
            ccom = ctmp(bwrd:ewrd)

            CALL ReadIndex(10,ccom,im,irr,NPer,NSta,SenInx,ctext)
            GOTO 480
470         CONTINUE
            CALL IncorrectInput(ctext)
480         CONTINUE
            GOTO 100
          ENDIF

          WRITE(6,2000)
          WRITE(6,2020) ctmp
          WRITE(6,2010)
          STOP

500       CONTINUE

          CLOSE(10)
        ENDIF
      ENDDO ! im

C     Total Number of data                               => NNT(1)
C     Number of data selected in DatInx of each responds => NN(:,:,1),NNT(2)
C     Number of data selected in SenInx of each responds => NN(:,:,2),NNT(3)

      nsum = 0
      DO im = 1,NMode
        nsum = nsum + NRes(im)*NPer(im)*NSta(im)
      ENDDO ! im
      NNT(1) = nsum

      NNT(2) = 0
      NNT(3) = 0
      DO im = 1,NMode
        DO ir = 1,NRes(im)
          NN(im,ir,1) = SumIndx(DatInx,NMODMX,NRESMX,NPERMX,NSTAMX,
     >                          im,ir,1,1,im,ir,NPer(im),NSta(im))
          NNT(2) = NNT(2)+NN(im,ir,1)
          NN(im,ir,2) = SumIndx(SenInx,NMODMX,NRESMX,NPERMX,NSTAMX,
     >                          im,ir,1,1,im,ir,NPer(im),NSta(im))
          NNT(3) = NNT(3)+NN(im,ir,2)
        ENDDO ! ir
      ENDDO ! im
      
      iunit = 61

      CALL Lenb(COut,length)
      DO im = 1,NMode
        IF (ModTyp(im).EQ.1) cout2 = COut(1:length)//'.six_tm'
        IF (ModTyp(im).EQ.2) cout2 = COut(1:length)//'.six_te'
        IF (ModTyp(im).EQ.3) cout2 = COut(1:length)//'.six_tp'
        OPEN(UNIT=iunit,FILE=cout2,STATUS='unknown')
        DO ir = 1,NRes(im)
          WRITE(crr,'(i1)') ir
          ctext = 'SENS_INCLUSION_NO_'//crr//' index '
          WRITE(iunit,*) ctext
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im) 
              WRITE(cindex(is:is),'(i1)') SenInx(im,ir,ip,is)
            ENDDO 
            WRITE(iunit,*) cindex(1:NSta(im))
          ENDDO ! ip 
        ENDDO ! ir
        CLOSE(iunit)
      ENDDO

2000  FORMAT('!!!!! ATTENTION, ERROR IN SENS INCLUSION FILE !!!!!')
2010  FORMAT('!!!!!  Please, correct and restart !!!!!')
2020  FORMAT('!!! Keyword ',a20,' is not recognized !!!')
      RETURN
      END ! ReadSenInc()
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE SetSknDepth(BackGroundRho,NMode,NPer,Period,SknDepth)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER NMode,NPer(*)
      REAL*8  BackGroundRho,SknDepth(NMODMX,NPERMX),
     >        Period(NMODMX,NPERMX)
C     for later version where skin depth are different from frequency
C     to frequency. 

      INTEGER im,ip
      REAL*8  omega,omue
 
C     calculate skin depth of each period 
C     assume that it is the distance of that data point(period,station)
C     from the surface.

      DO im = 1,NMode
        DO ip = 1,NPer(im)
          omega = (D2*PI)/Period(im,ip)
          omue  = omega*Mue
          SknDepth(im,ip) = DSQRT(D2*BackGroundRho/omue)
        ENDDO ! ip
      ENDDO ! im

     
      RETURN
      END ! SetSknDepth()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C     PART 7 : Static shift indices        => SSIndx(im,is),
C            : Static shift paras          => SSPara(im,is)

      SUBROUTINE ReadStaticInc(CDist,NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >           DatInx,SSIndx,SSPara,StsInx)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      CHARACTER*80 CDist(NMODMX)
      INTEGER NMode,NPer(*),NSta(*),NRes(*),ModTyp(*)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DatInx(NMODMX,NRESMX,NPERMX,NSTAMX)
      INTEGER SSIndx(NMODMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      REAL*8  SSPara(NMODMX,NSTAMX)


      INTEGER im,ir,ip,is,nst,len,maxcommand,bwrd,ewrd,ic,ist
      CHARACTER*80 ctmp,ccom
      CHARACTER*256 cindex

      INTEGER  FindStr,BegWrd,EndWrd
      EXTERNAL FindStr,BegWrd,EndWrd


      DO im = 1,NMode
        IF (CDist(im).EQ.'default') THEN
          DO is = 1,NSta(im)
            SSIndx(im,is) = 0
            SSPara(im,is) = D0
          ENDDO ! is
          DO ip = 1,NPer(im)
            DO is = 1,NSta(im)
              StsInx(im,ip,is) = 1
            ENDDO ! is
          ENDDO ! ip
        ELSE

          IF (ModTyp(im).EQ.3) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!!! No static distortion file for TP !!!!'
            WRITE(6,2010)
            STOP
          ENDIF
 
          OPEN(UNIT=10,FILE=CDist(im),STATUS='OLD')
          maxcommand = 3
          ic = 0

10        CONTINUE

          IF (ic.EQ.maxcommand) THEN
            CLOSE(10)
            GOTO 500
          ENDIF

          READ(10,'(a80)',END=20) ctmp
          GOTO 30

20        IF (ic.LT.maxcommand) THEN
            WRITE(6,2000)
            WRITE(6,*) '!!! End of file reached before completed !!!'
            WRITE(6,2010)
            STOP
          ENDIF

30        CONTINUE

c         blank line
          IF (FindStr(ctmp(1:10),'          ').GT.0) THEN
             GOTO 10
          ENDIF
c         comment line
          IF (FindStr(ctmp(1:10),'#').GT.0) THEN
             GOTO 10
          ENDIF
          ic = ic+1

          IF (FindStr(ctmp,'distortion_index').GT.0) THEN
            nst = 0
100         CONTINUE
            READ(10,'(a256)') cindex
            bwrd = BegWrd(cindex,1)
            ewrd = EndWrd(cindex,1)
            len  = ewrd - bwrd + 1
            ist = 0
            DO is = bwrd,ewrd
              ist = ist + 1
              READ(cindex(is:is),'(i1)',ERR=110) SSIndx(im,ist+nst)
            ENDDO
            nst = len + nst
            IF (nst.LT.NSta(im)) GOTO 100
            GOTO 120
110         CONTINUE
            CALL IncorrectInput('DISTORTION_INDEX    ')
120         CONTINUE
            GOTO 10
          ENDIF

          IF (FindStr(ctmp,'distortion_parameter').GT.0) THEN
            READ(10,*,ERR=150) (SSPara(im,is),is=1,NSta(im))
            GOTO 160
150         CONTINUE
            CALL IncorrectInput('DISTORTION_PARAMETER')
160         CONTINUE
            GOTO 10
          ENDIF
       
          IF (FindStr(ctmp,'distortion_inclusion').GT.0) THEN
            bwrd = BegWrd(ctmp,2)
            ewrd = EndWrd(ctmp,2)
            ccom = ctmp(bwrd:ewrd)

            IF (FindStr(ccom,'all').GT.0) THEN
              DO ip = 1,NPer(im)
                DO is = 1,NSta(im)
                  StsInx(im,ip,is) = 1
                ENDDO ! is
              ENDDO ! ip
              GOTO 260
            ENDIF
            IF (FindStr(ccom,'index').GT.0) THEN
             DO ip = 1,NPer(im)
              nst = 0
200           CONTINUE
              READ(10,'(a256)') cindex
              bwrd = BegWrd(cindex,1)
              ewrd = EndWrd(cindex,1)
              len  = ewrd - bwrd + 1
              ist  = 0
              DO is = bwrd,ewrd
               ist = ist + 1
               READ(cindex(is:is),'(i1)',ERR=250) StsInx(im,ip,ist+nst)
              ENDDO
              nst = len + nst
              IF (nst.LT.NSta(im)) GOTO 200
             ENDDO ! ip 
             GOTO 260
            ENDIF
            WRITE(6,2000)
            WRITE(6,2020) ctmp
            WRITE(6,2010)
            STOP
250         CONTINUE
            CALL IncorrectInput('DISTORTION_INCLUSION')
260         CONTINUE
            GOTO 10
          ENDIF

          WRITE(6,2000)
          WRITE(6,2020) ctmp
          WRITE(6,2010)
          STOP

500       CONTINUE
          CLOSE(10)
        ENDIF
      ENDDO ! im

C     Requirements and Rules for StsInx : 
C       - Wherever DatInx is omitted, StsInx is automatically omitted too.

      DO im = 1,NMode
        DO ir = 1,NRes(im) 
          IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN
            DO ip = 1,NPer(im)
              DO is = 1,NSta(im) 
                IF (DatInx(im,ir,ip,is).EQ.0) THEN
                  StsInx(im,ip,is) = 0
                ENDIF
              ENDDO ! is 
            ENDDO ! ip
            DO is = 1,NSta(im)
              IF (SSIndx(im,is).EQ.1) THEN
                nst = 0
                DO ip = 1,NPer(im)
                  nst = nst + StsInx(im,ip,is)
                ENDDO !ip
                IF (nst.EQ.0) THEN
                  WRITE(6,2000)
                  WRITE(6,2050) is
                  WRITE(6,2060) is
                  WRITE(6,2010)
                  STOP
                ENDIF
              ENDIF
            ENDDO !is 
          ENDIF
        ENDDO ! ir
      ENDDO ! im

2000  FORMAT('!!!!! ATTENTION, ERROR IN THE DISTORTION FILE !!!!!')
2010  FORMAT('!!!!!  Please, correct and restart  !!!!!') 
2020  FORMAT('!!! Keyword ',a20,' is not recognized !!!')

2050  FORMAT('!!! No data are selected at station ',i3,
     >       ' in distortion_inclusion !!!')
2060  FORMAT('!!! or  distortion_index at station ',i3,
     >       ' must be turned off !!!')

      RETURN
      END ! ReadStaticInc()
C                                                                    C
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

C     Smoothness Control File Format
C 
      SUBROUTINE SetSmooth(NMode,NSta,Nzb,Ny,ZDis,YDis,StaNod,StaPos,
     >           DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,HGamma,VGamma)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER DTIME_STEP,NMode,NSta(*),Ny,Nzb
      INTEGER StaNod(NMODMX,NSTAMX)
      REAL*8  StaPos(NMODMX,NSTAMX)
      REAL*8  DLENGTH_HOR,DLENGTH_VER,ZDis(*),YDis(*)
      REAL*8  HGamma(NZ0MX,NY0MX),VGamma(*),HXe,VXe

      INTEGER is,iz,iy,nst,snod(NMODMX*NSTAMX)
      REAL*8  hfx(NY0MX),staspace(NMODMX*NSTAMX)
      REAL*8  spos(NMODMX*NSTAMX)

C     Xe is the distance where diffusion energy drop by 1/e
C     Xe = sqrt(4*Gamma*2*DTIME_STEP) where Gamma is diffusion parameter
c     Gamma = (Xe**D2)/(D4*D2*DTIME_STEP)

      IF (DLENGTH_HOR.NE.D0) THEN
        DO iz = 1,Nzb
          DO iy = 1,Ny
            HGamma(iz,iy) = (DLENGTH_HOR**D2)/(D8*DTIME_STEP)
          ENDDO ! iy
        ENDDO ! iz
      ELSE
C       Determine total no. of distinct stations of every modes
C       and sort them out in ascending order
        CALL FindDistinctStation(NMode,NSta,StaPos,nst,spos)

        DO is = 1,nst-1
          staspace(is) = spos(is+1) - spos(is)
        ENDDO ! is
        DO is = 1,nst
          snod(is) = 0
          DO iy = 1,Ny
            IF (spos(is).EQ.YDis(iy)) THEN
              snod(is) = iy
              GOTO 200
            ENDIF
          ENDDO
200       CONTINUE
          IF (snod(is).EQ.0) THEN
            WRITE(6,*) 'ERROR snod '
            STOP
          ENDIF
        ENDDO ! is

        is = 1
        DO iy = 1,snod(is)-1
          hfx(iy) = MIN(staspace(is),YDis(snod(is)))
        ENDDO ! iy
        DO is = 2,nst
          DO iy = snod(is-1),snod(is)-1
            hfx(iy) = staspace(is-1)
          ENDDO ! iy
        ENDDO
        is = nst
        DO iy = snod(is),Ny
          hfx(iy) = MIN(staspace(is-1),
     >                  (YDis(Ny+1)-YDis(snod(is))))
        ENDDO ! iy

        DO iz = 1,Nzb
          DO iy = 1,Ny
            HXe  = MAX(hfx(iy),ZDis(iz+1))
            HGamma(iz,iy) = (HXe**D2)/(D8*DTIME_STEP)
          ENDDO ! iy
        ENDDO ! iz
      ENDIF


      IF (DLENGTH_VER.NE.D0) THEN
        DO iz = 1,Nzb
          VGamma(iz) = (DLENGTH_VER**D2)/(D8*DTIME_STEP)
        ENDDO ! iz
      ELSE
        DO iz = 1,Nzb
          VXe        = ZDis(iz+1)
          VGamma(iz) = (Vxe**D2)/(D8*DTIME_STEP)
        ENDDO ! iz
      ENDIF

      RETURN
      END ! SUBROUTINE SetSmooth

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
