
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE WriteModel(imodel,iter,Tol,MNorm,NMode,
     >           Nza,Nzb,Ny,Dza,Dzb,Dy,Rho)
      INCLUDE 'parameter.h'
 
      INTEGER Nzb,Nza,Ny,imodel,iter,NMode
      REAL*8  Rho(NZ0MX,NY0MX),Tol,MNorm 
      REAL*8  DY(NY0MX),DZB(NZ0MX),DZA(NZ0MX)
 
      INTEGER iy,iz
 
      IF (iter.EQ.-1) WRITE(imodel,303) 
      IF (iter.EQ.0)  WRITE(imodel,305) Tol
      IF (iter.GE.1)  WRITE(imodel,310) iter,Tol,MNorm
      WRITE(imodel,410) Ny
      WRITE(imodel,200) (Dy(iy),iy=1,Ny)

      WRITE(imodel,420) Nzb
      WRITE(imodel,200) (Dzb(iz),iz=1,Nzb)

      WRITE(imodel,430) Nza
      WRITE(imodel,200) (Dza(iz),iz=1,Nza)

      WRITE(imodel,500)
      DO iz = 1,Nzb
        WRITE(imodel,205) (Rho(iz,iy),iy=1,Ny)
      ENDDO
      CLOSE(imodel)

 
200   FORMAT(7F12.4)
205   FORMAT(7E12.4)
303   FORMAT('TITLE PRIOR_MODEL')
305   FORMAT('TITLE STARTING_MODEL_RMS= ',F12.4)
310   FORMAT('TITLE MODEL_ID#',I2,'_RMS= ',F12.4,' MODEL_NORM= ',F12.4)
400   FORMAT(7I4)
410   FORMAT('NY  ',i3)
420   FORMAT('NZB ',i3)
430   FORMAT('NZA ',i3)
500   FORMAT('RESISTIVITY_MODEL ')

999   CONTINUE

      RETURN
      END ! WriteModel

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE WriteResponse(iter,COut,num,
     >           NMode,NRes,NPer,NSta,FF,
     >           ModTyp,ResTyp,Period,StaPos,RMS)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*),
     >        ResTyp(NMODMX,NRESMX),imodel,iter
      REAL*8  FF(NMODMX,NRESMX,NPERMX,NSTAMX),Period(NMODMX,NPERMX),
     >        StaPos(NMODMX,NSTAMX),RMS(*)
      CHARACTER*3  num
      CHARACTER*80 COut,cout2
 
      INTEGER im,ir,ip,is,iss,isa,isb,ic,ncol,nc2,len
      REAL*8  per

      imodel = 31
      CALL Lenb(COut,len)

      DO im = 1,NMode
        IF (iter.EQ.-1) THEN
          IF (ModTyp(im).EQ.1) cout2 = COut(1:len)//'.data_tm'
          IF (ModTyp(im).EQ.2) cout2 = COut(1:len)//'.data_te'
          IF (ModTyp(im).EQ.3) cout2 = COut(1:len)//'.data_tp'
          GOTO 10
        ENDIF
        IF (iter.EQ.-2) THEN
          IF (ModTyp(im).EQ.1) cout2 = COut(1:len)//'.error_tm'
          IF (ModTyp(im).EQ.2) cout2 = COut(1:len)//'.error_te'
          IF (ModTyp(im).EQ.3) cout2 = COut(1:len)//'.error_tp'
          GOTO 10
        ENDIF
        IF (ModTyp(im).EQ.1) cout2 = COut(1:len)//'_resp_tm.'//num
        IF (ModTyp(im).EQ.2) cout2 = COut(1:len)//'_resp_te.'//num
        IF (ModTyp(im).EQ.3) cout2 = COut(1:len)//'_resp_tp.'//num

10      CONTINUE
        OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')
        IF (iter.EQ.-2) WRITE(imodel,180)
        IF (iter.EQ.-1) WRITE(imodel,190)
        IF (iter.EQ.0) WRITE(imodel,200) RMS(1)
        IF (iter.GE.1) WRITE(imodel,210) iter,RMS(1)
        IF (ModTyp(im).EQ.1) WRITE(imodel,220)
        IF (ModTyp(im).EQ.2) WRITE(imodel,230)
        IF (ModTyp(im).EQ.3) WRITE(imodel,240)

        WRITE(imodel,250) NRes(im)
        WRITE(imodel,260) NPer(im)
        WRITE(imodel,100) (Period(im,ip),ip=1,NPer(im))
        WRITE(imodel,280) NSta(im)
        WRITE(imodel,101) (StaPos(im,is),is=1,NSta(im))

        nc2  = MOD(NSta(im),7)
        IF (nc2.EQ.0) THEN
          ncol = IDINT(NSta(im)/D7)
        ELSE
          ncol = IDINT(NSta(im)/D7)+1
        ENDIF

        DO ir = 1,NRes(im)
           IF (ResTyp(im,ir).EQ.1) THEN
             IF (ir.EQ.1) WRITE(imodel,300) 'app   '
             IF (ir.EQ.2) WRITE(imodel,310) 'app   '
           ENDIF
           IF (ResTyp(im,ir).EQ.3) THEN
             IF (ir.EQ.1) WRITE(imodel,300) 'applog'
             IF (ir.EQ.2) WRITE(imodel,310) 'applog'
           ENDIF
           IF (ResTyp(im,ir).EQ.2) THEN
             IF (ir.EQ.1) WRITE(imodel,300) 'phsdeg'
             IF (ir.EQ.2) WRITE(imodel,310) 'phsdeg'
           ENDIF
           IF (ResTyp(im,ir).EQ.4) THEN
             IF (ir.EQ.1) WRITE(imodel,300) 'phsrad'
             IF (ir.EQ.2) WRITE(imodel,310) 'phsrad'
           ENDIF
           IF (ResTyp(im,ir).EQ.5) THEN
             IF (ir.EQ.1) WRITE(imodel,300) 'rel   '
             IF (ir.EQ.2) WRITE(imodel,310) 'rel   '
           ENDIF
           IF (ResTyp(im,ir).EQ.6) THEN
             IF (ir.EQ.1) WRITE(imodel,300) 'img   '
             IF (ir.EQ.2) WRITE(imodel,310) 'img   '
           ENDIF

           DO ip = 1,NPer(im)
             iss = 0
             DO ic = 1,ncol
               isa = iss+1
               isb = iss+7
               IF (isb.GT.NSta(im)) isb = NSta(im)

               IF (ResTyp(im,ir).EQ.1) THEN
                 IF (ic.EQ.1) THEN
                   WRITE(imodel,105) Period(im,ip),
     >                    (D10**FF(im,ir,ip,is),is=isa,isb)
                 ELSE
                   WRITE(imodel,110) ' ', 
     >                    (D10**FF(im,ir,ip,is),is=isa,isb)
                 ENDIF
               ENDIF
               IF (ResTyp(im,ir).EQ.3) THEN
                 IF (ic.EQ.1) THEN
                   WRITE(imodel,105) Period(im,ip),
     >                    (FF(im,ir,ip,is),is=isa,isb)
                 ELSE
                   WRITE(imodel,110) ' ', 
     >                    (FF(im,ir,ip,is),is=isa,isb)
                 ENDIF
               ENDIF
               IF (ResTyp(im,ir).EQ.2) THEN
                 IF (ic.EQ.1) THEN
                   WRITE(imodel,105) Period(im,ip),
     >                (-D1*(D180/pi)*FF(im,ir,ip,is),is=isa,isb)
                 ELSE
                   WRITE(imodel,110) ' ', 
     >                (-D1*(D180/pi)*FF(im,ir,ip,is),is=isa,isb)
                 ENDIF
               ENDIF
               IF (ResTyp(im,ir).EQ.4) THEN
                 IF (ic.EQ.1) THEN
                   WRITE(imodel,105) Period(im,ip),
     >                (-D1*FF(im,ir,ip,is),is=isa,isb)
                 ELSE
                   WRITE(imodel,110) ' ', 
     >                (-D1*FF(im,ir,ip,is),is=isa,isb)
                 ENDIF
               ENDIF
               IF ((ResTyp(im,ir).EQ.5).OR.(ResTyp(im,ir).EQ.6)) THEN
                 IF (ic.EQ.1) THEN
                   per = Period(im,ip)
                   WRITE(imodel,105) per,(FF(im,ir,ip,is),is=isa,isb)
                 ELSE
                   WRITE(imodel,110) ' ',(FF(im,ir,ip,is),is=isa,isb)
                 ENDIF
               ENDIF

               iss = iss+7
             ENDDO ! ic
           ENDDO ! ip
        ENDDO ! ir
        CLOSE(imodel)
      ENDDO ! im
 
100   FORMAT(7F12.4)
101   FORMAT(7F12.2)
105   FORMAT(8F12.4)
110   FORMAT(a12,7F12.4)

180   FORMAT('TITLE  ABSOLUTE_ERROR')
190   FORMAT('TITLE  DATA')
200   FORMAT('TITLE  STARTING_MODEL_RMS_= ',F12.4)
210   FORMAT('TITLE  MODEL_ID#',I2,'_RMS= ',F12.4)
220   FORMAT('MODE_TYPE               tm')
230   FORMAT('MODE_TYPE               te')
240   FORMAT('MODE_TYPE               tp')

250   FORMAT('NUMBER_OF_RESPONSE       ',i1)
260   FORMAT('NUMBER_OF_PERIOD       ',i3)
280   FORMAT('NUMBER_OF_STATION      ',i3)

300   FORMAT('DATA_RESPONSE_NO_1      ',a6)
310   FORMAT('DATA_RESPONSE_NO_2      ',a6)
 
      RETURN
      END ! WriteRespond()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE WriteErrorResponse(iter,COut,
     >           NMode,NRes,NPer,NSta,EE,
     >           ModTyp,ResTyp,Period,StaPos,ErrFlr)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*),
     >        ResTyp(NMODMX,NRESMX),imodel,iter
      REAL*8  EE(NMODMX,NRESMX,NPERMX,NSTAMX),Period(NMODMX,NPERMX),
     >        StaPos(NMODMX,NSTAMX),ErrFlr(NMODMX,NRESMX)
      CHARACTER*80 COut,cout2
 
      INTEGER im,ir,ip,is,iss,isa,isb,ic,ncol,nc2,len

      imodel = 31
      CALL Lenb(COut,len)

      DO im = 1,NMode
        IF (iter.EQ.-2) THEN
          IF (ModTyp(im).EQ.1) cout2 = COut(1:len)//'.error_tm'
          IF (ModTyp(im).EQ.2) cout2 = COut(1:len)//'.error_te'
          IF (ModTyp(im).EQ.3) cout2 = COut(1:len)//'.error_tp'
          GOTO 10
        ENDIF

10      CONTINUE
        OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')
        IF (iter.EQ.-2) WRITE(imodel,180)
        IF (ModTyp(im).EQ.1) WRITE(imodel,220)
        IF (ModTyp(im).EQ.2) WRITE(imodel,230)
        IF (ModTyp(im).EQ.3) WRITE(imodel,240)

        WRITE(imodel,250) NRes(im)
        WRITE(imodel,260) NPer(im)
        WRITE(imodel,100) (Period(im,ip),ip=1,NPer(im))
        WRITE(imodel,280) NSta(im)
        WRITE(imodel,101) (StaPos(im,is),is=1,NSta(im))

        nc2  = MOD(NSta(im),7)
        IF (nc2.EQ.0) THEN
          ncol = IDINT(NSta(im)/D7)
        ELSE
          ncol = IDINT(NSta(im)/D7)+1
        ENDIF

        DO ir = 1,NRes(im)
           IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN
             IF (ir.EQ.1) WRITE(imodel,320) 'applog'
             IF (ir.EQ.2) WRITE(imodel,330) 'applog'
           ENDIF
           IF ((ResTyp(im,ir).EQ.2).OR.(ResTyp(im,ir).EQ.4)) THEN
             IF (ir.EQ.1) WRITE(imodel,320) 'phsrad'
             IF (ir.EQ.2) WRITE(imodel,330) 'phsrad'
           ENDIF
           IF (ResTyp(im,ir).EQ.5) THEN
             IF (ir.EQ.1) WRITE(imodel,320) 'rel   '
             IF (ir.EQ.2) WRITE(imodel,330) 'rel   '
           ENDIF
           IF (ResTyp(im,ir).EQ.6) THEN
             IF (ir.EQ.1) WRITE(imodel,320) 'img   '
             IF (ir.EQ.2) WRITE(imodel,330) 'img   '
           ENDIF
    

           DO ip = 1,NPer(im)
             iss = 0
             DO ic = 1,ncol
               isa = iss+1
               isb = iss+7
               IF (isb.GT.NSta(im)) isb = NSta(im)

               IF (ic.EQ.1) THEN
                 WRITE(imodel,105) Period(im,ip),
     >                  (EE(im,ir,ip,is),is=isa,isb)
               ELSE
                 WRITE(imodel,110) ' ',(EE(im,ir,ip,is),is=isa,isb)
               ENDIF
               iss = iss+7
             ENDDO ! ic
           ENDDO ! ip
        ENDDO ! ir
        CLOSE(imodel)
      ENDDO ! im
 
100   FORMAT(7F12.4)
101   FORMAT(7F12.2)
105   FORMAT(8F12.4)
110   FORMAT(a12,7F12.4)

180   FORMAT('TITLE  ABSOLUTE_ERROR')
220   FORMAT('MODE_TYPE               tm')
230   FORMAT('MODE_TYPE               te')
240   FORMAT('MODE_TYPE               tp')

250   FORMAT('NUMBER_OF_RESPONSE       ',i1)
260   FORMAT('NUMBER_OF_PERIOD       ',i3)
280   FORMAT('NUMBER_OF_STATION      ',i3)

320   FORMAT('ERROR_RESPONSE_NO_1    ',a6)
330   FORMAT('ERROR_RESPONSE_NO_2    ',a6)
 
      RETURN
      END ! WriteErrorRespond()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C


      SUBROUTINE WriteIntro(iu,NMode,ModTyp,NRes,NPer,NSta,
     >           ResTyp,NNT,NN,Ny,Nzb,
     >           REQUIRED_RMS,SD_RMS,MAX_ITERATION,MAX_SMOOTHING,
     >           MIN_LM,STARTING_LM,MAX_LM,STEPSIZE_LM,SMOOTH_SZLM,
     >           CONT_HIGHER_RMS,CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,
     >           FIX_LM,FWD_Only,BackGround_Rho,
     >           DTIME_STEP,DLENGTH_HOR,DLENGTH_VER,
     >           CPriorModel,DFStatus)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'

      INTEGER iu,NMode,ModTyp(*),NRes(*),NPer(*),NSta(*),Ny,Nzb
      INTEGER NNT(3),NN(NMODMX,NRESMX,2)
      INTEGER ResTyp(NMODMX,NRESMX)
      INTEGER DTIME_STEP,DFStatus(*)
      INTEGER MAX_ITERATION, MAX_SMOOTHING
      INTEGER CONT_HIGHER_RMS,CONT_NOTFOUND_RMS,CONT_HIGHER_MNORM,
     >        FIX_LM,FWD_Only
      REAL*8  REQUIRED_RMS,SD_RMS,MIN_LM,STARTING_LM,MAX_LM,
     >        STEPSIZE_LM,SMOOTH_SZLM
      REAL*8  DLENGTH_HOR,DLENGTH_VER,BackGround_Rho
      CHARACTER*80 CPriorModel

      INTEGER im,ir
      REAL*8  ratio,n3,n2


      WRITE(iu,1000)
      WRITE(iu,999)
      WRITE(iu,1200)
      WRITE(iu,1210) Ny*Nzb
      WRITE(iu,1220) Ny
      WRITE(iu,1230) Nzb
      WRITE(iu,1000)

      WRITE(iu,999)
      WRITE(iu,1100)
      WRITE(iu,1150) NMode
      DO im = 1,NMode
        IF (ModTyp(im).EQ.1) WRITE(iu,1154) NRes(im)
        IF (ModTyp(im).EQ.2) WRITE(iu,1155) NRes(im)
        IF (ModTyp(im).EQ.3) WRITE(iu,1156) NRes(im)
        DO ir = 1,NRes(im)
          IF (ResTyp(im,ir).EQ.1) WRITE(iu,1160) ir
          IF (ResTyp(im,ir).EQ.3) WRITE(iu,1161) ir
          IF (ResTyp(im,ir).EQ.2) WRITE(iu,1162) ir
          IF (ResTyp(im,ir).EQ.4) WRITE(iu,1163) ir
          IF (ResTyp(im,ir).EQ.5) WRITE(iu,1164) ir
          IF (ResTyp(im,ir).EQ.6) WRITE(iu,1165) ir
        ENDDO
        WRITE(iu,1170) NPer(im)
        WRITE(iu,1175) NSta(im)
      ENDDO

      WRITE(iu,999)
      WRITE(iu,1110) NNT(1)
      WRITE(iu,1120) NNT(2)
      DO im = 1,NMode
        WRITE(iu,1125) im,(NN(im,ir,1),ir=1,NRes(im))
      ENDDO
      n3 = NNT(3)
      n2 = NNT(2)
      ratio = n3/n2
      WRITE(iu,1140) NNT(3),ratio
      DO im = 1,NMode
        WRITE(iu,1127) im,(NN(im,ir,2),ir=1,NRes(im))
      ENDDO
      WRITE(iu,1000)
 


      WRITE(iu,999)
      WRITE(iu,1600)
      IF (FWD_Only.EQ.1) THEN
        WRITE(iu,1660)
        WRITE(iu,1000)
      ELSE
        WRITE(iu,1610) REQUIRED_RMS-SD_RMS,REQUIRED_RMS+SD_RMS
        WRITE(iu,1615) MAX_ITERATION
        WRITE(iu,1620) MAX_SMOOTHING
c       IF (CONT_HIGHER_RMS.EQ.1) WRITE(iu,1630)
c       IF (CONT_HIGHER_RMS.EQ.0) WRITE(iu,1635)
c       IF (CONT_NOTFOUND_RMS.EQ.1) WRITE(iu,1640)
c       IF (CONT_NOTFOUND_RMS.EQ.0) WRITE(iu,1645)
c       IF (CONT_HIGHER_MNORM.EQ.1) WRITE(iu,1650)
c       IF (CONT_HIGHER_MNORM.EQ.0) WRITE(iu,1655)
        WRITE(iu,1000)

c       WRITE(iu,999)
c       WRITE(iu,1700)
c       WRITE(iu,1710) MIN_LM
c       WRITE(iu,1715) STARTING_LM
c       WRITE(iu,1720) MAX_LM
c       WRITE(iu,1730) STEPSIZE_LM
c       WRITE(iu,1740) SMOOTH_SZLM
c       IF (FIX_LM.EQ.0) WRITE(iu,1760)
c       IF (FIX_LM.EQ.1) WRITE(iu,1770)
c       WRITE(iu,1000)

c       WRITE(iu,999)
c       WRITE(iu,1300)
c       WRITE(iu,1320) DTIME_STEP
c       IF (DLENGTH_HOR.EQ.D0) WRITE(iu,1330)
c       IF (DLENGTH_HOR.NE.D0) WRITE(iu,1335) DLENGTH_HOR
c       IF (DLENGTH_VER.EQ.D0) WRITE(iu,1340)
c       IF (DLENGTH_VER.NE.D0) WRITE(iu,1345) DLENGTH_VER
c       IF (CPriorModel.NE.'default') THEN
c         WRITE(iu,1350)
c       ELSE 
c         WRITE(iu,1355)
c       ENDIF
c       DO jj = 1,Ny*Nzb
c         IF (DFStatus(jj).NE.1) THEN
c           WRITE(iu,1360)
c           IF (DTIME_STEP.EQ.1) THEN 
c              WRITE(iu,1365)
c              DTIME_STEP = 10
c           ENDIF
c           GOTO 700
c         ENDIF
c       ENDDO
c700     CONTINUE
c       WRITE(iu,1370) BackGround_Rho
c       WRITE(iu,1000)

      ENDIF

999   FORMAT(' ')
1000  FORMAT('------------------------------------------------------')

1100  FORMAT('          <<<<< DATA PARAMETERS  >>>>>')
1110  FORMAT('TOTAL NO. OF DATA                            : ',i6)
1120  FORMAT('TOTAL NO. OF DATA USED TO ESTIMATE RMS   <N> : ',i6)
1125  FORMAT('      NO. OF DATA (EST. RMS) FOR MODE # ',i4,' : ',2i6)
1127  FORMAT('      NO. OF DATA (CAL. REP) FOR MODE # ',i4,' : ',2i6)
1140  FORMAT('TOTAL NO. OF DATA USED TO COMPUTE REP.   <L> : ',i6,
     >       ' (',f4.2,'N)')

1150  FORMAT('TOTAL NO. OF MODE                            : ',i6)
1154  FORMAT(' <TM>     NO. OF RESPONSES                   : ',i6)
1155  FORMAT(' <TE>     NO. OF RESPONSES                   : ',i6)
1156  FORMAT(' <TP>     NO. OF RESPONSES                   : ',i6)

1160  FORMAT('                    RESPONSE # <',i4,
     >       '>        :     APP')
1161  FORMAT('                    RESPONSE # <',i4,
     >       '>        :     LOG10 APP')
1162  FORMAT('                    RESPONSE # <',i4,
     >       '>        :     PHS (DEG)')
1163  FORMAT('                    RESPONSE # <',i4,
     >       '>        :     PHS (RAD)')
1164  FORMAT('                    RESPONSE # <',i4,
     >       '>        :     REAL')
1165  FORMAT('                    RESPONSE # <',i4,
     >       '>        :     IMAG')

1170  FORMAT('          NO. OF PERIODS                     : ',i6)
1175  FORMAT('          NO. OF STATIONS                    : ',i6)

1200  FORMAT('          <<<<< MODEL PARAMETERS >>>>>')
1210  FORMAT('TOTAL NO. OF MODEL BLOCKS                <M> : ',i6)
1220  FORMAT('      NO. OF HORIZONTAL BLOCK (NY)           : ',i6)
1230  FORMAT('      NO. OF VERTICAL BLOCK (NZ) (NO AIR)    : ',i6)

1300  FORMAT('  <<<<< MODEL CONTROL AND MODEL COVARIANCE >>>>')
1320  FORMAT('           PSEUDO-DIFFUSION TIME STEP  = ',i6)
1330  FORMAT('  HORIZONTAL LENGTH SCALES VARY WITH',
     >        ' <STATION SPACES AND DEPTH>') 
1335  FORMAT('           HORIZONTAL LENGTH SCALES = ',e12.4)
1340  FORMAT('        VERTICAL LENGTH SCALES VARY WITH <DEPTH>')
1345  FORMAT('           VERTICAL LENGTH SCALES = ',e12.4)
1350  FORMAT('       PRIOR MODEL IS INCLUDED IN THE INVERSION')
1355  FORMAT('     NO PRIOR MODEL IS INCLUDED IN THE INVERSION')
1360  FORMAT('    MODEL CONSTRAINT IS ASSIGNED IN THE INERSION')
1365  FORMAT(' ???? FOR BEST RESULT, DTIME_STEP WILL',
     >       ' BE INCREASED TO 10')
1370  FORMAT(' BACKGROUND RESIS. VALUE',
     >       ' (USING IN INTERP.)  = ',f12.4,' OHM-M')
 
1600  FORMAT('          <<<<<  GOAL OF INVERSION >>>>>')
1610  FORMAT('      DESIRED RMS IN BETWEEN ',f5.2,' AND ',f5.2)
1615  FORMAT('             MAX NO. OF ITERATION = ',i4)
1620  FORMAT('       MAX NO. OF SMOOTHING ITERATION = ',i4)
1630  FORMAT(' <CONTINUE> EVEN FOUND HIGHER RMS THAN PREV ITER ',
     >       '(PHASE I)')
1635  FORMAT('   <STOP> IF FOUND HIGHER RMS THAN PREV ITER (PHASE I)')
1640  FORMAT(' <CONTINUE> EVEN FOUND HIGHER RMS THAN PREV ITER ',
     >       '(PHASE II)')
1645  FORMAT('   <STOP> IF FOUND HIGHER RMS THAN PREV ITER (PHASE II)')
1650  FORMAT(' <CONTINUE> EVEN FOUND HIGHER MODEL NORM THAN ',
     >       'PREV ITER (PHASE II)')
1655  FORMAT(' <STOP> IF FOUND HIGHER MODEL NORM THAN ',
     >       'PREV ITER (PHASE II)')
1660  FORMAT(' ONLY <FORWARD MODELING> OF THE STARTING MODEL ',
     >       ' WILL BE SOLVED.')

1700  FORMAT('<<<<< LOG10 LAGRANGE MULTIPLIER (LGM) INFORMATION >>>>>')
1710  FORMAT('                  MIN          = ',f5.2)
1715  FORMAT('                STARTING       = ',f5.2)
1720  FORMAT('                  MAX          = ',f5.2)
1730  FORMAT('           STEPSIZE (PHASE I)  = ',f5.2)
1740  FORMAT('           STEPSIZE (PHASE II) = ',f5.2)
1760  FORMAT('   <TRIAL AND ERROR SEARCH> OF LGM AFTER EACH ITERATION')
1770  FORMAT('             <FIX> LGM (UNLESS NO PROGRESS)')




      RETURN
      END! WriteIntro
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE WriteRMS(iter,COut,num, 
     >           NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >           RMS,RMSS,RMSP)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*),
     >        ResTyp(NMODMX,NRESMX),iter
      REAL*8  RMS(NMODMX*NRESMX+1),RMSS(NMODMX,NRESMX,NSTAMX)
      REAL*8  RMSP(NMODMX,NRESMX,NPERMX)
      CHARACTER*3  num
      CHARACTER*80 COut,cout2
 
      INTEGER im,ir,ip,is,len,imodel,imr
      CHARACTER*2  cmd
      CHARACTER*22 cresp

      imodel = 31
      CALL Lenb(COut,len)
      cout2 = COut(1:len)//'_rms.'//num
      OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')

      imr = 1
      WRITE(imodel,300) iter,RMS(1)
      DO im = 1,NMode
        IF (ModTyp(im).EQ.1) cmd = 'TM'
        IF (ModTyp(im).EQ.2) cmd = 'TE'
        IF (ModTyp(im).EQ.3) cmd = 'TP'

        DO ir = 1,NRes(im)
           imr = imr + 1
           IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN
             cresp = 'Apparent Resistivity :'
           ENDIF
           IF ((ResTyp(im,ir).EQ.2).OR.(ResTyp(im,ir).EQ.4)) THEN
             cresp = 'Phase                :'
           ENDIF
           IF (ResTyp(im,ir).EQ.5) THEN
             cresp = 'Real Part            :'
           ENDIF
           IF (ResTyp(im,ir).EQ.6) THEN
             cresp = 'Imaginary Part       :'
           ENDIF
           WRITE(imodel,320) cmd,cresp,RMS(imr)
           WRITE(imodel,330) cmd,cresp
           WRITE(imodel,100) (RMSS(im,ir,is),is=1,NSta(im))

           WRITE(imodel,340) cmd,cresp
           WRITE(imodel,100) (RMSP(im,ir,ip),ip=1,NPer(im))
        ENDDO !ir
      ENDDO !im
      CLOSE(imodel)

100   FORMAT(7F12.4)

300   FORMAT('ITERATION# ',i3,' OVERALL RMS = ',F12.4)
320   FORMAT(a3,' ',a22,' RMS = ',F12.4)
330   FORMAT(a3,' ',a22,' RMS BY SITES ')
340   FORMAT(a3,' ',a22,' RMS BY PERIODS')

      RETURN
      END ! WriteRMS()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C

      SUBROUTINE WriteDistortion(iter,COut,num, 
     >           NMode,NRes,NPer,NSta,ModTyp,ResTyp,
     >           SSIndx,SSPara,StsInx)
      INCLUDE 'parameter.h'
      INCLUDE 'constant.h'
 
      INTEGER NMode,NRes(*),NPer(*),NSta(*),ModTyp(*),
     >        ResTyp(NMODMX,NRESMX),iter
      INTEGER SSIndx(NMODMX,NSTAMX)
      INTEGER StsInx(NMODMX,NPERMX,NSTAMX)
      REAL*8  SSPara(NMODMX,NSTAMX)
      CHARACTER*3  num
      CHARACTER*80 COut,cout2
 
      INTEGER im,ir,ip,is,len,imodel
      CHARACTER*256 cindex

      imodel = 31
      CALL Lenb(COut,len)

      DO im = 1,NMode
        IF (ModTyp(im).EQ.1) cout2 = COut(1:len)//'_dist_tm.'//num
        IF (ModTyp(im).EQ.2) cout2 = COut(1:len)//'_dist_te.'//num
        IF (ModTyp(im).EQ.3) GOTO 10
        OPEN(UNIT=imodel,FILE=cout2,STATUS='unknown')

        DO ir = 1,NRes(im)
           IF ((ResTyp(im,ir).EQ.1).OR.(ResTyp(im,ir).EQ.3)) THEN
             WRITE(imodel,300)
             DO is = 1,NSta(im)
              WRITE(cindex(is:is),'(i1)') SSIndx(im,is)
             ENDDO
             WRITE(imodel,*) cindex(1:NSta(im))

             WRITE(imodel,310)
             WRITE(imodel,100) (SSPara(im,is),is=1,NSta(im))

             WRITE(imodel,320)
             DO ip = 1,NPer(im)
               DO is = 1,NSta(im)
                WRITE(cindex(is:is),'(i1)') StsInx(im,ip,is)
               ENDDO
               WRITE(imodel,*) cindex(1:NSta(im))
             ENDDO

           ENDIF
        ENDDO !ir
10      CONTINUE
        CLOSE(imodel)
      ENDDO !im

100   FORMAT(7F12.4)

300   FORMAT('DISTORTION_INDEX')
310   FORMAT('DISTORTION_PARAMETER')
320   FORMAT('DISTORTION_INCLUSION  index')

      RETURN
      END ! WriteDistortion()

C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
