      SUBROUTINE D02LAF(F,NEQ,T,TEND,Y,YP,YDP,RWORK,LRWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 16A REVISED. IER-972 (JUN 1993).
C
C     NOTE -
C     THE ARRAYS IN COMMON ARE LARGE ENOUGH TO ACCOMODATE THE
C     COEFFICIENTS OF THE 12(10) FORMULA PAIR.  THIS ADDS CLARITY
C     TO THE PROGRAM BUT SOME STORAGE COULD BE SAVED IF NECESSARY
C     (BY ADDITIONAL PROGRAMMING EFFORT).  TWO POSSIBLE WAYS OF
C     SAVING STORAGE ARE:
C     1) STORE THE MATRIX OF COEFFICIENTS  A  IN STRICTLY LOWER
C        TRIANGULAR FORM.  THIS WOULD SAVE  153  LOCATIONS.
C     2) INCLUDE THE STORAGE REQUIREMENT FOR THE COEFFICIENTS
C        IN THE WORKSPACE ARRAY  RWORK.
C
C
C     IF THE USE OF COMMON STATEMENTS IS TO BE AVOIDED, /AD02LA/ CAN
C     BE REPLACED BY CALLING THE APPROPRIATE COEFFICIENT ROUTINE
C     (D02LAY OR D02LAZ) ON EACH ENTRY TO D02LAF.  SIMILARLY, THE
C     VALUES IN /CD02LA/ COULD BE RECOMPUTED ON EACH ENTRY.
C     OTHER MEANS WILL HAVE TO BE FOUND TO REPLACE /BD02LA/.
C
C     .. Parameters ..
      INTEGER           ASK, OKAY, SET
      PARAMETER         (ASK=0,OKAY=1,SET=1)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02LAF')
      DOUBLE PRECISION  R, RMN1, RMN2, RMN3, SCALEL, SCALEH
      PARAMETER         (R=10.0D0,RMN1=1.0D0/R,RMN2=1.0D0/(R*R),
     *                  RMN3=1.0D0/(R*R*R),SCALEL=0.8D0,SCALEH=0.75D0)
      INTEGER           NOVHD, HUNDRD
      PARAMETER         (NOVHD=16,HUNDRD=100)
      DOUBLE PRECISION  ONE, ZERO, FIFTH, TWO, HALF, QUAR, ELVNTH, FOUR,
     *                  EIGHT, TEN, TWENTY, TWOHUN, TENPC
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0,FIFTH=0.2D0,TWO=2.0D0,
     *                  HALF=0.5D0,QUAR=0.25D0,ELVNTH=1.0D0/11.0D0,
     *                  FOUR=4.0D0,EIGHT=8.0D0,TEN=1.0D1,TWENTY=2.0D1,
     *                  TWOHUN=2.0D2,TENPC=0.1D0)
      INTEGER           SVINIH, SVHUSD, SVHNXT, SVOKST, SVFLST, SVMXST,
     *                  SVSTRT, SV1STP, SVOPTS, SVMETH, SVATST, SVTOL,
     *                  SVEPS, SVTOLD, SVNEQ, SVLRWK
      PARAMETER         (SVINIH=1,SVHUSD=2,SVHNXT=3,SVOKST=4,SVFLST=5,
     *                  SVMXST=6,SVSTRT=7,SV1STP=8,SVOPTS=9,SVMETH=10,
     *                  SVATST=11,SVTOL=12,SVEPS=13,SVTOLD=14,SVNEQ=15,
     *                  SVLRWK=16)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TEND
      INTEGER           IFAIL, LRWORK, NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), Y(NEQ), YDP(NEQ), YP(NEQ)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Scalars in Common ..
      DOUBLE PRECISION  DIR, EXPON, H, RANGE, RHIGH, SCALE, STBRAD,
     *                  TCHECK, TOL, TOOSML
      INTEGER           EFFCT1, EFFCT2, FINPH2, ISTP, NASTEP, NATT,
     *                  NBASE, NFLSTP, NRT1, NRT2, NRT3, NRT4, NSTAGE,
     *                  NSTMAX, NTHR, NTHRP, NWT, NWTP
      LOGICAL           FAILED, FIRST, HGIVEN, LAST, OPTRED, PHASE2,
     *                  SUCCES
C     .. Arrays in Common ..
      DOUBLE PRECISION  A(17,17), B(16), BHAT(15), BP(16), BPHAT(17),
     *                  C(17)
      INTEGER           PTR(16)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSH, ABSHP, ALPHA, BETA, DEL, ERR, ERRP, FDEL,
     *                  HSQ, HUSED, SCALE1, SCL, SCLP, SK, SPK, TAU,
     *                  THRES, THRESP, TOLD, TRY, WT, WTP, YDPDEL,
     *                  YDPNRM, YINTK, YPDEL, YPINTK, YPNRM
      INTEGER           ERSTAT, I, IER, J, JSTAGE, K, NREC, STATE
      LOGICAL           DFLT, DFLTP, INEFFC, LOOP, WASTE
C     .. Local Arrays ..
      CHARACTER*80      REC(5)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02LAX, D02LAY, D02LAZ, D02LXZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN, SIGN
C     .. Common blocks ..
      COMMON            /AD02LA/A, C, BPHAT, B, BP, BHAT, PTR
      COMMON            /BD02LA/DIR, RANGE, TCHECK, H, EXPON, RHIGH,
     *                  TOL, SCALE, STBRAD, TOOSML, NSTAGE, NSTMAX,
     *                  NASTEP, ISTP, NFLSTP, NATT, EFFCT1, EFFCT2,
     *                  FINPH2, LAST, FAILED, FIRST, HGIVEN, SUCCES,
     *                  OPTRED, PHASE2
      COMMON            /CD02LA/NRT1, NRT2, NRT3, NRT4, NBASE, NTHR,
     *                  NTHRP, NWT, NWTP
C     .. Save statement ..
      SAVE              /AD02LA/, /BD02LA/, /CD02LA/, LOOP
C     .. Data statements ..
      DATA              LOOP/.FALSE./
C     .. Executable Statements ..
C
      CALL D02LXZ(STATE,ASK)
      IF (STATE.NE.OKAY) THEN
         IER = 1
         IF (STATE.EQ.0) THEN
            WRITE (REC(1),FMT=99993)
            NREC = 1
         ELSE
            WRITE (REC(1),FMT=99992)
            WRITE (REC(2),FMT=99991) STATE - 1
            NREC = 2
         END IF
         GO TO 620
      END IF
C
      IER = 0
      NREC = 0
      NASTEP = 0
C
      IF (NEQ.NE.INT(RWORK(SVNEQ))) THEN
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99989)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99988) INT(RWORK(SVNEQ)), NEQ
         IER = 1
      END IF
      IF (LRWORK.NE.INT(RWORK(SVLRWK))) THEN
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99987)
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99986) INT(RWORK(SVLRWK)), LRWORK
         IER = 1
      END IF
      IF (T.EQ.TEND) THEN
         NREC = NREC + 1
         WRITE (REC(NREC),FMT=99998) T
         NREC = 1
         IER = 1
      END IF
      IF (IER.EQ.1) GO TO 620
C
      IF (RWORK(SVSTRT).EQ.ONE) THEN
C
C        BLOCK A.
C        FIRST CALL IN AN INTEGRATION.
C        EXTRACT DATA FROM SETUP ROUTINE D02LXF.
C
         TOLD = T
         H = RWORK(SVINIH)
         HUSED = ZERO
         RWORK(SVINIH) = ZERO
         NSTMAX = INT(RWORK(SVMXST)+FIFTH)
         RWORK(SVOPTS) = ZERO
         RWORK(SVSTRT) = ZERO
         TOL = RWORK(SVTOL)
C
C        INITIALIZATION.
C
         EFFCT1 = 0
         EFFCT2 = 0
         ISTP = 0
         NFLSTP = 0
         NATT = 0
         FIRST = .TRUE.
         LOOP = .FALSE.
         RANGE = ABS(TEND-T)
         DIR = SIGN(ONE,TEND-T)
C
C        THE NEXT THREE LINES RELATE TO THE CODE CHOICE OF H0.
C
         OPTRED = .FALSE.
         HGIVEN = H .NE. ZERO
         PHASE2 = FIRST .AND. .NOT. HGIVEN
C
C        SET UP POINTERS TO RWORK.
C
         NRT1 = NOVHD
         NRT2 = NOVHD + NEQ
         NRT3 = NOVHD + 2*NEQ
         NRT4 = NOVHD + 3*NEQ
         NTHR = NOVHD + 4*NEQ
         NTHRP = NOVHD + 5*NEQ
         NBASE = NOVHD + 6*NEQ
C
C        PICK UP CONSTANTS DEPENDING ON ORDER OF METHOD SELECTED AND SET
C        POINTERS. THE LOCAL ARRAY PTR CONTAINS POINTERS TO WHERE THE
C        STAGES ARE HELD IN THE ARRAY RWORK. THIS IS TO ALLOW SAVINGS IN
C        STORAGE FOR 12(10) PAIR AND BOTH PAIRS TO BE IMPLEMENTED IN THE
C        SAME BLOCK OF CODE.
C
         RHIGH = RWORK(SVMETH)
         IF (RWORK(SVMETH).EQ.ZERO) THEN
            EXPON = FIFTH
            STBRAD = FOUR
            TOOSML = TWENTY*RWORK(SVEPS)
            NSTAGE = 6
            SCALE = SCALEL
            CALL D02LAZ(A,C,BPHAT,B,BP)
            FINPH2 = 4
            NWT = NBASE + 3*NEQ
            NWTP = NBASE + 4*NEQ
            PTR(1) = 0
            PTR(2) = 1
            PTR(3) = 2
            PTR(4) = 3
            PTR(5) = 4
         ELSE
            EXPON = ELVNTH
            STBRAD = EIGHT
            TOOSML = TWOHUN*RWORK(SVEPS)
            NSTAGE = 17
            SCALE = SCALEH
            CALL D02LAY(A,C,BHAT,BPHAT,B,BP)
            FINPH2 = 14
            NWT = NBASE + 12*NEQ
            NWTP = NBASE + 13*NEQ
            PTR(1) = 0
            PTR(2) = 1
            PTR(3) = 2
            PTR(4) = 0
            PTR(5) = 3
            PTR(6) = 4
            PTR(7) = 5
            PTR(8) = 6
            PTR(9) = 7
            PTR(10) = 8
            PTR(11) = 9
            PTR(12) = 10
            PTR(13) = 11
            PTR(14) = 12
            PTR(15) = 13
            PTR(16) = 3
         END IF
C
         CALL F(NEQ,T,Y,YDP)
C
         IF ( .NOT. HGIVEN) THEN
C
C           PHASE 1 OF CALCULATION OF AN INITIAL STEPSIZE H0 BY THE CODE
C           IF || YP(T0) || .NE. 0 THEN H0 = (TOL**EXPON) / ||YP(T0)||
C                          ELSE H0 = TEND - T0
C           IF || YDP(T0) || .NE. 0 THEN
C              HP0 = (TOL**EXPON) / ||YDP(T0)||
C           ELSE HP0 = TEND - T0
C           H = MIN(H0, HP0)
C           NOTE: 1) EXPON = 1/5  FOR THE 6(4) PAIR,
C           EXPON = 1/11 FOR THE 12(10) PAIR.
C           2) THE STRATEGY GIVEN ABOVE HAS BEEN MODIFIED TO TRY TO
C           AVOID SELECTING TOO SMALL A VALUE FOR H0 WHEN THE
C           INITIAL VALUES OF Y AND YP ARE SMALL.
C
            YPNRM = ZERO
            YDPNRM = ZERO
            DFLT = .FALSE.
            DFLTP = .FALSE.
C
C           THE FOLLOWING LOOP WAS NOT VECTORIZABLE ON A  CRAY 1S.
C           THE FOLLOWING LOOP WAS NOT VECTORIZABLE ON A CDC
C           CYBER 205.
C
C           EVALUATE NORMS AS DESCRIBED ABOVE AND INITIALISE THE
C           WEIGHT VECTORS FOR USE IN PHASE2.
C
            DO 20 K = 1, NEQ
               THRES = RWORK(NTHR+K)
               WT = MAX(THRES,ABS(Y(K)))
               RWORK(NWT+K) = WT
               TRY = ABS(YP(K)/WT)
               IF (TRY.GT.YPNRM) THEN
                  YPNRM = TRY
                  DFLT = WT .EQ. THRES
               END IF
               THRESP = RWORK(NTHRP+K)
               WTP = MAX(THRESP,ABS(YP(K)))
               RWORK(NWTP+K) = WTP
               TRY = ABS(YDP(K)/WTP)
               IF (TRY.GT.YDPNRM) THEN
                  YDPNRM = TRY
                  DFLTP = WTP .EQ. THRESP
               END IF
   20       CONTINUE
            ABSH = RANGE
            ABSHP = RANGE
            TAU = TOL**EXPON
            IF (ABSH*YPNRM.GT.TAU .AND. .NOT. DFLT) ABSH = TAU/YPNRM
            IF (ABSHP*YDPNRM.GT.TAU .AND. .NOT. DFLTP)
     *          ABSHP = TAU/YDPNRM
            H = MIN(ABSH,ABSHP)
C
C           END OF PHASE 1.
C
         END IF
C
C        END OF BLOCK A.
C
      ELSE
C
C        BLOCK B.
C        CONTINUATION CALL.
C        CHECK DATA FOR CONSISTENCY.
C
         TOLD = RWORK(SVTOLD)
         HUSED = RWORK(SVHUSD)
         IF (DIR*(TEND-T).LT.ZERO) THEN
            WRITE (REC(1),FMT=99997)
            NREC = 1
            IER = 1
            GO TO 620
         END IF
C
C        CHECK FOR OPTIONAL INPUTS.
C
         IF (RWORK(SVOPTS).EQ.ONE) THEN
            H = RWORK(SVHNXT)
            NSTMAX = INT(RWORK(SVMXST)+FIFTH)
            RWORK(SVOPTS) = ZERO
            TOL = RWORK(SVTOL)
            IF (RHIGH.NE.RWORK(SVMETH)) THEN
               RHIGH = RWORK(SVMETH)
               IF (RWORK(SVMETH).EQ.ZERO) THEN
                  EXPON = FIFTH
                  STBRAD = FOUR
                  TOOSML = TWENTY*RWORK(SVEPS)
                  NSTAGE = 6
                  SCALE = SCALEL
                  FINPH2 = 4
                  NWT = NBASE + 3*NEQ
                  NWTP = NBASE + 4*NEQ
                  PTR(1) = 0
                  PTR(2) = 1
                  PTR(3) = 2
                  PTR(4) = 3
                  PTR(5) = 4
                  CALL D02LAZ(A,C,BPHAT,B,BP)
               ELSE
                  EXPON = ELVNTH
                  STBRAD = EIGHT
                  TOOSML = TWOHUN*RWORK(SVEPS)
                  NSTAGE = 17
                  SCALE = SCALEH
                  FINPH2 = 14
                  NWT = NBASE + 12*NEQ
                  NWTP = NBASE + 13*NEQ
                  PTR(1) = 0
                  PTR(2) = 1
                  PTR(3) = 2
                  PTR(4) = 0
                  PTR(5) = 3
                  PTR(6) = 4
                  PTR(7) = 5
                  PTR(8) = 6
                  PTR(9) = 7
                  PTR(10) = 8
                  PTR(11) = 9
                  PTR(12) = 10
                  PTR(13) = 11
                  PTR(14) = 12
                  PTR(15) = 13
                  PTR(16) = 3
                  CALL D02LAY(A,C,BHAT,BPHAT,B,BP)
               END IF
            END IF
         END IF
C
C        END OF BLOCK B.
C
      END IF
C
C     BLOCK C.
C     THE START OF A NEW STEP.
C
      H = SIGN(ABS(H),DIR)
C
   40 CONTINUE
      SUCCES = .FALSE.
   60 LAST = .FALSE.
C
C     THE CODE IS REQUIRED TO HIT TEND EXACTLY WITH A REASONABLE
C     STEPSIZE. CHECK THAT  T+H .LE. TEND  OR  T+2*H .LE. TEND
C     IN THE DIRECTION OF INTEGRATION. A CHECK IS INCORPORATED
C     LATER IN THE CODE TO CHECK WHEN H IS REDUCED TOO OFTEN BY
C     A SIGNIFICANT AMOUNT, HENCE IMPAIRING EFFICIENCY.
C
      IF (DIR*(T+H-TEND).GE.ZERO) THEN
         INEFFC = (TEND-T)/H .LT. HALF
         H = TEND - T
         LAST = .TRUE.
      ELSE
         IF (DIR*(T+TWO*H-TEND).GE.ZERO) THEN
            H = (TEND-T)/TWO
         END IF
      END IF
C
      FAILED = .FALSE.
C
C     RESTART HERE AFTER A STEP FAILURE.
C
   80 CONTINUE
C
C     THE CHECK FOR LIMITING PRECISION, USING  TOOSML,  IS BASED ON
C     THE DISTANCE BETWEEN COEFFICIENTS  C(I)  IN THE RKN FORMULAS,
C     AS DESCRIBED IN
C     SHAMPINE,L.F.,WATTS,H.A. (1979).  THE ART OF WRITING A
C     RUNGE-KUTTA CODE.  APPL. MATH. COMPUTATION 5, 93-121.
C
      IF (ABS(H).LT.TOOSML*MAX(ABS(T),ABS(T+H))) THEN
C
C        STEP TOO SMALL. CONTROL RETURNS TO USER THROUGH BOTTOM OF CODE.
C
         WRITE (REC(1),FMT=99996) H
         WRITE (REC(2),FMT=99995) T
         NREC = 2
         IER = 3
         GO TO 600
      END IF
      NASTEP = NASTEP + 1
      IF (NASTEP.GT.NSTMAX) THEN
         WRITE (REC(1),FMT=99994) NSTMAX
         NREC = 1
         IER = 2
         GO TO 600
      END IF
C
C     END OF BLOCK C.
C
C     BLOCK D.
C     FORM THE STAGES AND STORE THEM IN THE WORKING STORAGE ARRAY,
C     RWORK. TEMPORARY STORAGE:
C      RWORK(NRT1+K) CONTAINS THE SUM OF THE STAGES  SUM A(I,J)*G(J)
C          USED IN FORMING THE INTERNAL Y VALUES. AT THE END OF
C          THE BLOCK IT WILL CONTAIN   SUM BHAT(I) * G(I)   AND
C          WILL BE USED IN ERROR ESTIMATION FOR Y.
C      RWORK(NRT2+K) WILL CONTAIN  SUM BPHAT(I) * G(I)   AND WILL
C          BE USED IN ESTIMATING THE ERROR IN YP.
C     AT THE END OF THE BLOCK, THE NEW APPROXIMATION TO Y(K) WILL BE
C     IN RWORK(NRT3+K) AND THE APPROXIMATION TO YP(K) WILL BE IN
C     RWORK(NRT4+K).
C
C     NOTE.
C     MANY OF THE FOLLOWING LOOPS OVER  K=1,NEQ  HAVE CONSTANT ARRAY
C     VALUES INSIDE. THE CODE IS WRITTEN WITH CLARITY IN MIND. ANY
C     OPTIMIZING COMPILER WILL IDENTIFY THESE OCCURRENCES AND MOVE
C     VALUES OUTSIDE THE LOOPS. TWO LOOPS CONTAIN CHECKS FOR ZERO
C     MULTIPLIERS SO AS TO PREVENT NEEDLESS COMPUTATION.
C
      HSQ = H*H
      DO 100 K = 1, NEQ
         RWORK(NRT2+K) = BPHAT(1)*YDP(K)
  100 CONTINUE
      DO 260 I = 2, NSTAGE
C
         DO 120 K = 1, NEQ
            RWORK(NRT1+K) = A(I,1)*YDP(K)
  120    CONTINUE
C
         DO 160 J = 2, I - 1
            JSTAGE = NBASE + PTR(J-1)*NEQ
            IF (A(I,J).NE.ZERO) THEN
               DO 140 K = 1, NEQ
                  RWORK(NRT1+K) = RWORK(NRT1+K) + A(I,J)*RWORK(JSTAGE+K)
  140          CONTINUE
            END IF
  160    CONTINUE
C
         DO 180 K = 1, NEQ
            RWORK(NRT3+K) = Y(K) + C(I)*H*YP(K) + HSQ*RWORK(NRT1+K)
  180    CONTINUE
C
         JSTAGE = NBASE + PTR(I-1)*NEQ
         CALL F(NEQ,T+C(I)*H,RWORK(NRT3+1),RWORK(JSTAGE+1))
C
         IF (BPHAT(I).NE.ZERO) THEN
            DO 200 K = 1, NEQ
               RWORK(NRT2+K) = RWORK(NRT2+K) + BPHAT(I)*RWORK(JSTAGE+K)
  200       CONTINUE
         END IF
C
C        CONTROL USUALLY GOES FROM HERE TO THE START OF THE LOOP
C        (I=2,NSTAGE), UNLESS IT IS THE FIRST STEP AND THE CODE IS
C        COMPUTING H0, IN WHICH CASE THE NEXT BLOCK OF CODE IS EXECUTED
C
         IF (PHASE2) THEN
            YPNRM = ZERO
            YDPNRM = ZERO
            DEL = ZERO
            FDEL = ZERO
            YPDEL = ZERO
            YDPDEL = ZERO
C
C           THE FOLLOWING LOOP WAS NOT VECTORIZABLE ON A CDC
C           CYBER 205.
            DO 220 K = 1, NEQ
               YINTK = RWORK(NRT3+K)
               YPINTK = YP(K) + H*RWORK(NRT2+K)
               WT = MAX(RWORK(NWT+K),ABS(YINTK))
               WTP = MAX(RWORK(NWTP+K),ABS(YPINTK))
               RWORK(NWT+K) = WT
               RWORK(NWTP+K) = WTP
               YPNRM = MAX(YPNRM,ABS(YINTK)/WT,ABS(Y(K))/WT)
               YDPNRM = MAX(YDPNRM,ABS(YPINTK)/WTP,ABS(YP(K))/WTP)
               DEL = MAX(DEL,ABS(YINTK-Y(K))/WT)
               FDEL = MAX(FDEL,ABS(YPINTK-YP(K))/WT)
               YPDEL = MAX(YPDEL,ABS(YPINTK-YP(K))/WTP)
               YDPDEL = MAX(YDPDEL,ABS(RWORK(JSTAGE+K)-YDP(K))/WTP)
  220       CONTINUE
            DEL = MAX(DEL,ABS(C(I)*H)/RANGE)
            SCL = ONE
            SCLP = ONE
            IF (DEL.GT.TEN*RWORK(SVEPS)*YPNRM) THEN
               IF (ABS(H)*FDEL.GT.STBRAD*DEL)
     *             SCL = STBRAD/R*MAX(DEL/(ABS(H)*FDEL),RMN3)
            END IF
            IF (YPDEL.GT.TEN*RWORK(SVEPS)*YDPNRM) THEN
               IF (ABS(H)*YDPDEL.GT.STBRAD*YPDEL)
     *             SCLP = STBRAD/R*MAX(YPDEL/(ABS(H)*YDPDEL),RMN3)
            END IF
            SCALE1 = MIN(SCL,SCLP)
            IF (SCALE1.LT.ONE) THEN
               NATT = NATT + 1
               H = SCALE1*ABS(H)
               H = DIR*H
               LAST = .FALSE.
C
C              RESET THE WEIGHT VECTORS
C
               DO 240 K = 1, NEQ
                  RWORK(NWT+K) = MAX(ABS(Y(K)),RWORK(NTHR+K))
                  RWORK(NWTP+K) = MAX(ABS(YP(K)),RWORK(NTHRP+K))
  240          CONTINUE
               GO TO 80
            END IF
C
C           TO SAVE ON STORAGE THE WEIGHT VECTORS RWORK(NWT+...) AND
C           RWORK(NWTP+...) ARE IN LOCATIONS WHERE THE STAGES ARE
C           USUALLY KEPT AND HENCE PHASE2 OPERATES UPTO THE FOURTH
C           STAGE FOR THE 6(4) PAIR AND THE FOURTEENTH STAGE FOR THE
C           12(10) PAIR.
C
            PHASE2 = I .LT. FINPH2
         END IF
  260 CONTINUE
C
      DO 280 K = 1, NEQ
         RWORK(NRT4+K) = YP(K) + H*RWORK(NRT2+K)
  280 CONTINUE
C
C     THE 12(10) PAIR DOES NOT USE FSAL (FIRST STAGE ON STEP N+1 IS SAME
C     AS LAST STAGE ON STEP N), HENCE MUST ADD STAGES UP TO OBTAIN THE
C     APPROXIMATIONS. UTILIZE THE FACT THAT
C     BHAT(2)=BHAT(3)=BHAT(4)=BHAT(5)=BHAT(6)=BHAT(16)=BHAT(17)=0.
C
      IF (RHIGH.EQ.ONE) THEN
         DO 300 K = 1, NEQ
            RWORK(NRT1+K) = BHAT(1)*YDP(K)
  300    CONTINUE
         DO 340 J = 7, 15
            JSTAGE = NBASE + PTR(J-1)*NEQ
            DO 320 K = 1, NEQ
               RWORK(NRT1+K) = RWORK(NRT1+K) + BHAT(J)*RWORK(JSTAGE+K)
  320       CONTINUE
  340    CONTINUE
         DO 360 K = 1, NEQ
            RWORK(NRT3+K) = Y(K) + H*YP(K) + HSQ*RWORK(NRT1+K)
  360    CONTINUE
      END IF
C
C     END OF BLOCK D.
C
      PHASE2 = .FALSE.
C
C     BLOCK E.
C     ERROR ESTIMATION.
C     USE THE TEMPORARY STORAGE VECTORS RWORK(NRT1+K) AND RWORK(NRT2+K)
C     TO COMPUTE  YHAT - Y  AND  YPHAT - YP  RESPECTIVELY.
C     SEPARATE SECTIONS ARE USED FOR THE HIGH AND LOW ORDER FORMULAS
C     IN ORDER TO TAKE ADVANTAGE OF THE FACT THAT IN THE LOW ORDER
C     CASE, B(4) = B(5) = B(6) = ZERO, WHILE IN THE HIGH ORDER CASE
C     B(2) = B(3) = B(4) = B(5) = B(6) = B(17) = ZERO AND SIMILARLY FOR
C     BP.
C
      DO 380 K = 1, NEQ
         RWORK(NRT1+K) = RWORK(NRT1+K) - B(1)*YDP(K)
         RWORK(NRT2+K) = RWORK(NRT2+K) - BP(1)*YDP(K)
  380 CONTINUE
C
      IF (RHIGH.EQ.ZERO) THEN
C
C        ERROR ESTIMATES FOR THE 6(4) PAIR.
C
         DO 420 I = 2, 3
            JSTAGE = NBASE + PTR(I-1)*NEQ
            DO 400 K = 1, NEQ
               RWORK(NRT1+K) = RWORK(NRT1+K) - B(I)*RWORK(JSTAGE+K)
  400       CONTINUE
  420    CONTINUE
         DO 460 I = 2, NSTAGE
            JSTAGE = NBASE + PTR(I-1)*NEQ
            DO 440 K = 1, NEQ
               RWORK(NRT2+K) = RWORK(NRT2+K) - BP(I)*RWORK(JSTAGE+K)
  440       CONTINUE
  460    CONTINUE
C
      ELSE
C
C        ERROR ESTIMATES FOR 12(10) PAIR.
C
         DO 500 I = 7, 16
            JSTAGE = NBASE + PTR(I-1)*NEQ
            DO 480 K = 1, NEQ
               RWORK(NRT1+K) = RWORK(NRT1+K) - B(I)*RWORK(JSTAGE+K)
               RWORK(NRT2+K) = RWORK(NRT2+K) - BP(I)*RWORK(JSTAGE+K)
  480       CONTINUE
  500    CONTINUE
      END IF
C
C     FORM ERR, THE MAXIMUM OF THE WEIGHTED MAX NORMS OF THE ESTIMATED
C     LOCAL ERRORS IN Y AND YP.
C
      ERR = ZERO
      ERRP = ZERO
      DO 520 K = 1, NEQ
         SK = HSQ*RWORK(NRT1+K)
         SPK = H*RWORK(NRT2+K)
         WT = HALF*ABS(Y(K)) + QUAR*ABS(RWORK(NRT3+K)) +
     *        QUAR*ABS(RWORK(NRT3+K)-SK)
         WT = MAX(WT,RWORK(NTHR+K))
         WTP = HALF*ABS(YP(K)) + QUAR*ABS(RWORK(NRT4+K)) +
     *         QUAR*ABS(RWORK(NRT4+K)-SPK)
         WTP = MAX(WTP,RWORK(NTHRP+K))
         ERR = MAX(ERR,ABS(SK)/WT)
         ERRP = MAX(ERRP,ABS(SPK)/WTP)
  520 CONTINUE
      ERR = MAX(ERR,ERRP)
C
C     END OF BLOCK E.
C
      IF (ERR.GT.TOL) THEN
C
C        BLOCK F.
C        FAILED STEP.
C
         LAST = .FALSE.
         IF (FIRST .AND. .NOT. HGIVEN) THEN
C
C           PHASE 3 FOR CODE CHOICE OF H0.
C
            NATT = NATT + 1
            IF (SUCCES) THEN
C
C              THE CODE HAS DISCARDED AN INITIAL STEP IN FAVOUR
C              OF A MUCH LARGER PREDICTED STEP BUT THIS LARGER
C              STEP HAS FAILED.  THEREFORE CARRY OUT OPTIMAL
C              REDUCTION.
C
               OPTRED = .TRUE.
               ALPHA = SCALE*(TOL/ERR)**EXPON
               ALPHA = MAX(RMN2,ALPHA)
            ELSE
C
C              NO SUCCESSFUL STEP YET.  REDUCE H0 TO H0/R AND
C              RESET WEIGHT VECTORS.
C
               ALPHA = RMN1
               PHASE2 = .TRUE.
               DO 540 K = 1, NEQ
                  RWORK(NWT+K) = MAX(ABS(Y(K)),RWORK(NTHR+K))
                  RWORK(NWTP+K) = MAX(ABS(YP(K)),RWORK(NTHRP+K))
  540          CONTINUE
            END IF
         ELSE
C
C           NOT THE FIRST STEP (OR H0 WAS SUPPLIED BY THE USER),
C           SO USE THE NORMAL STEP REDUCTION ALGORITHM.
C
            NFLSTP = NFLSTP + 1
            IF (FAILED) THEN
               ALPHA = HALF
            ELSE
               FAILED = .TRUE.
               ALPHA = SCALE*(TOL/ERR)**EXPON
               ALPHA = MAX(RMN1,ALPHA)
            END IF
         END IF
         H = ALPHA*H
C
C        TRY THE STEP AGAIN.
C
         GO TO 80
C
C        END OF BLOCK F.
C
      ELSE
C
C        BLOCK G.
C        SUCCESSFUL STEP.
C
         SUCCES = .TRUE.
         BETA = ((ERR/TOL)**EXPON)/SCALE
         IF (FIRST) THEN
C
C           PHASE 3 FOR CODE CHOICE OF H0.
C
            ALPHA = ONE/MAX(RMN3,BETA)
         ELSE
C
C           NORMAL STEPSIZE PREDICTION.
C
            ALPHA = ONE/MAX(RMN1,BETA)
         END IF
         IF (FAILED) ALPHA = MIN(ONE,ALPHA)
         HUSED = H
         H = ALPHA*H
C
         IF (FIRST .AND. .NOT. LAST .AND. .NOT. OPTRED .AND. ALPHA.GT.R)
     *       THEN
C
C           PHASE 3 FOR CODE CHOICE OF H0. THE CODE HAS ATTEMPTED THE
C           FIRST STEP, IS NOT AT THE END OF THE INTEGRATION RANGE, HAS
C           NOT ENCOUNTERED A STEP FAILURE IN PHASE3, AND THE PREDICTED
C           INCREASE IS LARGER THAN THE MAXIMUM PERMITTED ON A NORMAL
C           STEP. RETAKE THE STEP.
C
            HUSED = ZERO
            NATT = NATT + 1
            GO TO 60
         END IF
         IF (FIRST) RWORK(SVINIH) = HUSED
         FIRST = .FALSE.
         TOLD = T
         T = T + HUSED
         EFFCT1 = EFFCT1 + 1
         IF (LAST) THEN
            T = TEND
            IF (INEFFC) EFFCT2 = EFFCT2 + 1
C
C           TEST FOR INEFFICIENT USE OF FORMULAS FOR OUTPUT.
C
            IF (EFFCT1.GT.HUNDRD) THEN
               WASTE = DBLE(EFFCT2)/DBLE(EFFCT1) .GT. TENPC
               IF (WASTE) THEN
                  WRITE (REC(1),FMT=99990)
                  NREC = 1
                  IER = 5
                  EFFCT1 = 0
                  EFFCT2 = 0
               END IF
            END IF
         END IF
         ISTP = ISTP + 1
C
C        ON A SUCCESSFUL STEP, COPY THE NEW APPROXIMATIONS INTO
C        Y AND YP. SAVE THE OLD VALUES OF Y, YP AND YDP IN RWORK,
C        ALONG WITH THE INTERMEDIATE STAGES FOR USE IN CASE DENSE
C        OUTPUT IS REQUIRED.
C
         DO 560 K = 1, NEQ
            RWORK(NRT1+K) = Y(K)
            Y(K) = RWORK(NRT3+K)
            RWORK(NRT2+K) = YP(K)
            YP(K) = RWORK(NRT4+K)
            RWORK(NRT3+K) = YDP(K)
  560    CONTINUE
C
C        SINCE THE 12(10) PAIR DOES NOT USE FSAL (SEE ABOVE), A FUNCTION
C        EVALUATION IS DONE HERE WHICH WILL BE USED ON THE NEXT STEP.
C        FOR THE 6(4) PAIR SET YDP FROM THE LAST STAGE.
C
         IF (RHIGH.EQ.ONE) THEN
            CALL F(NEQ,T,Y,YDP)
         ELSE
            JSTAGE = NBASE + PTR(5)*NEQ
            DO 580 K = 1, NEQ
               YDP(K) = RWORK(JSTAGE+K)
  580       CONTINUE
         END IF
C
C        IF THE CODE HAS NOT REACHED THE END OF THE INTEGRATION RANGE
C        OR IS OPERATING IN INTERVAL MODE, ATTEMPT THE NEXT STEP.
C
         IF ( .NOT. (LAST .OR. (RWORK(SV1STP).EQ.ONE))) GO TO 40
C
C        END OF BLOCK G.
C
      END IF
C
C     BLOCK H.
C     SET DIAGNOSTIC INFORMATION AND EXIT.
C
  600 RWORK(SVHUSD) = HUSED
      RWORK(SVHNXT) = H
      RWORK(SVOKST) = DBLE(ISTP)
      RWORK(SVFLST) = DBLE(NFLSTP)
      RWORK(SVATST) = DBLE(NATT)
      RWORK(SVTOLD) = TOLD
  620 CONTINUE
C
C     TEST FOR SUCCESSIVE FAILURES AT THE SAME VALUE OF T.
C
      IF (IER.GT.0) THEN
         IF ( .NOT. LOOP) THEN
            LOOP = .TRUE.
            TCHECK = T
         ELSE
            IF (TCHECK.EQ.T) THEN
               NREC = NREC + 1
               WRITE (REC(NREC),FMT=99999) T
               IER = 4
            ELSE
               LOOP = .FALSE.
            END IF
         END IF
      ELSE
         LOOP = .FALSE.
      END IF
C
C     END OF BLOCK H.
C
      ERSTAT = IER + 1
      CALL D02LAX(ERSTAT,SET)
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
C
99999 FORMAT (' ** Two successive errors detected at the current value',
     *       ' of T,',1P,D13.5,'.')
99998 FORMAT (' ** Routine called with T = TEND. The value of T is',1P,
     *       D13.5,'.')
99997 FORMAT (' ** TEND has been reset such that the direction of inte',
     *       'gration is reversed.')
99996 FORMAT (' ** To satisfy the accuracy requirements the step size,',
     *       1P,D13.5,',')
99995 FORMAT ('    at T =',1P,D13.5,', is too small for the machine pr',
     *       'ecision.')
99994 FORMAT (' ** The maximum number of steps,',I16,', has been attem',
     *       'pted.')
99993 FORMAT (' ** The setup routine D02LXF has not been called.')
99992 FORMAT (' ** The previous call to the setup routine D02LXF resul',
     *       'ted in the')
99991 FORMAT ('    error exit IFAIL = ',I2,' .')
99990 FORMAT (' ** Inefficiency detected in integrating exactly to val',
     *       'ues of TEND.')
99989 FORMAT (' ** The value of NEQ supplied to D02LAF is not that sup',
     *       'plied to D02LXF:')
99988 FORMAT ('    NEQ in D02LXF was ',I16,', NEQ in D02LAF is ',I16,
     *       '.')
99987 FORMAT (' ** The dimension of RWORK supplied to D02LAF is not th',
     *       'at supplied to D02LXF:')
99986 FORMAT ('    LRWORK in D02LXF was ',I16,', LRWORK in D02LAF is ',
     *       I16,'.')
      END
