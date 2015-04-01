      SUBROUTINE G02HMF(UCV,USERP,INDM,N,M,X,LDX,COV,A,WT,THETA,BL,BD,
     *                  MAXIT,NITMON,TOL,NIT,WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     ITERATIVE ALGORITHM FOR THE COMPUTATION OF THE MATRIX A
C     AND ROBUST ESIMATE OF COVARIANCE
C     STAHEL'S  METHOD
C     BASED ON ROUTINES IN ROBETH BY A. MARAZZI
C
C
C     PARAMETER CHECK AND INITIALIZATION
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02HMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BD, BL, TOL
      INTEGER           IFAIL, INDM, LDX, M, MAXIT, N, NIT, NITMON
C     .. Array Arguments ..
      DOUBLE PRECISION  A(M*(M+1)/2), COV(M*(M+1)/2), THETA(M),
     *                  USERP(*), WK(2*M), WT(N), X(LDX,M)
C     .. Subroutine Arguments ..
      EXTERNAL          UCV
C     .. Local Scalars ..
      DOUBLE PRECISION  BMAX, CX, SMAX, SV, SW, TEMP, U, UMAX, W, XN,
     *                  ZNR
      INTEGER           I, IERROR, IFLAG, IJ, IL, J, J1, J2, JJ, JJ1,
     *                  JJ2, JL, KK, L, MM, MM4, NOUT, NREC
      LOGICAL           AZERO
      CHARACTER*80      REC
C     .. Local Arrays ..
      CHARACTER*80      EREC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           P01ABF
      EXTERNAL          DNRM2, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02AAY, G02AAZ, DCOPY, DTPMV, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, MOD, DBLE
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (EREC(1),FMT=99999) M
      ELSE IF (N.LT.2) THEN
         WRITE (EREC(1),FMT=99998) N
      ELSE IF (N.LT.M) THEN
         WRITE (EREC(1),FMT=99997) N, M
      ELSE IF (LDX.LT.N) THEN
         WRITE (EREC(1),FMT=99996) LDX, N
      ELSE
         IERROR = 2
      END IF
      IF (IERROR.NE.1) THEN
         AZERO = .FALSE.
         IJ = 0
         DO 20 I = 1, M
            IJ = IJ + I
            IF (A(IJ).EQ.0.0D0) GO TO 40
   20    CONTINUE
         GO TO 60
C
   40    AZERO = .TRUE.
   60    IF (AZERO) THEN
            WRITE (EREC(1),FMT=99995) I
         ELSE IF (TOL.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99994) TOL
         ELSE IF (MAXIT.LE.0) THEN
            WRITE (EREC(1),FMT=99993) MAXIT
         ELSE IF (BL.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99992) BL
         ELSE IF (BD.LE.0.0D0) THEN
            WRITE (EREC(1),FMT=99991) BD
         ELSE
            IERROR = 0
         END IF
         IF (IERROR.NE.2) THEN
            DO 100 J = 1, M
               CX = X(1,J)
               DO 80 I = 2, N
                  IF (CX.NE.X(I,J)) GO TO 100
   80          CONTINUE
               GO TO 560
C
  100       CONTINUE
C
C           SET UP CONSTANTS ETC
C
            NIT = 0
            XN = DBLE(N)
            MM = (M+1)*M/2
C
C           SET UP FOR ITERATION MONITORING
C
            IF (NITMON.GT.0) THEN
               IFLAG = 0
               CALL X04ABF(IFLAG,NOUT)
               WRITE (REC,FMT=99988)
               CALL X04BAF(NOUT,REC)
            END IF
            DO 120 I = 1, N
               WT(I) = 0.0D0
  120       CONTINUE
  140       CONTINUE
            DO 160 I = 1, M
               WK(I) = 0.0D0
  160       CONTINUE
            DO 180 I = 1, MM
               COV(I) = 0.0D0
  180       CONTINUE
C
C           ITERATIONS
C
            NIT = NIT + 1
            UMAX = 0.0D0
C
C           COMPUTE AVERAGES
C
            SV = 0.0D0
            SW = 0.0D0
            DO 280 L = 1, N
               DO 200 J = 1, M
                  WK(M+J) = X(L,J) - THETA(J)
  200          CONTINUE
               CALL DTPMV('U','T','N',M,A,WK(M+1),1)
               ZNR = DNRM2(M,WK(M+1),1)
               CALL UCV(ZNR,USERP,U,W)
               IF (U.LT.0.0D0) THEN
                  GO TO 540
C
               ELSE IF (W.LT.0.0D0) THEN
                  GO TO 520
C
               ELSE
                  TEMP = ABS(U-WT(L))
                  IF (TEMP.GT.UMAX) UMAX = TEMP
                  WT(L) = U
                  SW = SW + W
                  SV = SV + U
                  IJ = 0
                  DO 240 I = 1, M
                     DO 220 J = 1, I
                        IJ = IJ + 1
                        COV(IJ) = WK(M+I)*U*WK(M+J) + COV(IJ)
  220                CONTINUE
  240             CONTINUE
                  DO 260 J = 1, M
                     WK(J) = WK(J) + W*(X(L,J)-THETA(J))
  260             CONTINUE
               END IF
  280       CONTINUE
C
C           FIND IMPROVEMENT MATRIX (COV) FOR A;
C           TRUNCATE IF NECESSARY; FIND MAXIMUM IMPROVEMENT
C
            IF (INDM.NE.1) SV = XN
            SMAX = 0.0D0
            IF (SV.EQ.0.0D0) GO TO 480
            IJ = 0
            DO 320 I = 1, M
               DO 300 J = 1, I - 1
                  IJ = IJ + 1
                  COV(IJ) = -MIN(MAX(COV(IJ)/SV,-BL),BL)
                  SMAX = MAX(SMAX,ABS(COV(IJ)))
  300          CONTINUE
               IJ = IJ + 1
               CX = -MIN(MAX((COV(IJ)/SV-1.0D0)/2.0D0,-BD),BD)
               SMAX = MAX(SMAX,ABS(CX))
               COV(IJ) = CX + 1.0D0
  320       CONTINUE
            IF (SW.GT.0.0D0) THEN
               BMAX = 0.0D0
               DO 340 J = 1, M
                  WK(J) = WK(J)/SW
                  TEMP = ABS(WK(J))/(1.0D0+ABS(THETA(J)))
                  IF (TEMP.GT.BMAX) BMAX = TEMP
  340          CONTINUE
            ELSE
               IERROR = 6
               WRITE (EREC(1),FMT=99980)
               GO TO 580
            END IF
            IF (UMAX.LT.TOL .AND. SMAX.LT.TOL .AND. BMAX.LT.TOL)
     *          GO TO 500
C
C           FIND NEW TRANSFORMATION MATRIX A=COV*A
C
            IL = 1
            DO 360 I = 1, M
               CALL DTPMV('U','N','N',I,A,COV(IL),1)
               IL = IL + I
  360       CONTINUE
            CALL DCOPY(MM,COV,1,A,1)
            DO 380 J = 1, M
               THETA(J) = THETA(J) + WK(J)
  380       CONTINUE
C
C           ITERATION MONITORING
C
            IF (NITMON.GT.0) THEN
               IF ((MOD(NIT,NITMON).EQ.0) .OR. (NIT.EQ.1)) THEN
                  WRITE (REC,FMT=99982)
                  CALL X04BAF(NOUT,REC)
                  WRITE (REC,FMT=99987) NIT, UMAX
                  CALL X04BAF(NOUT,REC)
                  WRITE (REC,FMT=99983)
                  CALL X04BAF(NOUT,REC)
                  DO 400 J = 1, M
                     WRITE (REC,FMT=99986) J, THETA(J)
                     CALL X04BAF(NOUT,REC)
  400             CONTINUE
                  WRITE (REC,FMT=99984)
                  CALL X04BAF(NOUT,REC)
                  MM4 = MIN(5,M)
                  J2 = 0
                  DO 420 J = 1, MM4
                     J1 = J2 + 1
                     J2 = J2 + J
                     WRITE (REC,FMT=99985) (A(KK),KK=J1,J2)
                     CALL X04BAF(NOUT,REC)
  420             CONTINUE
                  DO 460 J = 6, M
                     J1 = J2 + 1
                     J2 = J2 + J
                     JL = (J2-J1)/5 + 1
                     DO 440 JJ = 1, JL
                        JJ1 = (JJ-1)*5 + J1
                        JJ2 = MIN(J2,JJ1+4)
                        WRITE (REC,FMT=99985) (A(KK),KK=JJ1,JJ2)
                        CALL X04BAF(NOUT,REC)
  440                CONTINUE
  460             CONTINUE
               END IF
            END IF
C
C           STOP ITERATIONS IF GIVEN TOLERANCE OR MAXIMUM
C           NUMBER OF ITERATIONS IS REACHED
C
            IF (NIT.LT.MAXIT) GO TO 140
            IERROR = 5
            WRITE (EREC(1),FMT=99989)
            GO TO 500
  480       WRITE (EREC(1),FMT=99978)
            IERROR = 6
            GO TO 580
C
C           CALCULATE INV(A)
C
  500       CONTINUE
            CALL DCOPY(MM,A,1,COV,1)
            CALL G02AAZ('U','N',M,COV)
C
C           CALCULATE INV(A)INV(A)' AND STORE IN COV
C
            CALL G02AAY('U','N',M,COV)
            GO TO 580
C
  520       IERROR = 4
            WRITE (EREC(1),FMT=99981) ZNR, W
            GO TO 580
C
  540       IERROR = 4
            WRITE (EREC(1),FMT=99990) ZNR, U
            GO TO 580
C
  560       IERROR = 3
            WRITE (EREC(1),FMT=99979) J
         END IF
      END IF
  580 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,EREC)
C
99999 FORMAT (' ** On entry, M.lt.1 : M =',I16)
99998 FORMAT (' ** On entry, N.lt.2 : N =',I16)
99997 FORMAT (' ** On entry, N.lt.M : N =',I16,' M=',I16)
99996 FORMAT (' ** On entry, LDX.lt.N : LDX =',I16,' N =',I16)
99995 FORMAT (' ** On entry, diagonal element ',I4,' of A is 0.0')
99994 FORMAT (' ** On entry, TOL.le.0.0 : TOL =',D13.5)
99993 FORMAT (' ** On entry, MAXIT.le.0 : MAXIT =',I16)
99992 FORMAT (' ** On entry, BL.le.0 : BL =',D13.5)
99991 FORMAT (' ** On entry, BD.le.0 : BD =',D13.5)
99990 FORMAT (' ** U value returned by UCV.lt.0.0 : U(',D13.5,') =',
     *       D13.5)
99989 FORMAT (' ** Iterations to calculate weights failed to converge')
99988 FORMAT ('                    ** ITERATION MONITORING **')
99987 FORMAT (' ITERATION ',I16,'  MAX DELTA = ',D13.5)
99986 FORMAT (' ',I16,'  ',D13.5)
99985 FORMAT (5('  ',D13.5))
99984 FORMAT (' MATRIX A')
99983 FORMAT ('                I       THETA(I)')
99982 FORMAT ('  ')
99981 FORMAT (' ** W value returned by UCV.lt.0.0 : W(',D13.5,') =',
     *       D13.5)
99980 FORMAT (' ** Sum of w''s (D1) is zero')
99979 FORMAT (' ** ',I16,' TH column of X has constant value')
99978 FORMAT (' ** Sum of u''s (D2) is zero')
      END
