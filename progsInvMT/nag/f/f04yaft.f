      SUBROUTINE F04YAF(JOB,P,SIGMA,A,NRA,SVD,IRANK,SV,CJ,WORK,IFAIL)
C     MARK 11 RELEASE.  NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     F04YAF RETURNS ELEMENTS OF THE ESTIMATED VARIANCE-COVARIANCE
C     MATRIX OF THE SAMPLE REGRESSION COEFFICIENTS FOR THE SOLUTION
C     OF A LINEAR REGRESSION PROBLEM.
C
C     THE ROUTINE CAN BE USED TO FIND THE ESTIMATED VARIANCES OF THE
C     SAMPLE REGRESSION COEFFICIENTS.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 29-NOVEMBER-1982.  S.J. HAMMARLING.
C
C     NAG FORTRAN 66 ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04YAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGMA
      INTEGER           IFAIL, IRANK, JOB, NRA, P
      LOGICAL           SVD
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,P), CJ(P), SV(P), WORK(P)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, SCALE, SSQ, ZERO
      INTEGER           I, IDIAG, IERR, J, JJ, K, KK, KM1
      LOGICAL           OK
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, F06FCF, F06FDF, F06FJF, DCOPY, DAXPY,
     *                  F04YAY, X02ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ONE/1.0D+0/, ZERO/0.0D+0/
C     .. Executable Statements ..
C
C     CHECK THE INPUT PARAMETERS.
C
      OK = (P.GE.1) .AND. (SIGMA.GE.ZERO) .AND. (JOB.GE.(-1))
     *     .AND. (JOB.LE.P)
      IF (SVD) GO TO 20
      OK = (OK) .AND. (NRA.GE.P)
      GO TO 40
   20 CONTINUE
      OK = (OK) .AND. ((IRANK.GE.0) .AND. (IRANK.LE.P)
     *     .AND. (((JOB.GE.0) .AND. (NRA.GE.MAX(1,IRANK)))
     *     .OR. ((JOB.EQ.(-1)) .AND. (NRA.GE.P))))
   40 CONTINUE
      IF (OK) GO TO 60
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   60 CONTINUE
C
      IF ( .NOT. SVD) GO TO 100
      IF (IRANK.GT.0) GO TO 80
      IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      RETURN
   80 CONTINUE
  100 CONTINUE
C
      CALL X02ZAZ
C
C     DEAL WITH CASE WHERE SIGMA = ZERO. SET THE RELEVENT ELEMENTS
C     TO ZERO.
C
      IF (SIGMA.GT.ZERO) GO TO 180
      IF (JOB.EQ.(-1)) GO TO 120
      CALL F06FBF(P,ZERO,CJ,1)
      GO TO 160
  120 CONTINUE
      DO 140 J = 1, P
         CALL F06FBF(J,ZERO,A(1,J),1)
  140 CONTINUE
  160 CONTINUE
      IFAIL = 0
      RETURN
  180 CONTINUE
C
      IF (JOB.GT.(-1)) GO TO 420
C
C        FORM THE UPPER TRIANGULAR PART OF C. ( JOB = -1. )
C
      IF (SVD) GO TO 280
C
C           C = ( SIGMA**2 )*( U**( -1 ) )*( U**( -T ) )
C
C             = V*( V**T ),    V = SIGMA*( U**( -1 ) ).
C
C           FIRST FORM V.
C
      J = P
      DO 200 JJ = 1, P
         CALL F06FBF(J,ZERO,WORK,1)
         WORK(J) = SIGMA
C
         IDIAG = 1
         CALL F04YAY(1,J,A,NRA,WORK,IDIAG)
         IF (IDIAG.NE.0) GO TO 740
C
         CALL DCOPY(J,WORK,1,A(1,J),1)
         J = J - 1
  200 CONTINUE
C
C           NOW FORM C = V*( V**T ).
C
      DO 260 K = 1, P
         CALL DCOPY(K,A(1,K),1,WORK,1)
         IF (K.EQ.1) GO TO 240
         KM1 = K - 1
         DO 220 J = 1, KM1
            CALL DAXPY(J,WORK(J),WORK,1,A(1,J),1)
  220    CONTINUE
  240    CONTINUE
         CALL F06FDF(K,WORK(K),WORK,1,A(1,K),1)
  260 CONTINUE
      GO TO 400
  280 CONTINUE
C
C           C = ( SIGMA**2 )*P*( S**( -2 ) )*( P**T )
C
C             = ( V**T )*V,    V = SIGMA*( S**( -1 ) )*( P**T ).
C
C           FIRST FORM V.
C
      DO 300 I = 1, IRANK
         IF (SV(I).EQ.ZERO) GO TO 760
         WORK(I) = SIGMA/SV(I)
  300 CONTINUE
      DO 320 J = 1, P
         CALL F06FCF(IRANK,WORK,1,A(1,J),1)
  320 CONTINUE
C
C           NOW FORM C = ( V**T )*V.
C
      K = P
      DO 380 KK = 1, P
         CALL DCOPY(IRANK,A(1,K),1,WORK,1)
         IF (K.EQ.1) GO TO 360
         KM1 = K - 1
         DO 340 I = 1, KM1
            A(I,K) = DDOT(IRANK,A(1,I),1,WORK,1)
  340    CONTINUE
  360    CONTINUE
         A(K,K) = DDOT(IRANK,WORK,1,WORK,1)
         K = K - 1
  380 CONTINUE
  400 CONTINUE
      GO TO 720
  420 CONTINUE
C
      IF (JOB.GT.0) GO TO 620
C
C        FORM THE DIAGONAL ELEMENTS OF C. ( JOB = 0. )
C
      IF (SVD) GO TO 480
C
C           C = V*( V**T ),    V = SIGMA*( U**( -1 ) ).
C
C                         P
C           C( I, I ) =  SUM  ( V( I, J )**2 ).
C                       J = I
C
      CALL F06FBF(P,ZERO,CJ,1)
      DO 460 J = 1, P
         CALL F06FBF(J,ZERO,WORK,1)
         WORK(J) = SIGMA
C
         IDIAG = 1
         CALL F04YAY(1,J,A,NRA,WORK,IDIAG)
         IF (IDIAG.NE.0) GO TO 740
C
         DO 440 I = 1, J
            CJ(I) = CJ(I) + WORK(I)**2
  440    CONTINUE
  460 CONTINUE
      GO TO 600
  480 CONTINUE
C
C           C = ( V**T )*V,    V = SIGMA*( S**( -1 ) )*( P**T ).
C
C           C( J, J ) = ( V( J )**T )*V( J ).
C
      DO 580 J = 1, P
         IF (J.GT.1) GO TO 520
         DO 500 I = 1, IRANK
            IF (SV(I).EQ.ZERO) GO TO 760
            WORK(I) = A(I,1)/SV(I)
  500    CONTINUE
         GO TO 560
  520    CONTINUE
         DO 540 I = 1, IRANK
            WORK(I) = A(I,J)/SV(I)
  540    CONTINUE
  560    CONTINUE
         SCALE = ZERO
         SSQ = ONE
         CALL F06FJF(IRANK,WORK,1,SCALE,SSQ)
         CJ(J) = ((SIGMA*SCALE)**2)*SSQ
  580 CONTINUE
  600 CONTINUE
      GO TO 720
  620 CONTINUE
C
C        FORM THE JTH COLUMN OF C. ( JOB = J .GT. 0. )
C
      IF (SVD) GO TO 640
C
C           C = V*( V**T ),    V = SIGMA*( U**( -1 ) )
C
C           LET  V = ( W( 1 )  W( 2 ) ... W( P ) ).  THEN
C
C           ( U**T )*W( J ) = SIGMA*E( J ),   U*C( J ) = SIGMA*W( J ).
C
      CALL F06FBF(P,ZERO,CJ,1)
      CJ(JOB) = SIGMA**2
C
      IDIAG = 1
      CALL F04YAY(-1,P,A,NRA,CJ,IDIAG)
      IF (IDIAG.NE.0) GO TO 740
C
      IDIAG = 1
      CALL F04YAY(1,P,A,NRA,CJ,IDIAG)
      IF (IDIAG.NE.0) GO TO 740
C
      GO TO 700
  640 CONTINUE
C
C           C = ( V**T )*V,    V = SIGMA*( S**( -1 ) )*( P**T ).
C
C           C( J ) = ( V**T )*V( J ).
C
      DO 660 I = 1, IRANK
         IF (SV(I).EQ.ZERO) GO TO 760
         WORK(I) = ((SIGMA/SV(I))**2)*A(I,JOB)
  660 CONTINUE
      DO 680 J = 1, P
         CJ(J) = DDOT(IRANK,A(1,J),1,WORK,1)
  680 CONTINUE
  700 CONTINUE
  720 CONTINUE
      IFAIL = 0
      RETURN
C
  740 IERR = 3
      GO TO 780
  760 IERR = 4
  780 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
C
C     END OF F04YAF.
C
      END
