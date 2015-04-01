      DOUBLE PRECISION FUNCTION G01EBF(TAIL,T,DF,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G01EBF  --  RETURNS THE PROBABILITY ASSOCIATED WITH EITHER
C     (OR BOTH) TAIL(S) OF THE STUDENTS T-DISTRIBUTION WITH DF DEGREES
C     OF FREEDOM THROUGH THE FUNCTION NAME.
C     IF TAIL = 'U' or 'u' THE UPPER TAIL PROBABILITY IS RETURNED,
C     IF TAIL = 'S' or 's' THE TWO TAIL (SIGNIFICANCE LEVEL)
C                                          PROBABILITY IS RETURNED,
C     IF TAIL = 'C' or 'c' THE TWO TAIL (CONFIDENCE INTERVAL)
C                                          PROBABILITY IS RETURNED,
C     IF TAIL = 'L' or 'l' THE LOWER TAIL PROBABILITY IS RETURNED.
C
C     USES NAG LIBRARY ROUTINES P01AAF, S15ABF, G01EEF
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01EBF')
      DOUBLE PRECISION                 ZERO, HALF, ONE, TWO
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,
     *                                 TWO=2.0D0)
      DOUBLE PRECISION                 TWENTY
      PARAMETER                        (TWENTY=20.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF, T
      INTEGER                          IFAIL
      CHARACTER*1                      TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, B, BETA, BT, PDF, QBETA, TOL,
     *                                 TT, Y, Z
      INTEGER                          IERR, IFAIL2
      LOGICAL                          CCTAIL, LTAIL, STAIL, UTAIL
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S01BAF, S15ABF
      INTEGER                          P01ABF
      EXTERNAL                         S01BAF, S15ABF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G01EEF
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
      G01EBF = ZERO
      CCTAIL = .FALSE.
      LTAIL = .FALSE.
      STAIL = .FALSE.
      UTAIL = .FALSE.
      IF (TAIL.EQ.'C' .OR. TAIL.EQ.'c') THEN
         CCTAIL = .TRUE.
      ELSE IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
         LTAIL = .TRUE.
      ELSE IF (TAIL.EQ.'S' .OR. TAIL.EQ.'s') THEN
         STAIL = .TRUE.
      ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
         UTAIL = .TRUE.
      ELSE
         IERR = 1
         WRITE (REC,FMT=99998) TAIL
         GO TO 20
      END IF
      IF (DF.GE.ONE) THEN
         IERR = 0
         IFAIL2 = 1
         G01EBF = ZERO
         IF (DF.LT.TWENTY .AND. T.LT.ZERO) THEN
C
C           USE BETA FUNCTION FOR T.LT.0 AND DF .LT. 20
C
            TOL = 1.0D-5
            BT = DF/(DF+T**2)
            CALL G01EEF(BT,DF/TWO,HALF,TOL,BETA,QBETA,PDF,IFAIL2)
            IF (STAIL) THEN
               G01EBF = BETA
            ELSE IF (UTAIL) THEN
               G01EBF = HALF + HALF*QBETA
            ELSE IF (LTAIL) THEN
               G01EBF = HALF*BETA
            ELSE
               G01EBF = QBETA
            END IF
         ELSE IF (DF.LT.TWENTY .AND. T.GE.ZERO) THEN
C
C           USE BETA FUNCTION FOR T.GE.0 AND DF .LT. 20
C
            TOL = 1.0D-5
            BT = DF/(DF+T**2)
            CALL G01EEF(BT,DF/TWO,HALF,TOL,BETA,QBETA,PDF,IFAIL2)
            IF (STAIL) THEN
               G01EBF = BETA
            ELSE IF (UTAIL) THEN
               G01EBF = HALF*BETA
            ELSE IF (LTAIL) THEN
               G01EBF = HALF + HALF*QBETA
            ELSE
               G01EBF = QBETA
            END IF
C
C   USE HILL'S METHOD FOR DF .GE. 20
C
         ELSE
            Z = ONE
            TT = T*T
            Y = TT/DF
            Y = S01BAF(Y,IFAIL)
            A = DF - HALF
            B = 48D0*A*A
            Y = Y*A
            Y = (((((-0.4D0*Y-3.3D0)*Y-24.0D0)*Y-85.5D0)
     *          /(0.8D0*Y*Y+100.0D0+B)+Y+3.0D0)/B+ONE)*SQRT(Y)
            G01EBF = S15ABF(-Y,IFAIL2)
            IF (LTAIL .AND. T.GT.ZERO) THEN
               G01EBF = ONE - G01EBF
            ELSE IF (UTAIL .AND. T.LT.ZERO) THEN
               G01EBF = ONE - G01EBF
            ELSE IF (STAIL) THEN
               G01EBF = TWO*G01EBF
            ELSE IF (CCTAIL) THEN
               G01EBF = ONE - TWO*G01EBF
            END IF
         END IF
      ELSE
         WRITE (REC,FMT=99999) DF
         IERR = 2
      END IF
   20 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, DF.lt.1.0: DF = ',1P,D13.5)
99998 FORMAT (1X,'** On entry, TAIL is not valid. TAIL = ',A1)
      END
