      SUBROUTINE G08AHF(N1,X,N2,Y,TAIL,U,UNOR,P,TIES,RANKS,WRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G08AHF performs the Mann-Whitney U test on two independent
C     samples of possibly unequal size. The tail probaility is
C     returned via P and corresponds to the TAIL option chosen.
C     P is based on the normal approximation. To calculate the
C     exact probability routine G08AHY or G08AHZ must be used.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, U, UNOR
      INTEGER           IFAIL, N1, N2
      LOGICAL           TIES
      CHARACTER*1       TAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  RANKS(N1+N2), WRK(N1+N2), X(N1), Y(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  FN1, FN2, FNT, FNX, R1, S, TS, WN, Z
      INTEGER           I, IERROR, IFA, NT
      LOGICAL           LOWER, TOTAIL, UPPER
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          S15ABF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AGY
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, DBLE, SQRT
C     .. Executable Statements ..
      TOTAIL = TAIL .EQ. 'T' .OR. TAIL .EQ. 't'
      UPPER = TAIL .EQ. 'U' .OR. TAIL .EQ. 'u'
      LOWER = TAIL .EQ. 'L' .OR. TAIL .EQ. 'l'
      IF (N1.LT.1 .OR. N2.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N1, N2
      ELSE IF ( .NOT. TOTAIL .AND. .NOT. UPPER .AND. .NOT. LOWER) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) TAIL
      ELSE
         IERROR = 0
         NT = N1 + N2
         DO 20 I = 1, N1
            WRK(I) = X(I)
   20    CONTINUE
         DO 40 I = 1, N2
            WRK(N1+I) = Y(I)
   40    CONTINUE
         TIES = .FALSE.
         CALL G08AGY(WRK,RANKS,NT,1,TS)
         TIES = TS .GT. 0.0D0
         R1 = 0.0D0
         DO 60 I = 1, N1
            R1 = R1 + RANKS(I)
   60    CONTINUE
         FNT = DBLE(NT)
         FN1 = DBLE(N1)
         FN2 = DBLE(N2)
         FNX = FN1*FN2
         U = R1 - DBLE(N1*(N1+1)/2)
         IF (TIES) THEN
            S = (FNX*(FNT+1.0D0)/3.0D0) - 4.0D0*FNX*TS/(FNT*(FNT-1.0D0))
            IF (S.LT.X02AMF()) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99997) S
               GO TO 80
            ELSE
               S = SQRT(S)
            END IF
         ELSE
            S = SQRT(FNX*(FNT+1.0D0)/3.0D0)
         END IF
         WN = 2.0D0*U - FNX
C
C        Normalised statistic with continuity correction
C
         IF (WN.GT.0.0D0) THEN
            UNOR = (WN-1.0D0)/S
         ELSE IF (WN.LT.0.0D0) THEN
            UNOR = (WN+1.0D0)/S
         ELSE
            UNOR = 0.0D0
         END IF
C
C           Normal approximation to tail probability
C
         Z = UNOR
         IF (LOWER .AND. WN.LT.0.0D0) THEN
            P = S15ABF(Z,IFA)
         ELSE IF (UPPER .AND. WN.GT.0.0D0) THEN
            P = 1.0D0 - S15ABF(Z,IFA)
         ELSE IF (LOWER .AND. WN.GE.0.0D0) THEN
            IFA = 1
            Z = (WN+1.0D0)/S
            P = S15ABF(Z,IFA)
         ELSE IF (UPPER .AND. WN.LE.0.0D0) THEN
            IFA = 1
            Z = (WN-1.0D0)/S
            P = 1.0D0 - S15ABF(Z,IFA)
         ELSE IF (TOTAIL) THEN
            P = S15ABF(Z,IFA)
            P = 2.0D0*(MIN(P,1.0D0-P))
         END IF
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
99999 FORMAT (1X,'** On entry, N1.lt.1 or N2.lt.1, N1 = ',I16,
     *       ' , N2 = ',I16)
99998 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99997 FORMAT (1X,'** Variance of the sample = ',D13.5)
      END
