      SUBROUTINE S17DCF(FNU,Z,N,SCALE,CY,NZ,CWRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-759 (DEC 1989).
C
C     Original name: CBESY
C
C     PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C     DESCRIPTION
C     ===========
C
C         ON SCALE='U', S17DCF COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON SCALE='S', S17DCF RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0E0
C           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
C                   SCALE= 'U' OR 'u' RETURNS
C                             CY(I)=Y(FNU+I-1,Z), I=1,...,N
C                        = 'S' OR 's' RETURNS
C                             CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
C                             WHERE Y=AIMAG(Z)
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           CWRK   - A COMPLEX WORK VECTOR OF DIMENSION AT LEAST N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=Y(FNU+I-1,Z)  OR
C                    CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON PARAMETER SCALE.
C           NZ     - NZ=0 , A NORMAL RETURN
C                    NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
C                    UNDERFLOW (GENERALLY ON SCALE='S')
C           IFAIL  - ERROR FLAG
C                    IFAIL=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IFAIL=1, INPUT ERROR   - NO COMPUTATION
C                    IFAIL=2, OVERFLOW      - NO COMPUTATION, CABS(Z) IS
C                             TOO SMALL
C                    IFAIL=3, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                             TOO LARGE
C                    IFAIL=4, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DON
C                             BUT LOSSES OF SIGNIFICANCE BY ARGUMENT
C                             REDUCTION PRODUCE LESS THAN HALF OF MACHIN
C                             ACCURACY
C                    IFAIL=5, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                             TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI
C                             CANCE BY ARGUMENT REDUCTION
C                    IFAIL=6, ERROR              - NO COMPUTATION,
C                             ALGORITHM TERMINATION CONDITION NOT MET
C
C     LONG DESCRIPTION
C     ================
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C         Y(FNU,Z)=0.5*(H(1,FNU,Z)-H(2,FNU,Z))/I
C
C         WHERE I**2 = -1 AND THE HANKEL BESSEL FUNCTIONS H(1,FNU,Z)
C         AND H(2,FNU,Z) ARE CALCULATED IN S17DLF.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD
C         INTEGERS THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE
C         POSITIVE HALF ODD INTEGER,THE MAGNITUDE OF Y(-FNU,Z)=J(FNU,Z)*
C         SIN(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS
C         NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
C         LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM
C         CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT. THUS,
C         WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
C         ODD INTEGER. HERE, LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IFAIL=4 IS TRIGGERED WHERE UR=X02AJF()=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IFAIL=5. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
C         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
C         INTEGER, U3=X02BBF(). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
C         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
C         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
C         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
C         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
C         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
C         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
C         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
C         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
C
C         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
C         BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
C         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
C         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
C         ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
C         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
C         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
C         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
C         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
C         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
C         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
C         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
C         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
C         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
C         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
C         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
C         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
C         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
C         OR -PI/2+P.
C
C     REFERENCES
C     ==========
C               HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ
C                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
C                 COMMERCE, 1955.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
C
C               COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
C
C               A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
C                 1018, MAY, 1985
C
C               A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX
C                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
C                 MATH. SOFTWARE, 1986
C
C     DATE WRITTEN   830501   (YYMMDD)
C     REVISION DATE  830501   (YYMMDD)
C     AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='S17DCF')
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  FNU
      INTEGER           IFAIL, N, NZ
      CHARACTER*1       SCALE
C     .. Array Arguments ..
      COMPLEX*16        CWRK(N), CY(N)
C     .. Local Scalars ..
      COMPLEX*16        C1, C2, EX, HCI, ZU, ZV
      DOUBLE PRECISION  AA, ASCLE, ATOL, BB, ELIM, EY, R1, R2, RTOL,
     *                  TAY, TOL, XX, YY
      INTEGER           I, IERR, IFAIL1, IFAIL2, K, K1, K2, KODE, NZ1,
     *                  NZ2
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF, X02BHF, X02BKF, X02BLF
      EXTERNAL          P01ABF, X02AJF, X02BHF, X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          S17DLF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, DCONJG, COS, EXP, LOG10,
     *                  MIN, DBLE, SIN, MAX
C     .. Executable Statements ..
      XX = DBLE(Z)
      YY = DIMAG(Z)
      IERR = 0
      NZ = 0
      IF (SCALE.EQ.'U' .OR. SCALE.EQ.'u') THEN
         KODE = 1
      ELSE IF (SCALE.EQ.'S' .OR. SCALE.EQ.'s') THEN
         KODE = 2
      ELSE
         KODE = -1
      END IF
      IF (XX.EQ.0.0D0 .AND. YY.EQ.0.0D0) THEN
         IERR = 1
         WRITE (REC,FMT=99999)
      ELSE IF (FNU.LT.0.0D0) THEN
         IERR = 1
         WRITE (REC,FMT=99998) FNU
      ELSE IF (KODE.EQ.-1) THEN
         IERR = 1
         WRITE (REC,FMT=99997) SCALE
      ELSE IF (N.LT.1) THEN
         IERR = 1
         WRITE (REC,FMT=99996) N
      END IF
      IF (IERR.EQ.0) THEN
         IF (IFAIL.EQ.1) THEN
            IFAIL1 = 1
            IFAIL2 = 1
         ELSE
            IFAIL1 = -13
            IFAIL2 = -13
         END IF
         HCI = DCMPLX(0.0D0,0.5D0)
         CALL S17DLF(1,FNU,Z,N,SCALE,CY,NZ1,IFAIL1)
         IF (IFAIL1.EQ.0 .OR. IFAIL1.EQ.4) THEN
            CALL S17DLF(2,FNU,Z,N,SCALE,CWRK,NZ2,IFAIL2)
            IF (IFAIL2.EQ.0 .OR. IFAIL2.EQ.4) THEN
               IF (IFAIL1.EQ.4) IFAIL2 = 4
               NZ = MIN(NZ1,NZ2)
               IF (KODE.EQ.2) THEN
                  TOL = MAX(X02AJF(),1.0D-18)
                  K1 = X02BKF()
                  K2 = X02BLF()
                  K = MIN(ABS(K1),ABS(K2))
C                 ------------------------------------------------------
C                 ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND
C                 OVERFLOW LIMIT
C                 ------------------------------------------------------
                  ELIM = 2.303D0*(K*LOG10(DBLE(X02BHF()))-3.0D0)
                  R1 = COS(XX)
                  R2 = SIN(XX)
                  EX = DCMPLX(R1,R2)
                  EY = 0.0D0
                  TAY = ABS(YY+YY)
                  IF (TAY.LT.ELIM) EY = EXP(-TAY)
                  IF (YY.LT.0.0D0) THEN
                     C1 = EX
                     C2 = DCONJG(EX)*DCMPLX(EY,0.0D0)
                  ELSE
                     C1 = EX*DCMPLX(EY,0.0D0)
                     C2 = DCONJG(EX)
                  END IF
                  NZ = 0
                  RTOL = 1.0D0/TOL
                  ASCLE = EXP(-ELIM)*RTOL
                  DO 20 I = 1, N
C                    CY(I) = HCI*(C2*CWRK(I)-C1*CY(I))
                     ZV = CWRK(I)
                     AA = DBLE(ZV)
                     BB = DIMAG(ZV)
                     ATOL = 1.0D0
                     IF (MAX(ABS(AA),ABS(BB)).LE.ASCLE) THEN
                        ZV = ZV*RTOL
                        ATOL = TOL
                     END IF
                     ZV = ZV*C2*HCI
                     ZV = ZV*ATOL
                     ZU = CY(I)
                     AA = DBLE(ZU)
                     BB = DIMAG(ZU)
                     ATOL = 1.0D0
                     IF (MAX(ABS(AA),ABS(BB)).LE.ASCLE) THEN
                        ZU = ZU*RTOL
                        ATOL = TOL
                     END IF
                     ZU = ZU*C1*HCI
                     ZU = ZU*ATOL
                     CY(I) = ZV - ZU
                     IF (CY(I).EQ.DCMPLX(0.0D0,0.0D0) .AND. EY.EQ.0.0D0)
     *                   NZ = NZ + 1
   20             CONTINUE
               ELSE
                  DO 40 I = 1, N
                     CY(I) = HCI*(CWRK(I)-CY(I))
   40             CONTINUE
               END IF
            END IF
            IFAIL = P01ABF(IFAIL,IFAIL2,SRNAME,0,REC)
         ELSE
            IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
         END IF
      ELSE
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      END IF
      RETURN
C
99999 FORMAT (1X,'** On entry, Z = (0.0,0.0)')
99998 FORMAT (1X,'** On entry, FNU .LT. 0: FNU = ',D13.5)
99997 FORMAT (1X,'** On entry, SCALE has an illegal value: SCALE = ''',
     *       A,'''')
99996 FORMAT (1X,'** On entry, N .LE. 0: N = ',I16)
      END
