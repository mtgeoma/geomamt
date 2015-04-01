      SUBROUTINE S17DEF(FNU,Z,N,SCALE,CY,NZ,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-760 (DEC 1989).
C
C     Original name: CBESJ
C
C     PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
C
C     DESCRIPTION
C     ===========
C
C         ON SCALE='U', S17DEF COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(I)=J(FNU+I-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
C         -PI.LT.ARG(Z).LE.PI. ON SCALE='S', S17DEF RETURNS THE SCALED
C         FUNCTIONS
C
C         CY(I)=EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C
C         WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
C         LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
C         ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS
C         (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),  -PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0E0
C           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
C                    SCALE= 'U' OR 'u' RETURNS
C                              CY(I)=J(FNU+I-1,Z), I=1,...,N
C                         = 'S' OR 's' RETURNS
C                              CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=J(FNU+I-1,Z)  OR
C                    CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
C                    DEPENDING ON SCALE, Y=AIMAG(Z).
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , LAST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I = N-NZ+1,...,N
C           IFAIL  - ERROR FLAG
C                   IFAIL=0, NORMAL RETURN - COMPUTATION COMPLETED
C                   IFAIL=1, INPUT ERROR   - NO COMPUTATION
C                   IFAIL=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
C                            TOO LARGE ON SCALE='U'
C                   IFAIL=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                   IFAIL=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                   IFAIL=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C     LONG DESCRIPTION
C     ================
C
C         THE COMPUTATION IS CARRIED OUT BY THE FORMULA
C
C         J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z).GE.0.0
C
C         J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z).LT.0.0
C
C         WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR NEGATIVE ORDERS,THE FORMULA
C
C              J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)
C
C         CAN BE USED. HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE
C         THE FUNCTION CHANGES RADICALLY. WHEN FNU IS A LARGE POSITIVE
C         INTEGER,THE MAGNITUDE OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A
C         LARGE NEGATIVE POWER OF TEN. BUT WHEN FNU IS NOT AN INTEGER,
C         Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
C         TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY
C         UNIT ROUNDOFF FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN
C         OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU. HERE,
C         LARGE MEANS FNU.GT.CABS(Z).
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IFAIL=3 IS TRIGGERED WHERE UR=X02AJF()=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IFAIL=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
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
      PARAMETER         (SRNAME='S17DEF')
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  FNU
      INTEGER           IFAIL, N, NZ
      CHARACTER*1       SCALE
C     .. Array Arguments ..
      COMPLEX*16        CY(N)
C     .. Local Scalars ..
      COMPLEX*16        CI, CSGN, ZN
      DOUBLE PRECISION  AA, ALIM, ARG, ASCLE, ATOL, AZ, BB, DIG, ELIM,
     *                  FN, FNUL, HPI, R1, R1M5, R2, RL, RTOL, TOL, YY
      INTEGER           I, IERR, INU, INUH, IR, K, K1, K2, KODE, NL,
     *                  NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AHF, X02AJF
      INTEGER           P01ABF, X02BBF, X02BHF, X02BJF, X02BKF, X02BLF
      EXTERNAL          X02AHF, X02AJF, P01ABF, X02BBF, X02BHF, X02BJF,
     *                  X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          S17DEZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, DCONJG, COS, INT, LOG10,
     *                  MAX, MIN, MOD, DBLE, SIN, SQRT, EXP
C     .. Data statements ..
      DATA              HPI/1.57079632679489662D0/
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      NZ = 0
      IF (SCALE.EQ.'U' .OR. SCALE.EQ.'u') THEN
         KODE = 1
      ELSE IF (SCALE.EQ.'S' .OR. SCALE.EQ.'s') THEN
         KODE = 2
      ELSE
         KODE = -1
      END IF
      IF (FNU.LT.0.0D0) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) FNU
      ELSE IF (KODE.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) SCALE
      ELSE IF (N.LT.1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99997) N
      END IF
      IF (IERR.EQ.0) THEN
C        ---------------------------------------------------------------
C        SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C        TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C        ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C        EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C        EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C        UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C        RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE
C        Z.
C        DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C        FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE
C        FNU.
C        ---------------------------------------------------------------
         TOL = MAX(X02AJF(),1.0D-18)
         K1 = X02BKF()
         K2 = X02BLF()
         R1M5 = LOG10(DBLE(X02BHF()))
         K = MIN(ABS(K1),ABS(K2))
         ELIM = 2.303D0*(K*R1M5-3.0D0)
         K1 = X02BJF() - 1
         AA = R1M5*K1
         DIG = MIN(AA,18.0D0)
         AA = AA*2.303D0
         ALIM = ELIM + MAX(-AA,-41.45D0)
         RL = 1.2D0*DIG + 3.0D0
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
         CI = DCMPLX(0.0D0,1.0D0)
         YY = DIMAG(Z)
         AZ = ABS(Z)
C        ---------------------------------------------------------------
C        TEST FOR RANGE
C        ---------------------------------------------------------------
         AA = 0.5D0/TOL
         BB = X02BBF(1.0D0)*0.5D0
         AA = MIN(AA,BB,X02AHF(1.0D0))
         FN = FNU + N - 1
         IF (AZ.LE.AA) THEN
            IF (FN.LE.AA) THEN
               AA = SQRT(AA)
               IF (AZ.GT.AA) THEN
                  IERR = 3
                  NREC = 1
                  WRITE (REC,FMT=99996) AZ, AA
               ELSE IF (FN.GT.AA) THEN
                  IERR = 3
                  NREC = 1
                  WRITE (REC,FMT=99995) FN, AA
               END IF
C              ---------------------------------------------------------
C              CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF
C              SIGNIFICANCE WHEN FNU IS LARGE
C              ---------------------------------------------------------
               INU = INT(FNU)
               INUH = INU/2
               IR = INU - 2*INUH
               ARG = (FNU-INU+IR)*HPI
               R1 = COS(ARG)
               R2 = SIN(ARG)
               CSGN = DCMPLX(R1,R2)
               IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
C              ---------------------------------------------------------
C              ZN IS IN THE RIGHT HALF PLANE
C              ---------------------------------------------------------
               ZN = -Z*CI
               IF (YY.LT.0.0D0) THEN
                  ZN = -ZN
                  CSGN = DCONJG(CSGN)
                  CI = DCONJG(CI)
               END IF
               CALL S17DEZ(ZN,FNU,KODE,N,CY,NZ,RL,FNUL,TOL,ELIM,ALIM)
               IF (NZ.GE.0) THEN
                  NL = N - NZ
                  IF (NL.NE.0) THEN
                     RTOL = 1.0D0/TOL
                     ASCLE = EXP(-ELIM)*RTOL
                     DO 20 I = 1, NL
C                       CY(I) = CY(I)*CSGN
                        ZN = CY(I)
                        AA = DBLE(ZN)
                        BB = DIMAG(ZN)
                        ATOL = 1.0D0
                        IF (MAX(ABS(AA),ABS(BB)).LE.ASCLE) THEN
                           ZN = ZN*RTOL
                           ATOL = TOL
                        END IF
                        ZN = ZN*CSGN
                        CY(I) = ZN*ATOL
                        CSGN = CSGN*CI
   20                CONTINUE
                  END IF
               ELSE IF (NZ.EQ.-3) THEN
                  NZ = 0
                  IERR = 4
                  NREC = 1
                  WRITE (REC,FMT=99991) AZ, AA
               ELSE IF (NZ.EQ.(-2)) THEN
                  NZ = 0
                  IERR = 5
                  NREC = 1
                  WRITE (REC,FMT=99994)
               ELSE
                  NZ = 0
                  IERR = 2
                  NREC = 1
                  WRITE (REC,FMT=99993) YY, ELIM
               END IF
               IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
               RETURN
            ELSE
               NZ = 0
               IERR = 4
               NREC = 1
               WRITE (REC,FMT=99992) FN, AA
            END IF
         ELSE
            NZ = 0
            IERR = 4
            NREC = 1
            WRITE (REC,FMT=99991) AZ, AA
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, FNU .LT. 0: FNU = ',D13.5)
99998 FORMAT (1X,'** On entry, SCALE has an illegal value: SCALE = ''',
     *       A,'''')
99997 FORMAT (1X,'** On entry, N .LE. 0: N = ',I16)
99996 FORMAT (1X,'** Results lack precision because abs(Z) =',1P,D13.5,
     *       ' .GT.',D13.5)
99995 FORMAT (1X,'** Results lack precision because FNU+N-1 =',1P,D13.5,
     *       ' .GT.',D13.5)
99994 FORMAT (1X,'** No computation - algorithm termination condition ',
     *       'not met.')
99993 FORMAT (1X,'** No computation because imag(Z) =',1P,D13.5,' .GT.',
     *       D12.5,', SCALE = ''U''.')
99992 FORMAT (1X,'** No computation because FNU+N-1 =',1P,D13.5,' .GT.',
     *       D13.5)
99991 FORMAT (1X,'** No computation because abs(Z) =',1P,D13.5,' .GT.',
     *       D13.5)
      END
