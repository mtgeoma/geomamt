      SUBROUTINE S17DHF(DERIV,Z,SCALE,BI,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-780 (DEC 1989).
C
C     Original name: CBIRY
C
C     PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
C     DESCRIPTION
C     ===========
C
C         ON SCALE='U', S17DHF COMPUTES THE COMPLEX AIRY FUNCTION BI(Z)
C         OR ITS DERIVATIVE DBI(Z)/DZ ON DERIV='F' OR DERIV='D'
C         RESPECTIVELY. ON
C         SCALE='S', A SCALING OPTION CEXP(-AXZTA)*BI(Z) OR CEXP(-AXZTA)
C         *DBI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN
C         BOTH THE LEFT AND RIGHT HALF PLANES WHERE
C         ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA) AND AXZTA=ABS(XZTA).
C         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
C         MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y)
C           DERIV  - RETURN FUNCTION (DERIV='F') OR DERIVATIVE
C                    (DERIV='D')
C           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
C                    SCALE = 'U' OR SCALE = 'u' RETURNS
C                           BI=BI(Z)                  ON DERIV='F' OR
C                           BI=DBI(Z)/DZ              ON DERIV='D'
C                    SCALE = 'S' OR SCALE = 's' RETURNS
C                           BI=CEXP(-AXZTA)*BI(Z)     ON DERIV='F' OR
C                           BI=CEXP(-AXZTA)*DBI(Z)/DZ ON DERIV='D' WHERE
C                           ZTA=(2/3)*Z*CSQRT(Z)=CMPLX(XZTA,YZTA)
C                           AND AXZTA=ABS(XZTA)
C         OUTPUT
C           BI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR DERIV
C                    AND SCALE
C           IFAIL  - ERROR FLAG
C                   IFAIL=0, NORMAL RETURN - COMPUTATION COMPLETED
C                   IFAIL=1, INPUT ERROR   - NO COMPUTATION
C                   IFAIL=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
C                            TOO LARGE WITH SCALE = 'U'
C                   IFAIL=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
C                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
C                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
C                   IFAIL=4, CABS(Z) TOO LARGE  - NO COMPUTATION
C                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
C                            REDUCTION
C                   IFAIL=5, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C     LONG DESCRIPTION
C     ================
C
C         BI AND DBI ARE COMPUTED FOR CABS(Z).GT.1.0 FROM THE I BESSEL
C         FUNCTIONS BY
C
C                BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
C               DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
C                               C=1.0/SQRT(3.0)
C                               ZTA=(2/3)*Z**(3/2)
C
C         WITH THE POWER SERIES FOR CABS(Z).LE.1.0.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
C         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
C         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
C         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
C         FLAG IFAIL=3 IS TRIGGERED WHERE UR=X02AJF()=UNIT ROUNDOFF.
C         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
C         ALL SIGNIFICANCE IS LOST AND IFAIL=4. IN ORDER TO USE THE INT
C         FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
C         LARGEST INTEGER, U3=X02BBF(). THUS, THE MAGNITUDE OF ZETA
C         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
C         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
C         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
C         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
C         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
C         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
C         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
C         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
C         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
C         PRECISION ARITHMETIC.
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
      PARAMETER         (SRNAME='S17DHF')
C     .. Scalar Arguments ..
      COMPLEX*16        BI, Z
      INTEGER           IFAIL
      CHARACTER         DERIV, SCALE
C     .. Local Scalars ..
      COMPLEX*16        CONE, CSQ, S1, S2, TRM1, TRM2, Z3, ZTA
      DOUBLE PRECISION  AA, AD, AK, ALAZ, ALIM, ATRM, AZ, AZ3, BB, BK,
     *                  C1, C2, CK, COEF, D1, D2, DIG, DK, ELIM, FID,
     *                  FMR, FNU, FNUL, PI, R1M5, RL, SAVAA, SFAC, TOL,
     *                  TTH, Z3I, Z3R, ZI, ZR
      INTEGER           ID, IERR, K, K1, K2, KODE, NREC, NZ
C     .. Local Arrays ..
      COMPLEX*16        CY(2)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AHF, X02AJF
      INTEGER           P01ABF, X02BBF, X02BHF, X02BJF, X02BKF, X02BLF
      EXTERNAL          X02AHF, X02AJF, P01ABF, X02BBF, X02BHF, X02BJF,
     *                  X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          S17DEZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, LOG, LOG10, MAX,
     *                  MIN, DBLE, SIN, SQRT
C     .. Data statements ..
      DATA              TTH, C1, C2, COEF, PI/6.66666666666666667D-01,
     *                  6.14926627446000736D-01,
     *                  4.48288357353826359D-01,
     *                  5.77350269189625765D-01,
     *                  3.14159265358979324D+00/
      DATA              CONE/(1.0D0,0.0D0)/
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      NZ = 0
      IF (DERIV.EQ.'F' .OR. DERIV.EQ.'f') THEN
         ID = 0
      ELSE IF (DERIV.EQ.'D' .OR. DERIV.EQ.'d') THEN
         ID = 1
      ELSE
         ID = -1
      END IF
      IF (SCALE.EQ.'U' .OR. SCALE.EQ.'u') THEN
         KODE = 1
      ELSE IF (SCALE.EQ.'S' .OR. SCALE.EQ.'s') THEN
         KODE = 2
      ELSE
         KODE = -1
      END IF
      IF (ID.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999) DERIV
      ELSE IF (KODE.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) SCALE
      END IF
      IF (IERR.EQ.0) THEN
         AZ = ABS(Z)
         TOL = MAX(X02AJF(),1.0D-18)
         FID = ID
         IF (AZ.GT.1.0D0) THEN
C           ------------------------------------------------------------
C           CASE FOR CABS(Z).GT.1.0
C           ------------------------------------------------------------
            FNU = (1.0D0+FID)/3.0D0
C           ------------------------------------------------------------
C           SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C           TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C           ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW
C           LIMIT.
C           EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C           EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS
C           NEAR UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC
C           IS DONE.
C           RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR
C           LARGE Z.
C           DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C           FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR
C           LARGE FNU.
C           ------------------------------------------------------------
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
            ALAZ = LOG(AZ)
C           ------------------------------------------------------------
C           TEST FOR RANGE
C           ------------------------------------------------------------
            AA = 0.5D0/TOL
            BB = X02BBF(1.0D0)*0.5D0
            AA = MIN(AA,BB,X02AHF(1.0D0))
            AA = AA**TTH
            IF (AZ.GT.AA) THEN
               IERR = 4
               NREC = 1
               NZ = 0
               WRITE (REC,FMT=99997) AZ, AA
            ELSE
               AA = SQRT(AA)
               SAVAA = AA
               IF (AZ.GT.AA) THEN
                  IERR = 3
                  NREC = 1
                  WRITE (REC,FMT=99996) AZ, AA
               END IF
               CSQ = SQRT(Z)
               ZTA = Z*CSQ*DCMPLX(TTH,0.0D0)
C              ---------------------------------------------------------
C              RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS
C              SMALL
C              ---------------------------------------------------------
               SFAC = 1.0D0
               ZI = DIMAG(Z)
               ZR = DBLE(Z)
               AK = DIMAG(ZTA)
               IF (ZR.LT.0.0D0) THEN
                  BK = DBLE(ZTA)
                  CK = -ABS(BK)
                  ZTA = DCMPLX(CK,AK)
               END IF
               IF (ZI.EQ.0.0D0 .AND. ZR.LE.0.0D0) ZTA = DCMPLX(0.0D0,AK)
               AA = DBLE(ZTA)
               IF (KODE.NE.2) THEN
C                 ------------------------------------------------------
C                 OVERFLOW TEST
C                 ------------------------------------------------------
                  BB = ABS(AA)
                  IF (BB.GE.ALIM) THEN
                     BB = BB + 0.25D0*LOG(AZ)
                     SFAC = TOL
                     IF (BB.GT.ELIM) GO TO 20
                  END IF
               END IF
               FMR = 0.0D0
               IF (AA.LT.0.0D0 .OR. ZR.LE.0.0D0) THEN
                  FMR = PI
                  IF (ZI.LT.0.0D0) FMR = -PI
                  ZTA = -ZTA
               END IF
C              ---------------------------------------------------------
C              AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
C              KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM S17DEZ
C              ---------------------------------------------------------
               CALL S17DEZ(ZTA,FNU,KODE,1,CY,NZ,RL,FNUL,TOL,ELIM,ALIM)
               IF (NZ.GE.0) THEN
                  AA = FMR*FNU
                  Z3 = DCMPLX(SFAC,0.0D0)
                  S1 = CY(1)*DCMPLX(COS(AA),SIN(AA))*Z3
                  FNU = (2.0D0-FID)/3.0D0
                  CALL S17DEZ(ZTA,FNU,KODE,2,CY,NZ,RL,FNUL,TOL,ELIM,
     *                        ALIM)
                  CY(1) = CY(1)*Z3
                  CY(2) = CY(2)*Z3
C                 ------------------------------------------------------
C                 BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
C                 ------------------------------------------------------
                  S2 = CY(1)*DCMPLX(FNU+FNU,0.0D0)/ZTA + CY(2)
                  AA = FMR*(FNU-1.0D0)
                  S1 = (S1+S2*DCMPLX(COS(AA),SIN(AA)))*DCMPLX(COEF,
     *                 0.0D0)
                  IF (ID.EQ.1) THEN
                     S1 = Z*S1
                     BI = S1*DCMPLX(1.0D0/SFAC,0.0D0)
                  ELSE
                     S1 = CSQ*S1
                     BI = S1*DCMPLX(1.0D0/SFAC,0.0D0)
                  END IF
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               ELSE IF (NZ.EQ.(-3)) THEN
                  NZ = 0
                  IERR = 4
                  NREC = 1
                  WRITE (REC,FMT=99997) AZ, SAVAA
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               ELSE IF (NZ.NE.(-1)) THEN
                  NZ = 0
                  IERR = 5
                  NREC = 1
                  WRITE (REC,FMT=99995)
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               END IF
   20          NZ = 0
               IERR = 2
               NREC = 1
               WRITE (REC,FMT=99994) ZR
            END IF
         ELSE
C           ------------------------------------------------------------
C           POWER SERIES FOR CABS(Z).LE.1.
C           ------------------------------------------------------------
            S1 = CONE
            S2 = CONE
            IF (AZ.LT.TOL) THEN
               AA = C1*(1.0D0-FID) + FID*C2
               BI = DCMPLX(AA,0.0D0)
            ELSE
               AA = AZ*AZ
               IF (AA.GE.TOL/AZ) THEN
                  TRM1 = CONE
                  TRM2 = CONE
                  ATRM = 1.0D0
                  Z3 = Z*Z*Z
                  AZ3 = AZ*AA
                  AK = 2.0D0 + FID
                  BK = 3.0D0 - FID - FID
                  CK = 4.0D0 - FID
                  DK = 3.0D0 + FID + FID
                  D1 = AK*DK
                  D2 = BK*CK
                  AD = MIN(D1,D2)
                  AK = 24.0D0 + 9.0D0*FID
                  BK = 30.0D0 - 9.0D0*FID
                  Z3R = DBLE(Z3)
                  Z3I = DIMAG(Z3)
                  DO 40 K = 1, 25
                     TRM1 = TRM1*DCMPLX(Z3R/D1,Z3I/D1)
                     S1 = S1 + TRM1
                     TRM2 = TRM2*DCMPLX(Z3R/D2,Z3I/D2)
                     S2 = S2 + TRM2
                     ATRM = ATRM*AZ3/AD
                     D1 = D1 + AK
                     D2 = D2 + BK
                     AD = MIN(D1,D2)
                     IF (ATRM.LT.TOL*AD) THEN
                        GO TO 60
                     ELSE
                        AK = AK + 18.0D0
                        BK = BK + 18.0D0
                     END IF
   40             CONTINUE
               END IF
   60          IF (ID.EQ.1) THEN
                  BI = S2*DCMPLX(C2,0.0D0)
                  IF (AZ.GT.TOL) BI = BI + Z*Z*S1*DCMPLX(C1/(1.0D0+FID),
     *                                0.0D0)
                  IF (KODE.NE.1) THEN
                     ZTA = Z*SQRT(Z)*DCMPLX(TTH,0.0D0)
                     AA = DBLE(ZTA)
                     AA = -ABS(AA)
                     BI = BI*DCMPLX(EXP(AA),0.0D0)
                  END IF
               ELSE
                  BI = S1*DCMPLX(C1,0.0D0) + Z*S2*DCMPLX(C2,0.0D0)
                  IF (KODE.NE.1) THEN
                     ZTA = Z*SQRT(Z)*DCMPLX(TTH,0.0D0)
                     AA = DBLE(ZTA)
                     AA = -ABS(AA)
                     BI = BI*DCMPLX(EXP(AA),0.0D0)
                  END IF
               END IF
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, DERIV has illegal value: DERIV = ''',A,
     *       '''')
99998 FORMAT (1X,'** On entry, SCALE has illegal value: SCALE = ''',A,
     *       '''')
99997 FORMAT (1X,'** No computation because abs(Z) =',1P,D13.5,' .GT.',
     *       D13.5)
99996 FORMAT (1X,'** Results lack precision because abs(Z) =',1P,D13.5,
     *       ' .GT.',D13.5)
99995 FORMAT (1X,'** No computation - algorithm termination condition ',
     *       'not met.')
99994 FORMAT (1X,'** No computation because real(Z) =',1P,D13.5,' is t',
     *       'oo large when SCALE = ''U''.')
      END
