      SUBROUTINE S17DLF(M,FNU,Z,N,SCALE,CY,NZ,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-781 (DEC 1989).
C
C     Original name: CBESH
C
C     PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
C
C     DESCRIPTION
C     ===========
C
C         ON SCALE='U', S17DLF COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
C         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
C         Z.NE.CMPLX(0.0E0,0.0E0) IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI.
C         ON SCALE='S', S17DLF COMPUTES THE SCALED HANKEL FUNCTIONS
C
C         CY(I)=H(M,FNU+J-1,Z)*EXP(-MM*Z*I)       MM=3-2M,      I**2=-1.
C
C         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER
C         AND LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN
C         THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0E0
C           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
C                    SCALE = 'U' OR SCALE = 'u' RETURNS
C                             CY(J)=H(M,FNU+J-1,Z),      J=1,...,N
C                          = 'S' OR SCALE = 's' RETURNS
C                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
C                                  J=1,...,N  ,  I**2=-1
C           M      - KIND OF HANKEL FUNCTION, M=1 OR 2
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(J)=H(M,FNU+J-1,Z)  OR
C                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
C                    DEPENDING ON SCALE, I**2=-1.
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
C                              J=1,...,NZ WHEN Y.GT.0.0 AND M=1 OR
C                              Y.LT.0.0 AND M=2. FOR THE COMPLMENTARY
C                              HALF PLANES, NZ STATES ONLY THE NUMBER
C                              OF UNDERFLOWS.
C           IERR    -ERROR FLAG
C                    IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
C                    IERR=1, INPUT ERROR   - NO COMPUTATION
C                    IERR=2, OVERFLOW      - NO COMPUTATION,
C                            CABS(Z) TOO SMALL
C                    IERR=3  OVERFLOW      - NO COMPUTATION,
C                            FNU+N-1 TOO LARGE
C                    IERR=4, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                    IERR=5, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                    IERR=6, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C     LONG DESCRIPTION
C     ================
C
C         THE COMPUTATION IS CARRIED OUT BY THE RELATION
C
C         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
C             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
C
C         FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
C         RIGHT HALF PLANE RE(Z).GE.0.0. THE K FUNCTION IS CONTINUED
C         TO THE LEFT HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
C         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
C         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
C         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
C         WHOLE Z PLANE FOR Z TO INFINITY.
C
C         FOR NEGATIVE ORDERS,THE FORMULAE
C
C               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
C               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
C                         I**2=-1
C
C         CAN BE USED.
C
C         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
C         MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
C         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
C         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
C         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
C         IERR=4 IS TRIGGERED WHERE UR=X02AJF()=UNIT ROUNDOFF. ALSO
C         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
C         LOST AND IERR=5. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
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
      PARAMETER         (SRNAME='S17DLF')
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  FNU
      INTEGER           IFAIL, M, N, NZ
      CHARACTER*1       SCALE
C     .. Array Arguments ..
      COMPLEX*16        CY(N)
C     .. Local Scalars ..
      COMPLEX*16        CSGN, ZN, ZT
      DOUBLE PRECISION  AA, ALIM, ALN, ARG, ASCLE, ATOL, AZ, BB, CPN,
     *                  DIG, ELIM, FMM, FN, FNUL, HPI, R1M5, RHPI, RL,
     *                  RTOL, SGN, SPN, TOL, UFL, XN, XX, YN, YY
      INTEGER           I, IERR, INU, INUH, IR, K, K1, K2, KODE, MM, MR,
     *                  NN, NREC, NUF, NW
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AHF, X02AJF
      INTEGER           P01ABF, X02BBF, X02BHF, X02BJF, X02BKF, X02BLF
      EXTERNAL          X02AHF, X02AJF, P01ABF, X02BBF, X02BHF, X02BJF,
     *                  X02BKF, X02BLF
C     .. External Subroutines ..
      EXTERNAL          S17DEV, S17DGX, S17DLY, S17DLZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, INT, LOG, LOG10,
     *                  MAX, MIN, MOD, DBLE, SIGN, SIN, SQRT
C     .. Data statements ..
C
      DATA              HPI/1.57079632679489662D0/
C     .. Executable Statements ..
      NZ = 0
      NREC = 0
      XX = DBLE(Z)
      YY = DIMAG(Z)
      IERR = 0
      IF (SCALE.EQ.'U' .OR. SCALE.EQ.'u') THEN
         KODE = 1
      ELSE IF (SCALE.EQ.'S' .OR. SCALE.EQ.'s') THEN
         KODE = 2
      ELSE
         KODE = -1
      END IF
      IF (XX.EQ.0.0D0 .AND. YY.EQ.0.0D0) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99999)
      ELSE IF (FNU.LT.0.0D0) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99998) FNU
      ELSE IF (KODE.EQ.-1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99997) SCALE
      ELSE IF (N.LT.1) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99996) N
      ELSE IF (M.LT.1 .OR. M.GT.2) THEN
         IERR = 1
         NREC = 1
         WRITE (REC,FMT=99995) M
      END IF
      IF (IERR.EQ.0) THEN
         NN = N
C        ---------------------------------------------------------------
C        SET PARAMETERS RELATED TO MACHINE CONSTANTS.
C        TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
C        ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
C        EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND
C        EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
C        UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
C        RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR
C        LARGE Z.
C        DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
C        FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE
C        FNU
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
         FNUL = 10.0D0 + 6.0D0*(DIG-3.0D0)
         RL = 1.2D0*DIG + 3.0D0
         FN = FNU + NN - 1
         MM = 3 - M - M
         FMM = MM
         ZN = Z*DCMPLX(0.0D0,-FMM)
         XN = DBLE(ZN)
         YN = DIMAG(ZN)
         AZ = ABS(Z)
C        ---------------------------------------------------------------
C        TEST FOR RANGE
C        ---------------------------------------------------------------
         AA = 0.5D0/TOL
         BB = X02BBF(1.0D0)*0.5D0
         AA = MIN(AA,BB,X02AHF(1.0D0))
         IF (AZ.LE.AA) THEN
            IF (FN.LE.AA) THEN
               AA = SQRT(AA)
               IF (AZ.GT.AA) THEN
                  IERR = 4
                  NREC = 1
                  WRITE (REC,FMT=99994) AZ, AA
               ELSE IF (FN.GT.AA) THEN
                  IERR = 4
                  NREC = 1
                  WRITE (REC,FMT=99993) FN, AA
               END IF
C              ---------------------------------------------------------
C              OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
C              ---------------------------------------------------------
               UFL = EXP(-ELIM)
               IF (AZ.GE.UFL) THEN
                  IF (FNU.GT.FNUL) THEN
C                    ---------------------------------------------------
C                    UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL
C                    ---------------------------------------------------
                     MR = 0
                     IF ((XN.LT.0.0D0) .OR. (XN.EQ.0.0D0 .AND. YN.LT.
     *                   0.0D0 .AND. M.EQ.2)) THEN
                        MR = -MM
                        IF (XN.EQ.0.0D0 .AND. YN.LT.0.0D0) ZN = -ZN
                     END IF
                     CALL S17DLY(ZN,FNU,KODE,MR,NN,CY,NW,TOL,ELIM,ALIM)
                     IF (NW.LT.0) THEN
                        GO TO 40
                     ELSE
                        NZ = NZ + NW
                     END IF
                  ELSE
                     IF (FN.GT.1.0D0) THEN
                        IF (FN.GT.2.0D0) THEN
                           CALL S17DEV(ZN,FNU,KODE,2,NN,CY,NUF,TOL,ELIM,
     *                                 ALIM)
                           IF (NUF.LT.0) THEN
                              GO TO 60
                           ELSE
                              NZ = NZ + NUF
                              NN = NN - NUF
C                             ------------------------------------------
C                             HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1
C                             ON RETURN FROM S17DEV
C                             IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C                             ------------------------------------------
                              IF (NN.EQ.0) THEN
                                 IF (XN.LT.0.0D0) THEN
                                    GO TO 60
                                 ELSE
                                    IFAIL = P01ABF(IFAIL,IERR,SRNAME,
     *                                      NREC,REC)
                                    RETURN
                                 END IF
                              END IF
                           END IF
                        ELSE IF (AZ.LE.TOL) THEN
                           ARG = 0.5D0*AZ
                           ALN = -FN*LOG(ARG)
                           IF (ALN.GT.ELIM) GO TO 60
                        END IF
                     END IF
                     IF ((XN.LT.0.0D0) .OR. (XN.EQ.0.0D0 .AND. YN.LT.
     *                   0.0D0 .AND. M.EQ.2)) THEN
C                       ------------------------------------------------
C                       LEFT HALF PLANE COMPUTATION
C                       ------------------------------------------------
                        MR = -MM
                        CALL S17DLZ(ZN,FNU,KODE,MR,NN,CY,NW,RL,FNUL,TOL,
     *                              ELIM,ALIM)
                        IF (NW.LT.0) THEN
                           GO TO 40
                        ELSE
                           NZ = NW
                        END IF
                     ELSE
C                       ------------------------------------------------
C                       RIGHT HALF PLANE COMPUTATION, XN.GE.0. .AND.
C                       (XN.NE.0. .OR.  YN.GE.0. .OR. M=1)
C                       ------------------------------------------------
                        CALL S17DGX(ZN,FNU,KODE,NN,CY,NZ,TOL,ELIM,ALIM)
                     END IF
                  END IF
C                 ------------------------------------------------------
C                 H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
C
C                 ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
C                 ------------------------------------------------------
                  SGN = SIGN(HPI,-FMM)
C                 ------------------------------------------------------
C                 CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF
C                 SIGNIFICANCE WHEN FNU IS LARGE
C                 ------------------------------------------------------
                  INU = INT(FNU)
                  INUH = INU/2
                  IR = INU - 2*INUH
                  ARG = (FNU-INU+IR)*SGN
                  RHPI = 1.0D0/SGN
                  CPN = RHPI*COS(ARG)
                  SPN = RHPI*SIN(ARG)
C                 ZN = CMPLX(-SPN,CPN)
                  CSGN = DCMPLX(-SPN,CPN)
C                 IF (MOD(INUH,2).EQ.1) ZN = -ZN
                  IF (MOD(INUH,2).EQ.1) CSGN = -CSGN
                  ZT = DCMPLX(0.0D0,-FMM)
                  RTOL = 1.0D0/TOL
                  ASCLE = UFL*RTOL
                  DO 20 I = 1, NN
C                    CY(I) = CY(I)*ZN
C                    ZN = ZN*ZT
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
                     CSGN = CSGN*ZT
   20             CONTINUE
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
   40             IF (NW.EQ.(-3)) THEN
                     NZ = 0
                     IERR = 5
                     NREC = 1
                     WRITE (REC,FMT=99988) AZ, AA
                     IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                     RETURN
                  ELSE IF (NW.NE.(-1)) THEN
                     NZ = 0
                     IERR = 6
                     NREC = 1
                     WRITE (REC,FMT=99992)
                     IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                     RETURN
                  END IF
   60             IERR = 3
                  NZ = 0
                  NREC = 1
                  WRITE (REC,FMT=99991) FN
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               ELSE
                  IERR = 2
                  NZ = 0
                  NREC = 1
                  WRITE (REC,FMT=99990) AZ, UFL
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               END IF
            ELSE
               NZ = 0
               IERR = 5
               NREC = 1
               WRITE (REC,FMT=99989) FN, AA
            END IF
         ELSE
            NZ = 0
            IERR = 5
            NREC = 1
            WRITE (REC,FMT=99988) AZ, AA
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, Z = (0.0,0.0)')
99998 FORMAT (1X,'** On entry, FNU .LT. 0: FNU = ',D13.5)
99997 FORMAT (1X,'** On entry, SCALE has an illegal value: SCALE = ''',
     *       A,'''')
99996 FORMAT (1X,'** On entry, N .LE. 0: N = ',I16)
99995 FORMAT (1X,'** On entry, M has illegal value: M = ',I16)
99994 FORMAT (1X,'** Results lack precision because abs(Z) =',1P,D13.5,
     *       ' .GT.',D13.5)
99993 FORMAT (1X,'** Results lack precision, FNU+N-1 =',1P,D13.5,
     *       ' .GT.',D13.5)
99992 FORMAT (1X,'** No computation - algorithm termination condition ',
     *       'not met.')
99991 FORMAT (1X,'** No computation because FNU+N-1 =',1P,D13.5,' is t',
     *       'oo large.')
99990 FORMAT (1X,'** No computation because abs(Z) =',1P,D13.5,' .LT. ',
     *       D13.5)
99989 FORMAT (1X,'** No computation because FNU+N-1 =',1P,D13.5,' .GT.',
     *       D13.5)
99988 FORMAT (1X,'** No computation because abs(Z) =',1P,D13.5,' .GT.',
     *       D13.5)
      END
