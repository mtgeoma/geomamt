      SUBROUTINE S18DCF(FNU,Z,N,SCALE,CY,NZ,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-784 (DEC 1989).
C
C     Original name: CBESK
C
C     PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
C
C     DESCRIPTION
C     ===========
C
C         ON SCALE='U', S18DCF COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
C         BESSEL FUNCTIONS CY(J)=K(FNU+J-1,Z) FOR REAL, NONNEGATIVE
C         ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z.NE.CMPLX(0.0,0.0)
C         IN THE CUT PLANE -PI.LT.ARG(Z).LE.PI. ON SCALE='S', S18DCF
C         RETURNS THE SCALED K FUNCTIONS,
C
C         CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,
C
C         WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
C         RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
C         NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
C         FUNCTIONS (REF. 1).
C
C         INPUT
C           Z      - Z=CMPLX(X,Y),Z.NE.CMPLX(0.,0.),-PI.LT.ARG(Z).LE.PI
C           FNU    - ORDER OF INITIAL K FUNCTION, FNU.GE.0.0E0
C           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
C           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
C                    SCALE = 'U' OR SCALE = 'u' RETURNS
C                             CY(I)=K(FNU+I-1,Z), I=1,...,N
C                    SCALE = 'S' OR SCALE = 's' RETURNS
C                             CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C
C         OUTPUT
C           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
C                    VALUES FOR THE SEQUENCE
C                    CY(I)=K(FNU+I-1,Z), I=1,...,N OR
C                    CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
C                    DEPENDING ON SCALE
C           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
C                    NZ= 0   , NORMAL RETURN
C                    NZ.GT.0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
C                              DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
C                              I=1,...,N WHEN X.GE.0.0. WHEN X.LT.0.0
C                              NZ STATES ONLY THE NUMBER OF UNDERFLOWS
C                              IN THE SEQUENCE.
C           IFAIL  - ERROR FLAG
C                   IFAIL=0, NORMAL RETURN - COMPUTATION COMPLETED
C                   IFAIL=1, INPUT ERROR   - NO COMPUTATION
C                   IFAIL=2, OVERFLOW      - NO COMPUTATION, CABS(Z) IS
C                            TOO SMALL
C                   IFAIL=3, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
C                            TOO LARGE
C                   IFAIL=4, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
C                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
C                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
C                            ACCURACY
C                   IFAIL=5, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
C                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
C                            CANCE BY ARGUMENT REDUCTION
C                   IFAIL=6, ERROR              - NO COMPUTATION,
C                            ALGORITHM TERMINATION CONDITION NOT MET
C
C     LONG DESCRIPTION
C     ================
C
C         EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS
C         DNU AND DNU+1.0 IN THE RIGHT HALF PLANE X.GE.0.0. FORWARD
C         RECURRENCE GENERATES HIGHER ORDERS. K IS CONTINUED TO THE LEFT
C         HALF PLANE BY THE RELATION
C
C         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
C         MP=MR*PI*I, MR=+1 OR -1, RE(Z).GT.0, I**2=-1
C
C         WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.
C
C         FOR LARGE ORDERS, FNU.GT.FNUL, THE K FUNCTION IS COMPUTED
C         BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.
C
C         FOR NEGATIVE ORDERS, THE FORMULA
C
C                       K(-FNU,Z) = K(FNU,Z)
C
C         CAN BE USED.
C
C         S18DCF ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS
C         AVAILABLE.
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
C                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983.
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
      PARAMETER         (SRNAME='S18DCF')
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  FNU
      INTEGER           IFAIL, N, NZ
      CHARACTER*1       SCALE
C     .. Array Arguments ..
      COMPLEX*16        CY(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, ALIM, ALN, ARG, AZ, BB, DIG, ELIM, FN, FNUL,
     *                  R1M5, RL, TOL, UFL, XX, YY
      INTEGER           IERR, K, K1, K2, KODE, MR, NN, NREC, NUF, NW
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
      INTRINSIC         ABS, DIMAG, EXP, LOG, LOG10, MAX, MIN, DBLE,
     *                  SQRT
C     .. Executable Statements ..
      IERR = 0
      NREC = 0
      NZ = 0
      XX = DBLE(Z)
      YY = DIMAG(Z)
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
C        FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR
C        LARGE FNU
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
         AZ = ABS(Z)
         FN = FNU + NN - 1
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
                  WRITE (REC,FMT=99995) AZ, AA
               ELSE IF (FN.GT.AA) THEN
                  IERR = 4
                  NREC = 1
                  WRITE (REC,FMT=99994) FN, AA
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
                     IF (XX.LT.0.0D0) THEN
                        MR = 1
                        IF (YY.LT.0.0D0) MR = -1
                     END IF
                     CALL S17DLY(Z,FNU,KODE,MR,NN,CY,NW,TOL,ELIM,ALIM)
                     IF (NW.GE.0) THEN
                        NZ = NZ + NW
                        IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                        RETURN
                     END IF
                  ELSE
                     IF (FN.GT.1.0D0) THEN
                        IF (FN.GT.2.0D0) THEN
                           CALL S17DEV(Z,FNU,KODE,2,NN,CY,NUF,TOL,ELIM,
     *                                 ALIM)
                           IF (NUF.LT.0) THEN
                              GO TO 20
                           ELSE
                              NZ = NZ + NUF
                              NN = NN - NUF
C                             ------------------------------------------
C                             HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1
C                             ON RETURN FROM S17DEV
C                             IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
C                             ------------------------------------------
                              IF (NN.EQ.0) THEN
                                 IF (XX.LT.0.0D0) THEN
                                    GO TO 20
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
                           IF (ALN.GT.ELIM) GO TO 20
                        END IF
                     END IF
                     IF (XX.LT.0.0D0) THEN
C                       ------------------------------------------------
C                       LEFT HALF PLANE COMPUTATION
C                       PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2.
C                       ------------------------------------------------
                        IF (NZ.NE.0) THEN
                           GO TO 20
                        ELSE
                           MR = 1
                           IF (YY.LT.0.0D0) MR = -1
                           CALL S17DLZ(Z,FNU,KODE,MR,NN,CY,NW,RL,FNUL,
     *                                 TOL,ELIM,ALIM)
                           IF (NW.GE.0) THEN
                              NZ = NW
                              IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                              RETURN
                           END IF
                        END IF
                     ELSE
C                       ------------------------------------------------
C                       RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0.
C                       ------------------------------------------------
                        CALL S17DGX(Z,FNU,KODE,NN,CY,NW,TOL,ELIM,ALIM)
                        IF (NW.GE.0) THEN
                           NZ = NW
                           IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                           RETURN
                        END IF
                     END IF
                  END IF
                  IF (NW.EQ.(-3)) THEN
                     NZ = 0
                     IERR = 5
                     NREC = 1
                     WRITE (REC,FMT=99989) AZ, AA
                     IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                     RETURN
                  ELSE IF (NW.NE.(-1)) THEN
                     NZ = 0
                     IERR = 6
                     NREC = 1
                     WRITE (REC,FMT=99993)
                     IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                     RETURN
                  END IF
   20             IERR = 3
                  NZ = 0
                  NREC = 1
                  WRITE (REC,FMT=99992) FN
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               ELSE
                  IERR = 2
                  NZ = 0
                  NREC = 1
                  WRITE (REC,FMT=99991) AZ, UFL
                  IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
                  RETURN
               END IF
            ELSE
               NZ = 0
               IERR = 5
               NREC = 1
               WRITE (REC,FMT=99990) FN, AA
            END IF
         ELSE
            NZ = 0
            IERR = 5
            NREC = 1
            WRITE (REC,FMT=99989) AZ, AA
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
99995 FORMAT (1X,'** Results lack precision because abs(Z) =',1P,D13.5,
     *       ' .GT.',D13.5)
99994 FORMAT (1X,'** Results lack precision because FNU+N-1 =',1P,D13.5,
     *       ' .GT.',D13.5)
99993 FORMAT (1X,'** No computation - algorithm termination condition ',
     *       'not met.')
99992 FORMAT (1X,'** No computation because FNU+N-1 =',1P,D13.5,' is t',
     *       'oo large.')
99991 FORMAT (1X,'** No computation because abs(Z) =',1P,D13.5,' .LT.',
     *       D13.5)
99990 FORMAT (1X,'** No computation because FNU+N-1 =',1P,D13.5,' .GT.',
     *       D13.5)
99989 FORMAT (1X,'** No computation because abs(Z) =',1P,D13.5,' .GT.',
     *       D13.5)
      END
