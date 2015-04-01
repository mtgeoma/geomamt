      SUBROUTINE D01AJY(N,EPSTAB,RESULT,ABSERR,RES3LA,NRES)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QEXT.
C     ................................................................
C
C           PURPOSE
C              THE ROUTINE DETERMINES THE LIMIT OF A GIVEN SEQUENCE OF
C              APPROXIMATIONS, BY MEANS OF THE EPSILON ALGORITHM
C              OF P. WYNN.
C              AN ESTIMATE OF THE ABSOLUTE ERROR IS ALSO GIVEN.
C              THE CONDENSED EPSILON TABLE IS COMPUTED. ONLY THOSE
C              ELEMENTS NEEDED FOR THE COMPUTATION OF THE NEXT DIAGONAL
C              ARE PRESERVED.
C
C           PARAMETERS
C              N      - INTEGER
C                       EPSTAB(N) CONTAINS THE NEW ELEMENT IN THE
C                       FIRST COLUMN OF THE EPSILON TABLE.
C
C              EPSTAB - REAL
C                       VECTOR OF DIMENSION 52 CONTAINING THE ELEMENTS
C                       OF THE TWO LOWER DIAGONALS OF THE TRIANGULAR
C                       EPSILON TABLE
C                       THE ELEMENTS ARE NUMBERED STARTING AT THE
C                       RIGHT-HAND CORNER OF THE TRIANGLE.
C
C              RESULT - REAL
C                       RESULTING APPROXIMATION TO THE INTEGRAL
C
C              ABSERR - REAL
C                       ESTIMATE OF THE ABSOLUTE ERROR COMPUTED FROM
C                       RESULT AND THE 3 PREVIOUS RESULTS
C
C              RES3LA - REAL
C                       VECTOR OF DIMENSION 3 CONTAINING THE LAST 3
C                       RESULTS
C
C              NRES   - INTEGER
C                       NUMBER OF CALLS TO THE ROUTINE
C                       (SHOULD BE ZERO AT FIRST CALL)
C
C     ..................................................................
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ABSERR, RESULT
      INTEGER           N, NRES
C     .. Array Arguments ..
      DOUBLE PRECISION  EPSTAB(52), RES3LA(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  DELTA1, DELTA2, DELTA3, E0, E1, E1ABS, E2, E3,
     *                  EPMACH, EPSINF, ERR1, ERR2, ERR3, ERROR, OFLOW,
     *                  RES, SS, TOL1, TOL2, TOL3, UFLOW
      INTEGER           I, IB, IB2, IE, INDX, K1, K2, K3, LIMEXP,
     *                  NEWELM, NUM
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           E0     - THE 4 ELEMENTS ON WHICH THE
C           E1       COMPUTATION OF A NEW ELEMENT IN
C           E2       THE EPSILON TABLE IS BASED
C           E3                 E0
C                        E3    E1    NEW
C                              E2
C           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
C                    DIAGONAL
C           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
C           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
C                    OF ERROR
C
C           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON TABLE
C           CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER DIAGONAL
C           OF THE EPSILON TABLE IS DELETED.
C
      EPMACH = X02AJF()
      UFLOW = X02AMF()
      OFLOW = 1.0D+00/UFLOW
      NRES = NRES + 1
      ABSERR = OFLOW
      RESULT = EPSTAB(N)
      IF (N.LT.3) GO TO 200
      LIMEXP = 50
      EPSTAB(N+2) = EPSTAB(N)
      NEWELM = (N-1)/2
      EPSTAB(N) = OFLOW
      NUM = N
      K1 = N
      DO 80 I = 1, NEWELM
         K2 = K1 - 1
         K3 = K1 - 2
         RES = EPSTAB(K1+2)
         E0 = EPSTAB(K3)
         E1 = EPSTAB(K2)
         E2 = RES
         E1ABS = ABS(E1)
         DELTA2 = E2 - E1
         ERR2 = ABS(DELTA2)
         TOL2 = MAX(ABS(E2),E1ABS)*EPMACH
         DELTA3 = E1 - E0
         ERR3 = ABS(DELTA3)
         TOL3 = MAX(E1ABS,ABS(E0))*EPMACH
         IF (ERR2.GT.TOL2 .OR. ERR3.GT.TOL3) GO TO 20
C
C           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE ACCURACY,
C           CONVERGENCE IS ASSUMED.
C           RESULT = E2
C           ABSERR = ABS(E1-E0)+ABS(E2-E1)
C
         RESULT = RES
         ABSERR = ERR2 + ERR3
C        ***JUMP OUT OF DO-LOOP
         GO TO 200
   20    E3 = EPSTAB(K1)
         EPSTAB(K1) = E1
         DELTA1 = E1 - E3
         ERR1 = ABS(DELTA1)
         TOL1 = MAX(E1ABS,ABS(E3))*EPMACH
C
C           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT A PART
C           OF THE TABLE BY ADJUSTING THE VALUE OF N
C
         IF (ERR1.LE.TOL1 .OR. ERR2.LE.TOL2 .OR. ERR3.LE.TOL3)
     *       GO TO 40
         SS = 1.0D+00/DELTA1 + 1.0D+00/DELTA2 - 1.0D+00/DELTA3
         EPSINF = ABS(SS*E1)
C
C           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
C           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
C           OF N.
C
         IF (EPSINF.GT.1.0D-04) GO TO 60
   40    N = I + I - 1
C        ***JUMP OUT OF DO-LOOP
         GO TO 100
C
C           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST THE VALUE OF
C           RESULT.
C
   60    RES = E1 + 1.0D+00/SS
         EPSTAB(K1) = RES
         K1 = K1 - 2
         ERROR = ERR2 + ABS(RES-E2) + ERR3
         IF (ERROR.GT.ABSERR) GO TO 80
         ABSERR = ERROR
         RESULT = RES
   80 CONTINUE
C
C           SHIFT THE TABLE.
C
  100 IF (N.EQ.LIMEXP) N = 2*(LIMEXP/2) - 1
      IB = 1
      IF ((NUM/2)*2.EQ.NUM) IB = 2
      IE = NEWELM + 1
      DO 120 I = 1, IE
         IB2 = IB + 2
         EPSTAB(IB) = EPSTAB(IB2)
         IB = IB2
  120 CONTINUE
      IF (NUM.EQ.N) GO TO 160
      INDX = NUM - N + 1
      DO 140 I = 1, N
         EPSTAB(I) = EPSTAB(INDX)
         INDX = INDX + 1
  140 CONTINUE
  160 IF (NRES.GE.4) GO TO 180
      RES3LA(NRES) = RESULT
      ABSERR = OFLOW
      GO TO 200
C
C           COMPUTE ERROR ESTIMATE
C
  180 ABSERR = ABS(RESULT-RES3LA(3)) + ABS(RESULT-RES3LA(2)) +
     *         ABS(RESULT-RES3LA(1))
      RES3LA(1) = RES3LA(2)
      RES3LA(2) = RES3LA(3)
      RES3LA(3) = RESULT
  200 ABSERR = MAX(ABSERR,0.5D+00*EPMACH*ABS(RESULT))
      RETURN
      END
