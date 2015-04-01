      SUBROUTINE D03UBF(N1,N2,N3,N1M,N2M,A,B,C,D,E,F,G,APARAM,IT,R,
     *                  WRKSP1,WRKSP2,WRKSP3,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10A REVISED. IER-387 (OCT 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     **************************************************************
C
C            DAVID A. H.  JACOBS         22/5/79
C     D03UBF PERFORMS ONE ITERATION OF THE STRONGLY IMPLICIT
C     PROCEDURE AT EACH CALL TO CALCULATE THE SUCCESSIVE
C     APPROXIMATE CORRECTIONS TO THE SOLUTION OF A SYSTEM
C     OF SIMULTANEOUS ALGEBRAIC EQUATIONS FOR WHICH THE
C     ITERATIVE UPDATE MATRIX IS OF THE SEVEN POINT
C     MOLECULE FORM ON A TOPOLOGICALLY THREE-DIMENSIONAL
C     RECTANGULAR MESH.
C
C     INPUTS
C     ------
C
C     N1       NUMBER OF NODES IN THE FIRST COORDINATE DIRECTION.
C
C     N2       NUMBER OF NODES IN THE SECOND COORDINATE DIRECTION.
C
C     N3       NUMBER OF NODES IN THE THIRD COORDINATE DIRECTION.
C
C     N1M      FIRST DIMENSION OF ALL THREE DIMENSIONAL ARRAYS.
C
C     N2M      SECOND DIMENSION OF ALL THREE DIMENSIONAL ARRAYS.
C
C     A,B,C,   ARRAYS ALL OF DIMENSION (N1M,N2M,N3) STORING
C     D,E,F,G  THE COEFFICIENTS OF THE ITERATIVE UPDATE
C              EQUATIONS  -
C
C       A(I,J,K)*S(I,J,K-1)+B(I,J,K)*S(I,J-1,K)+C(I,J,K)*S(I-1,J,K)+
C       D(I,J,K)*S(I,J,K)+E(I,J,K)*S(I+1,J,K)+F(I,J,K)*S(I,J+1,K)+
C       G(I,J,K)*S(I,J,K+1)=R(I,J,K)
C
C              WITH I=1,2,...,N1 , J=1,2,...,N2 AND K=1,2,...,N3
C              AND WHERE S(I,J,K) IS THE ARRAY WHOSE
C              APPROXIMATE VALUES ARE SOUGHT (AND WHICH
C              OVERWRITE THE RESIDUALS R(I,J,K) STORED ON
C              INPUT IN THE ARRAY R). ANY VALUES OF S
C              OUTSIDE THE  (1-N1,1-N2,1-N3) ARRAY ARE
C              TAKEN AS ZERO.
C
C     R        IS INPUT AS THE ARRAY OF RESIDUALS,
C              CORRESPONDING TO R ABOVE, CHANGED ON OUTPUT
C              TO THE VALUES  OF THE UPDATE ARRAY
C              CORRESPONDING TO S ABOVE, DIMENSIONED
C              R(N1M,N2M,N3)
C
C     APARAM   IS AN ITERATION ACCELERATION PARAMETER FACTOR
C              TYPICALLY SET TO 1.0. IF CONVERGENCE IS
C              SLOW, IT CAN BE DECREASED. IF DIVERGENCE IS
C              OBTAINED, IT SHOULD BE INCREASED. IN EITHER
C              CASE IT MUST NOT GO OUTSIDE THE BOUNDS
C              PRESCRIBED. (SEE IFAIL PARAMETER DETAILS.)
C
C     IT        IS THE ITERATION COUNTER SET AND
C              INCREMENTED BY THE CALLING ROUTINE, IT IS
C              USED TO DETERMINE THE APPROPRIATE ITERATION
C              PARAMETERS.
C
C     WRKSP1   IS A WORKSPACE ARRAY OF DIMENSION (N1M,N2M,N3)
C
C     WRKSP2    ... DITTO ...
C
C     WRKSP3    ... DITTO ...
C
C     IFAIL    IS AN ERROR PARAMETER INDICATOR SET BY THE USER TO
C              INDICATE THE TYPE OF FAILURE IF AN ERROR IS
C              ENCOUNTERED.
C
C     ALL ABOVE PASSED AS ARGUMENTS
C
C     PROCESS
C     -------
C
C     SET ERROR PARAMETER
C     CHECK INPUT INTEGERS
C     SET FREQUENTLY REQUIRED INTEGER VARIABLES
C     COMMENCEMENT OF CALCULATIONAL PROCEDURE
C     SET ODD/EVEN COUNTERS   KS=1  FOR ODD ITERATIONS
C                            KS=2  FOR EVEN ITERATIONS
C     DETERMINE THE NUMBER OF THE ACCELERATION PARAMETER TO
C     BE USED (THE SAME PARAMETER IS USED FOR 2 ITERATIONS,
C     THERE ARE 9 PARAMETERS IN ALL).
C     FIRST CALCULATE THE TERM, ALM, IN THE LARGEST PARAMETER
C     CHECK THE VALUES OF ALM
C     DETERMINE THE ITERATION PARAMETER
C     START OF APPROXIMATE LU FACTORIZATION DETERMINING THE
C     VALUES OF THE ELEMENTS SA,SB,SC IN THE LOWER
C     TRIANGULAR MATRIX L AND WRKSP1, WRKSP2 AND WRKSP3 IN
C     THE UPPER TRIANGULAR MATRIX U. PROGRESS FORWARDS
C     FIRST  AND INVERT THE LOWER TRIANGULAR MATRIX L AS
C     ONE PROCEEDS SO THAT THE ELEMENTS SA,SB,SC DO NOT
C     HAVE TO BE STORED IN ARRAYS.THE ELEMENTS WRKSP1,
C     WRKSP2 AND WRKSP3 HAVVE TO BE STORED FOR THE
C     SUBSEQUENT  INVERSION OF THE MATRIX U BY BACK
C     SUBSTITUTION.
C     DO KL=1,N3
C       (KL IS ONLY A COUNTER IN THE THIRD COORDINATE DIRECTION
C        K IS THE INDEX OF THE THIRD COORDINATE AND ON ALTERNATE
C        ITERATIONS SCANS FIRST INCREASING AND THEN DECREASING).
C       DO JL=1,N2
C          (JL IS ONLY A COUNTER IN THE SECOND COORDINATE DIRECTION
C           J IS THE INDEX OF THE SECOND COORDINATE AND ON ALTERNATE
C           ITERATIONS SCANS FIRST INCREASING AND THEN DECREASING).
C          DO I=1,N1
C             (I IS THE INDEX IN THE FIRST COORDINATE DIRECTION)
C             STORE THE VALUE OF THE CENTRAL COEFFICIENT
C             DETERMINE THE TYPE OF THE DIFFERENCE EQUATION
C                FOR A SEVEN POINT MOLECULE EQUATION -
C                   SET VALUES OF ADJACENT ARRAY ELEMENTS
C                    DEPENDENT ON NODAL POSITION
C
C                   CALCULATE THE OFF-DIAGONAL ELEMENTS OF
C                   L, NAMELY SA,SB,SC
C                   CALCULATE OFTEN USED EXPRESIONS
C                   CALCULATE THE VALUE OF THE DIAGONAL ELEMENT OF L
C                   CALCULATE AND STORE THE ELEMENTS OF U,
C                    NAMELY WRKSP1, WRKSP2 AND WRKSP3
C
C                   CALCULATE THE ELEMENTS OF L**(-1) *
C                   RESIDUAL ARRAY
C                FOR THE EXPLICIT EQUATION DO THE SAME
C     END OF SCAN FOR THE FORWARD ELIMINATION
C     BACK SUBSTITUTION TO MULTIPLY BY  U**(-1)
C     DO KL=1,N3
C       K RUNS BACKWARDS AND FORWARDS ALTERNATELY
C       DO JL=1,N2
C          J RUNS BACKWARDS AND FORWARDS ALTERNATELY
C          DO I=N1,N1-1,...,2,1
C             R IS INITIALLY L**(-1) * RESIDUAL, IT BECOMES
C              (U**(-1) * (L**(-1) * RESIDUAL) = CHANGE
C     END OF SCAN FOR THE BACK SUBSTITUTION
C     RETURN
C
C     OUTPUTS
C     -------
C
C     R         ARRAY STORING THE APPROXIMATE SOLUTION TO
C               THE SYSTEM OF EQUATIONS PROVIDED, AFTER
C               THE ONE ITERATION. NOTE THAT IT HAS
C               OVERWRITTEN THE INPUT RESIDUAL.
C
C     IFAIL     ERROR INDICATOR
C          = 0  CORRECT RETURN
C          = 1  EITHER N1.LE.1 OR N2.LE.1 OR N3.LE.1
C          = 2  N1M IS LESS THAN N1 OR N2M IS LESS THAN N2
C          = 3  APARAM.LE.0.0
C          = 4  APARAM.GT.((N1-1)**2+(N2-1)**2+(N3-1)**2)/3.
C
C     BOTH PASSED AS ARGUMENTS
C
C     ROUTINES USED
C     -------------
C
C     P01AAF     ERROR HANDLING ROUTINE
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03UBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  APARAM
      INTEGER           IFAIL, IT, N1, N1M, N2, N2M, N3
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N1M,N2M,N3), B(N1M,N2M,N3), C(N1M,N2M,N3),
     *                  D(N1M,N2M,N3), E(N1M,N2M,N3), F(N1M,N2M,N3),
     *                  G(N1M,N2M,N3), R(N1M,N2M,N3),
     *                  WRKSP1(N1M,N2M,N3), WRKSP2(N1M,N2M,N3),
     *                  WRKSP3(N1M,N2M,N3)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALM, ALPHA, DS, EA, EG, EW, RIBJK, RIJBK, RIJKB,
     *                  SA, SB, SC, SD, SEIBJK, SEIJBK, SEIJKB, SFIBJK,
     *                  SFIJBK, SFIJKB, SGIBJK, SGIJBK, SGIJKB, XKS1,
     *                  XKS2
      INTEGER           I, IB, IERROR, IL, IS, J, JB, JL, K, KB, KL, KS,
     *                  KS1, KS2, N1M1, N1P1, N2M1, N2P1, N3M1, N3P1
C     .. Local Arrays ..
      DOUBLE PRECISION  ALP(9)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE
C     .. Data statements ..
      DATA              ALP(1), ALP(2), ALP(3), ALP(4), ALP(5), ALP(6),
     *                  ALP(7), ALP(8), ALP(9)/1.0D0, 0.625D0, 0.25D0,
     *                  0.875D0, 0.5D0, 0.125D0, 0.75D0, 0.375D0, 0.0D0/
C     .. Executable Statements ..
C
C     SET ERROR PARAMETER
C
      IERROR = 0
C
C     CHECK INPUT INTEGERS
C
      IF (N1.LE.1) GO TO 300
      IF (N2.LE.1) GO TO 300
      IF (N3.LE.1) GO TO 300
C
      IF (N1M.LT.N1) GO TO 320
      IF (N2M.LT.N2) GO TO 320
C
C     SET FREQUENTLY REQUIRED INTEGER VARIABLES
C
      N1M1 = N1 - 1
      N2M1 = N2 - 1
      N3M1 = N3 - 1
      N1P1 = N1 + 1
      N2P1 = N2 + 1
      N3P1 = N3 + 1
C
C     COMMENCEMENT OF CALCULATIONAL PROCEDURE
C     ---------------------------------------
C
C
C     SET ODD/EVEN COUNTER    KS=1  FOR ODD ITERATIONS
C                            KS=2  FOR EVEN ITERATIONS
C
      KS = MOD(IT-1,2) + 1
      IF (KS.EQ.0) KS = 2
C
C     SET FREQUENTLY USED ODD/EVEN PARAMETER DERIVATIVES
C
      KS1 = 2 - KS
      XKS1 = DBLE(KS1)
      KS2 = KS - 1
      XKS2 = DBLE(KS2)
C
C     DETERMINE THE NUMBER OF THE ACCELERATION PARAMETER TO
C     BE USED, THE SAME PARAMETER IS USED FOR 2 ITERATIONS,
C     THERE ARE 9 PARAMETERS IN ALL.
      IS = MOD(IT-1,18)
      IF (IS.LT.0) IS = IS + 18
      IS = IS/2 + 1
C
C       CALCULATE THE ITERATION ACCELERATION PARAMETER
C
C       (1) CALCULATE THE TERM IN THE LARGEST PARAMETER
C
      ALM = 3.D0*APARAM/(DBLE(N1M1*N1M1+N2M1*N2M1+N3M1*N3M1))
C
C       (2) CHECK THE VALUES OF ALM
C
      IF (APARAM.LE.0.0D0) GO TO 340
C
      IF (ALM.GT.1.0D0) GO TO 360
C
C       (3) THEN DETERMINE THE ITERATION PARAMETER
C
      ALPHA = 1.D0 - ALM**ALP(IS)
C
C     START OF APPROXIMATE LU FACTORIZATION DETERMINING THE
C     VALUES OF THE ELEMENTS SA,SB,SC,SD IN THE LOWER
C     TRIANGULAR MATRIX L AND WRKSP1, WRKSP2 AND WRKSP3 IN
C     THE UPPER TRIANGULAR MATRIX U. PROGRESS FORWARDS FIRST
C     AND INVERT THE LOWER TRIANGULAR MATRIX L AS ONE
C     PROCEEDS SO THAT  THE ELEMENTS SA,SB,SC,SD DO NOT HAVE
C     TO BE STORED IN ARRAYS.THE ELEMENTS WRKSP1, WRKSP2 AND
C     WRKSP3 HAVE TO BE STORED FOR THE SUBSEQUENT INVERSION
C     OF THE MATRIX U BY BACK SUBSTITUTION.
      DO 220 KL = 1, N3
C
C        KL IS ONLY A COUNTER IN THE THIRD COORDINATE DIRECTION
C        K IS THE INDEX OF THE THIRD COORDINATE AND ON ALTERNATE
C        ITERATIONS SCANS FIRST INCREASING THE DECREASING
C
         K = KS1*KL + KS2*(N3P1-KL)
         KB = K - KS1 + KS2
C
         DO 200 JL = 1, N2
C
C           JL IS ONLY A COUNTER IN THE SECOND COORDINATE DIRECTION
C           J IS THE INDEX OF THE THIRD COORDINATE AND ON ALTERNATE
C           ITERATIONS SCANS FIRST INCREASING THEN DECREASING
C
            J = KS1*JL + KS2*(N2P1-JL)
            JB = J - KS1 + KS2
C
            DO 180 I = 1, N1
C
C              I IS THE INDEX IN THE FIRST COORDINATE DIRECTION
C
               IB = I - 1
C
C              STORE THE VALUE OF THE CENTRAL COEFFICIENT
C
               DS = D(I,J,K)
C
C              DETERMINE THE TYPE OF THE DIFFERENCE EQUATION
C
               IF (DS.EQ.0.0D0) GO TO 140
C
C                SEVEN POINT MOLECULE EQUATION
C
C                SET VALUES OF ADJACENT ARRAY ELEMENTS
C                 DEPENDENT ON NODAL POSITION
C
               IF ((KB.EQ.0) .OR. (KB.EQ.N3P1)) GO TO 20
C
C                   (1A) IF NOT ON BOUNDARY SO THAT  KB NE 0 OR N3P1
C
               SEIJKB = WRKSP1(I,J,KB)
               SFIJKB = WRKSP2(I,J,KB)
               SGIJKB = WRKSP3(I,J,KB)
               RIJKB = R(I,J,KB)
               GO TO 40
C
   20          CONTINUE
C
C                   (1B) BOUNDARY VALUES EXTERIOR TO ARRAY
C                   SET TO ZERO
               SEIJKB = 0.0D0
               SFIJKB = 0.0D0
               SGIJKB = 0.0D0
               RIJKB = 0.0D0
   40          CONTINUE
C
               IF ((JB.EQ.0) .OR. (JB.EQ.N2P1)) GO TO 60
C
C                   (2A) IF NOT ON  BOUNDARY  SO THAT
C                   JB.NE.0 OR N2P1
               SEIJBK = WRKSP1(I,JB,K)
               SFIJBK = WRKSP2(I,JB,K)
               SGIJBK = WRKSP3(I,JB,K)
               RIJBK = R(I,JB,K)
               GO TO 80
C
   60          CONTINUE
C
C                   (2B) THE VALUES EXTERIOR TO THE ARRAYS
C                   ARE SET TO ZERO
               SEIJBK = 0.0D0
               SFIJBK = 0.0D0
               SGIJBK = 0.0D0
               RIJBK = 0.0D0
   80          CONTINUE
C
               IF (I.EQ.1) GO TO 100
C
C                   (3A) IF NOT ON I=1 SO THAT  IB.NE.0
C
               SEIBJK = WRKSP1(IB,J,K)
               SFIBJK = WRKSP2(IB,J,K)
               SGIBJK = WRKSP3(IB,J,K)
               RIBJK = R(IB,J,K)
               GO TO 120
C
  100          CONTINUE
C
C                   (3B) ON  I=1  THE VALUES EXTERIOR TO
C                   THE ARRAY ARE SET TO 0
               SEIBJK = 0.0D0
               SFIBJK = 0.0D0
               SGIBJK = 0.0D0
               RIBJK = 0.0D0
  120          CONTINUE
C
C                CALCULATE THE ELEMENTS OF THE LOWER
C                TRIANGULAR MATRIX
               SA = (XKS1*A(I,J,K)+XKS2*G(I,J,K))
     *              /(1.D0+ALPHA*(SEIJKB+SFIJKB))
               SB = (XKS1*B(I,J,K)+XKS2*F(I,J,K))
     *              /(1.D0+ALPHA*(SGIJBK+SEIJBK))
               SC = C(I,J,K)/(1.D0+ALPHA*(SGIBJK+SFIBJK))
C
C                CALCULATE OFTEN USED EXPRESIONS
C
C
               EA = SA*SEIJKB + SB*SEIJBK
               EG = SC*SFIBJK + SA*SFIJKB
               EW = SB*SGIJBK + SC*SGIBJK
C
C                CALCULATE THE VALUE OF THE DIAGONAL ELEMENT OF L
C
               SD = 1.D0/(DS+ALPHA*(EA+EG+EW)
     *              -(SC*SEIBJK+SB*SFIJBK+SA*SGIJKB))
C
C                CALCULATE AND STORE THE ELEMENTS OF THE
C                UPPER TRIANGULAR MATRIX
               WRKSP1(I,J,K) = (E(I,J,K)-ALPHA*EA)*SD
               WRKSP2(I,J,K) = (XKS1*F(I,J,K)+XKS2*B(I,J,K)-ALPHA*EG)*SD
               WRKSP3(I,J,K) = (XKS1*G(I,J,K)+XKS2*A(I,J,K)-ALPHA*EW)*SD
C
C                CALCULATION OF THE ELEMENTS OF THE INVERSE
C                OF L * RES VECTOR
               R(I,J,K) = SD*(R(I,J,K)-SA*RIJKB-SB*RIJBK-SC*RIBJK)
C
               GO TO 160
C
  140          CONTINUE
C
C                CALCULATION FOR THE EXPLICIT EQUATION
C
               WRKSP1(I,J,K) = 0.0D0
               WRKSP2(I,J,K) = 0.0D0
               WRKSP3(I,J,K) = 0.0D0
C
  160          CONTINUE
C
  180       CONTINUE
  200    CONTINUE
  220 CONTINUE
C
C     END OF SCAN IN THE THREE COORDINATE DIRECTIONS FOR THE
C     FORWARD ELIMINATION
C     PROGRESS BACKWARDS TO MULTIPLY BY THE INVERSE OF THE MATRIX U
C
      DO 280 KL = 1, N3
C
C        K RUNS BACKWARDS AND FORWARDS ALTERNATELY
C
         K = KS1*(N3P1-KL) + KS2*KL
         KB = K + KS1 - KS2
C
         DO 260 JL = 1, N2
C
C           J RUNS BACKWARDS AND FORWARDS ALTERNATELY
C
            J = KS1*(N2P1-JL) + KS2*JL
            JB = J + KS1 - KS2
C
            DO 240 IL = 1, N1
C
C              THE INDEX FOR THE FIRST COORDINATE ALWAYS
C              RUNS BACKWARDS SINCE IN THE FIRST LOOP IT
C              ALWAYS RAN FORWARDS
               I = N1P1 - IL
C
C              R IS INITIALLY THE INVERSE OF L * RES AND
C              BECOMES THE INVERSE OF U * INVERSE OF L * RES
C              ,  NAMELY THE CHANGE VECTOR
               IF ((KB.NE.0) .AND. (KB.NE.N3P1)) R(I,J,K) = R(I,J,K) -
     *             WRKSP3(I,J,K)*R(I,J,KB)
               IF ((JB.NE.0) .AND. (JB.NE.N2P1)) R(I,J,K) = R(I,J,K) -
     *             WRKSP2(I,J,K)*R(I,JB,K)
               IF (I.NE.N1) R(I,J,K) = R(I,J,K) - WRKSP1(I,J,K)
     *                                 *R(I+1,J,K)
C
  240       CONTINUE
  260    CONTINUE
  280 CONTINUE
C
C     CORRECT RETURN
C
      IFAIL = 0
      RETURN
C
C     ERROR RETURNS
C
  300 CONTINUE
      IERROR = 1
      GO TO 380
C
  320 CONTINUE
      IERROR = 2
      GO TO 380
C
  340 CONTINUE
      IERROR = 3
      GO TO 380
C
  360 CONTINUE
      IERROR = 4
C
C
  380 CONTINUE
C
C     ERROR CONDITION
C
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END