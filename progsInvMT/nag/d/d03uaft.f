      SUBROUTINE D03UAF(N1,N2,N1M,A,B,C,D,E,APARAM,IT,R,WRKSP1,WRKSP2,
     *                  IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 10A REVISED. IER-386 (OCT 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C     D03UAF PERFORMS 1 ITERATION OF THE STRONGLY IMPLICIT PROCEDURE
C     AT EACH CALL TO CALCULATE THE SUCCESSIVE APPROXIMATE
C     CORRECTIONS TO THE SOLUTION OF A SYSTEM OF SIMULTANEOUS
C     ALGEBRAIC EQUATIONS FOR WHICH THE ITERATIVE UPDATE MATRIX IS
C     OF THE FIVE POINT MOLECULE FORM ON A TOPOLOGICALLY TWO-
C     DIMENSIONAL RECTANGULAR MESH.
C
C     STRONGLY IMPLICIT PROCEDURE ROUTINE FOR 2 DIMENSIONAL 05
C     POINT MOLECULES.
C
C     INPUTS
C
C     N1      NUMBER OF NODES IN THE FIRST COORDINATE DIRECTION.
C     N2      NUMBER OF NODES IN THE SECOND COORDINATE DIRECTION.
C     N1M     FIRST DIMENSION OF ALL THE TWO-DIMENSIONAL ARRAYS.
C     A       ARRAY OF DIMENSION (N1M,N2) STORING THE COEFFICIENT
C             OF THE ITERATIVE UPDATE EQUATIONS AS SHOWN BELOW -
C     B       ... DIMENSION (N1M,N2) ...
C     C       ... DIMENSION (N1M,N2) ...
C     D       ... DIMENSION (N1M,N2) ...
C     E       ... DIMENSION (N1M,N2) ...
C
C     A(I,J)*S(I,J-1)+B(I,J)*S(I-1,J)+C(I,J)*S(I,J)+D(I,J)*
C      S(I+1,J)+E(I,J)*S(I,J+1)=R(I,J)
C
C             WITH I=1,2,...,N1 AND J=1,2,...,N2 AND WHERE S(I,J)
C             IS THE ARRAY WHOSE APPROXIMATE VALUES ARE SOUGHT
C             (AND WHICH OVERWRITE THE RESIDUALS R(I,J) STORED ON
C             INPUT IN THE ARRAY R). ANY VALUES OF S OUTSIDE THE
C             (1-N1,1-N2) ARRAY WHICH DEFINES THE PROBLEM REGION
C             ARE TAKEN AS ZERO.
C     R       IS INPUT AS THE ARRAY OF RESIDUALS, CORRESPONDING
C             TO R ABOVE, CHANGED ON OUTPUT TO THE VALUES OF THE
C             UPDATE ARRAY CORRESPONDING TO S ABOVE, DIMENSIONED
C             (N1M,N2).
C     APARAM  IS AN ITERATION ACCELERATION PARAMETER FACTOR
C             TYPICALLY SET TO 1.0. IF CONVERGENCE IS SLOW, IT
C             CAN BE DECREASED. IF DIVERGENCE IS OBTAINED, IT
C             SHOULD BE INCREASED. IN EITHER CASE IT MUST NOT GO
C             OUTSIDE THE BOUNDS PRESCRIBED (SEE IFAIL PARAMETER
C             FOR D03UAF).
C     IT      IS THE ITERATION COUNTER SET AND INCREMENTED BY
C             THE CALLING ROUTINE, IT IS USED TO DETERMINE THE
C             APPROPRIATE ITERATION PARAMETERS.
C     WRKSP1  IS A WORKSPACE ARRAY OF DIMENSION (N1M,N2).
C     WRKSP2  ... DITTO ...
C     IFAIL   IS AN ERROR PARAMETER INDICATOR SET BY THE USER TO
C             INDICATE THE TYPE OF FAILURE IF AN ERROR IS
C             ENCOUNTERED.
C
C     PROCESS
C
C     SET ERROR PARAMETER
C     CHECK INPUT INTEGERS
C     SET FREQUENTLY REQUIRED INTEGER VARIABLES
C     COMMENCEMENT OF CALCULATIONAL PROCEDURE
C     SET ODD/EVEN COUNTERS  KS=1 FOR ODD ITERATIONS
C                            KS=2 FOR EVEN ITERATIONS
C     DETERMINE THE NUMBER OF THE ACCELERATION PARAMETER TO BE USED
C     (THE SAME PARAMETER IS USED FOR 2 ITERATIONS, THERE ARE 9
C     PARAMETERS IN ALL).
C     FIRST CALCULATE THE TERM, ALM, IN THE LARGEST PARAMETER
C     CHECK THE VALUES OF ALM
C     DETERMINE THE ITERATION PARAMETER
C     START THE APPROXIMATE LU FACTORIZATION DETERMINING THE VALUES
C     OF THE ELEMENTS SB,SC IN THE LOWER TRIANGULAR MATRIX L AND
C     WRKSP1 AND WRKSP2 IN THE UPPER TRIANGULAR MATRIX U. PROGRESS
C     FORWARDS FIRST AND INVERT THE LOWER TRIANGULAR MATRIX L AS
C     ONE PROCEEDS SO THAT THE ELEMENTS SB,SC DO NOT HAVE TO BE
C     STORED IN ARRAYS. THE ELEMENTS WRKSP1 AND WRKSP2 HAVE TO BE
C     STORED FOR THE SUBSEQUENT INVERSION OF THE MATRIX U BY BACK
C     SUBSTITUTION.
C     DO JL = 1,N2
C        (JL IS ONLY A COUNTER IN THE SECOND COORDINATE DIRECTION
C         J IS THE INDEX OF THE SECOND COORDINATE AND ON
C         ALTERNATE ITERATIONS SCANS FIRST INCREASING AND THEN
C         DECREASING).
C         DO I = 1,N1
C            (I IS THE INDEX IN THE FIRST COORDINATE DIRECTION)
C             STORE TH VALUE OF THE CENTRAL COEFFICIENT
C             DETERMINE THE TYPE OF THE DIFFERENCE EQUATION
C               FOR A FIVE POINT MOLECULE EQUATION -
C                  SET VALUES OF ADJACENT ARRAY ELEMENTS DEPENDENT
C                  ON NODAL POSITION
C                  CALCULATE THE OFF-DIAGONAL ELEMENTS OF L,
C                  NAMELY SB AND SC
C                  CALCULATE THE OFTEN USED EXPRESSIONS
C                  CALCULATE THE VALUE OF THE DIAGONAL ELEMENT OF L
C                  CALCULATE AND STORE THE ELEMENTS OF U, NAMELY
C                  WRKSP1 AND WRKSP2
C                  CALCULATE THE ELEMENTS OF L**(-1) * RESIDUAL
C                  ARRAY
C               FOR THE EXPLICIT EQUATION DO THE SAME
C     END OF SCAN FOR FORWARD ELIMINATION
C     BACK SUBSTITUTION TO MULTIPLY BY U**(-1)
C     DO JL = 1,N2
C        J RUNS BACKWARDS AND FORWARDS ALTERNATELY
C        DO I = N1,N1-1,....,2,1
C           R IS INITIALLY L**(-1) * RESIDUAL, IT BECOMES
C           (U**(-1) * (L**(-1) *RESIDUAL) = CHANGE
C     END OF SCAN FOR BACK SUBSTITUTION
C     RETURN
C
C     OUTPUTS
C
C     R       ARRAY STORING THE APPROXIMATE SOLUTION TO THE SYSTEM
C             OF EQUATIONS PROVIDED, AFTER ONE ITERATION. NOTE
C             THAT IT HAS OVERWRITTEN THE INPUT RESIDUAL.
C     IFAIL   ERROR INDICATOR
C             =0 CORRECT RETURN
C             =1 EITHER N1.LE.1 OR N2.LE.1
C             =2 N1M IS LESS THAN N1
C             =3 APARAM.LE.0.0
C             =4 APARAM.GT.((N1-1)**2+N2-1)**2)/2.
C
C     ROUTINES USED
C
C     P01AAF  ERROR HANDLING ROUTINE
C
C     **************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03UAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  APARAM
      INTEGER           IFAIL, IT, N1, N1M, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N1M,N2), B(N1M,N2), C(N1M,N2), D(N1M,N2),
     *                  E(N1M,N2), R(N1M,N2), WRKSP1(N1M,N2),
     *                  WRKSP2(N1M,N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALM, ALPHA, CS, RIBJL, RIJB, SB, SBSEIJ, SC,
     *                  SCSFIJ, SD, SEIBJL, SEIJB, SFIBJL, SFIJB, XKS1,
     *                  XKS2
      INTEGER           I, IB, IERROR, IL, IS, J, JB, JL, KS, KS1, KS2,
     *                  N1M1, N1P1, N2M1, N2P1
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
      IERROR = 0
C
C     CHECK INPUT INTEGERS
C
      IF (N1.LE.1) GO TO 220
      IF (N2.LE.1) GO TO 220
C
      IF (N1M.LT.N1) GO TO 240
C
C     SET FREQUENTLY REQUIRED INTEGER VARIABLES
C
      N1M1 = N1 - 1
      N2M1 = N2 - 1
      N1P1 = N1 + 1
      N2P1 = N2 + 1
C
C     COMMENCEMENT OF CALCULATIONAL PROCEDURE
C     ---------------------------------------
C
C
C     SET ODD/EVEN COUNTER    KS=1  FOR ODD ITERATIONS
C     KS=2  FOR EVEN ITERATIONS
      KS = MOD(IT-1,2) + 1
      IF (KS.EQ.0) KS = 2
C
C
C     SET FREQUENTLY USED ODD/EVEN PARAMETER DERIVATIVES
C
      KS1 = 2 - KS
      XKS1 = DBLE(KS1)
      KS2 = KS - 1
      XKS2 = DBLE(KS2)
C
C     DETERMINE THE NUMBER OF THE ACCELERATION PARAMETER TO BE USED,
C     THE SAME PARAMETER IS USED FOR 2 ITERATIONS, THERE ARE 9
C     PARAMETERS IN ALL
C
      IS = MOD(IT-1,18)
      IF (IS.LT.0) IS = IS + 18
      IS = IS/2 + 1
C
C     CALCULATION OF THE ITERATION ACCELERATION PARAMETER
C
C     (1) CALCULATE THE TERM IN THE LARGEST PARAMETER
C
      ALM = 2.D0*APARAM/(DBLE(N1M1*N1M1+N2M1*N2M1))
C
C     (2) CHECK THE VALUES OF ALM
C
      IF (APARAM.LE.0.0D0) GO TO 260
C
      IF (ALM.GT.1.0D0) GO TO 280
C
C     (3) THEN DETERMINE THE ITERATION PARAMETER
C
      ALPHA = 1.D0 - ALM**ALP(IS)
C
C     START OF APPROXIMATE LU FACTORIZATION DETERMINING THE VALUES
C     OF THE ELEMENTS SB,SC IN THE LOWER TRIANGULAR MATRIX L AND
C     WRKSP1 AND WRKSP2 IN THE UPPER TRIANGULAR MATRIX U. PROGRESS
C     FORWARDS FIRST AND INVERT THE LOWER TRIANGULAR MATRIX L AS ONE
C     PROCEEDS SO THAT THE ELEMENTS SB,SC DO NOT HAVE TO BE STORED
C     IN ARRAYS. THE ELEMENTS WRKSP1 AND WRKSP2 HAVE TO BE STORED
C     FOR THE SUBSEQUENT INVERSION OF THE MATRIX U BY BACK
C     SUBSTITUTION.
C
      DO 160 JL = 1, N2
C
C        JL IS ONLY A COUNTER IN THE SECOND COORDINATE DIRECTION
C        J IS THE INDEX OF THE SECOND COORDINATE AND ON ALTERNATE
C        ITERATIONS SCANS FIRST INCREASING AND THEN DECREASING
C
         J = KS1*JL + KS2*(N2P1-JL)
         JB = J - KS1 + KS2
C
         DO 140 I = 1, N1
C
C           I IS THE INDEX IN THE FIRST COORDINATE DIRECTION
C
            IB = I - 1
C
C           STORE THE VALUE OF THE CENTRAL COEFFICIENT
C
            CS = C(I,J)
C
C           DETERMINE THE TYPE OF THE DIFFERENCE EQUATION
C
            IF (CS.EQ.0.0D0) GO TO 100
C
C           FIVE POINT MOLECULE EQUATION
C
C           SET VALUES OF ADJACENT ARRAY ELEMENTS DEPENDENT ON NODAL
C           POSITION
C
            IF ((JB.EQ.0) .OR. (JB.EQ.N2P1)) GO TO 20
C
C           (1A) IF NOT ON  BOUNDARY  SO THAT JB.NE.0 OR N2P1
C
            SEIJB = WRKSP1(I,JB)
            SFIJB = WRKSP2(I,JB)
            RIJB = R(I,JB)
            GO TO 40
C
   20       CONTINUE
C
C           (1B) THE VALUES EXTERIOR TO THE ARRAYS ARE SET TO ZERO
C
            SEIJB = 0.0D0
            SFIJB = 0.0D0
            RIJB = 0.0D0
   40       CONTINUE
C
            IF (I.EQ.1) GO TO 60
C
C           (2A) IF NOT ON I=1  I-1.NE.0  SO THAT
C
            SEIBJL = WRKSP1(IB,J)
            SFIBJL = WRKSP2(IB,J)
            RIBJL = R(IB,J)
            GO TO 80
C
   60       CONTINUE
C
C           (2B) ON  I=1  THE VALUES EXTERIOR TO THE ARRAY ARE SET TO 0
C
            SEIBJL = 0.0D0
            SFIBJL = 0.0D0
            RIBJL = 0.0D0
   80       CONTINUE
C
C           CALCULATE THE ELEMENTS OF THE LOWER TRIANGULAR MATRIX
C
            SB = (XKS1*A(I,J)+XKS2*E(I,J))/(1.D0+ALPHA*SEIJB)
            SC = B(I,J)/(1.D0+ALPHA*SFIBJL)
C
C           CALCULATE OFTEN USED EXPRESIONS
C
            SBSEIJ = SB*SEIJB
            SCSFIJ = SC*SFIBJL
C
C           CALCULATE THE VALUE OF THE DIAGONAL ELEMENT OF L
C
            SD = 1.D0/(-SB*SFIJB-SC*SEIBJL+CS+ALPHA*(SBSEIJ+SCSFIJ))
C
C           CALCULATE AND STORE THE ELEMENTS OF THE UPPER TRIANGULAR
C           MATRIX
C
            WRKSP1(I,J) = (D(I,J)-ALPHA*SBSEIJ)*SD
            WRKSP2(I,J) = (XKS1*E(I,J)+XKS2*A(I,J)-ALPHA*SCSFIJ)*SD
C
C           CALCULATION OF THE ELEMENTS OF THE INVERSE OF L * RES VECTOR
C
            R(I,J) = (R(I,J)-SB*RIJB-SC*RIBJL)*SD
C
            GO TO 120
C
  100       CONTINUE
C
C           CALCULATION FOR THE EXPLICIT EQUATION
C
            WRKSP1(I,J) = 0.0D0
            WRKSP2(I,J) = 0.0D0
C
  120       CONTINUE
C
  140    CONTINUE
  160 CONTINUE
C
C     END OF SCAN IN THE TWO COORDINATE DIRECTIONS FOR THE FORWARD
C     ELIMINATION
C
C     PROGRESS BACKWARDS TO MULTIPLY BY THE INVERSE OF THE MATRIX U
C
      DO 200 JL = 1, N2
C
C        J RUNS BACKWARDS AND FORWARDS ALTERNATELY
C
         J = KS1*(N2P1-JL) + KS2*JL
         JB = J + KS1 - KS2
C
         DO 180 IL = 1, N1
C
C           THE INDEX FOR THE FIRST COORDINATE ALWAYS RUNS BACKWARDS
C           SINCE IN THE FIRST LOOP IT ALWAYS RAN FORWARDS
C
            I = N1P1 - IL
C
C           R IS INITIALLY THE INVERSE OF L * RESIDUAL
C           IT BECOMES THE INVERSE OF U * INVERSE OF L * RESIDUAL
C           NAMELY THE CHANGE VECTOR
C
            IF ((JB.NE.0) .AND. (JB.NE.N2P1)) R(I,J) = R(I,J) -
     *          WRKSP2(I,J)*R(I,JB)
C
            IF (I.NE.N1) R(I,J) = R(I,J) - WRKSP1(I,J)*R(I+1,J)
C
  180    CONTINUE
  200 CONTINUE
C
C     CORRECT RETURN
C
      IFAIL = 0
      RETURN
C
C     ERROR RETURNS
C
  220 CONTINUE
      IERROR = 1
      GO TO 300
C
  240 CONTINUE
      IERROR = 2
      GO TO 300
C
  260 CONTINUE
      IERROR = 3
      GO TO 300
C
  280 CONTINUE
      IERROR = 4
C
C     ERROR CONDITION
C
  300 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
C
      RETURN
C
      END
