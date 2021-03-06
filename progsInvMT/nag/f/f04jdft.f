      SUBROUTINE F04JDF(M,N,A,NRA,B,TOL,SIGMA,IRANK,WORK,LWORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 17 REVISED. IER-1676 (JUN 1995).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (SVDLS2)
C
C     F04JDF RETURNS THE N ELEMENT VECTOR X, OF MINIMAL
C     LENGTH, THAT MINIMIZES THE EUCLIDEAN LENGTH OF THE M
C     ELEMENT VECTOR R GIVEN BY
C
C     R = B-A*X ,
C
C     WHERE A IS AN M*N (M.LE.N) MATRIX AND B IS AN M ELEMENT
C     VECTOR. X IS OVERWRITTEN ON B.
C
C     THE SOLUTION IS OBTAINED VIA A SINGULAR VALUE
C     DECOMPOSITION (SVD) OF THE MATRIX A GIVEN BY
C
C     A = Q*(D 0)*(P**T) ,
C
C     WHERE Q AND P ARE ORTHOGONAL AND D IS A DIAGONAL MATRIX WITH
C     NON-NEGATIVE DIAGONAL ELEMENTS, THESE BEING THE SINGULAR
C     VALUES OF A.
C
C     INPUT PARAMETERS.
C
C     M     - NUMBER OF ROWS OF A. M MUST BE AT LEAST UNITY.
C
C     N     - NUMBER OF COLUMNS OF A. N MUST BE AT LEAST M.
C
C     A     - AN M*N REAL MATRIX.
C
C     NRA   - ROW DIMENSION OF A AS DECLARED IN THE CALLING PROGRAM.
C             NRA MUST BE AT LEAST M.
C
C     B     - AN N ELEMENT REAL VECTOR.
C             THE FIRST M ELEMENTS OF B MUST CONTAIN THE
C             VECTOR B OF THE LEAST SQUARES PROBLEM.
C
C
C     TOL   - A RELATIVE TOLERANCE USED TO DETERMINE THE RANK OF A.
C             TOL SHOULD BE CHOSEN AS APPROXIMATELY THE
C             LARGEST RELATIVE ERROR IN THE ELEMENTS OF A.
C             FOR EXAMPLE IF THE ELEMENTS OF A ARE CORRECT
C             TO ABOUT 4 SIGNIFICANT FIGURES THEN TOL
C             SHOULD BE CHOSEN AS ABOUT 5.0*10.0**(-4).
C
C     IFAIL - THE USUAL FAILURE PARAMETER. IF IN DOUBT SET
C             IFAIL TO ZERO BEFORE CALLING THIS ROUTINE.
C
C     OUTPUT PARAMETERS.
C
C     A     - A WILL CONTAIN THE FIRST M ROWS OF THE
C             ORTHOGONAL MATRIX P**T OF THE SVD.
C
C     B     - B WILL CONTAIN THE MINIMAL LEAST SQUARES
C             SOLUTION VECTOR X.
C
C     SIGMA - IF M IS GREATER THAN IRANK THEN SIGMA WILL CONTAIN THE
C             STANDARD ERROR GIVEN BY
C             SIGMA=L(R)/SQRT(M-IRANK), WHERE L(R) DENOTES
C             THE EUCLIDEAN LENGTH OF THE RESIDUAL VECTOR
C             R. IF M=IRANK THEN SIGMA IS RETURNED AS ZERO.
C
C     IRANK - THE RANK OF THE MATRIX A.
C
C     IFAIL - ON NORMAL RETURN IFAIL WILL BE ZERO.
C             IN THE UNLIKELY EVENT THAT THE QR-ALGORITHM
C             FAILS TO FIND THE SINGULAR VALUES IN 50*M
C             ITERATIONS THEN IFAIL IS SET TO 2.
C             IF AN INPUT PARAMETER IS INCORRECTLY SUPPLIED
C             THEN IFAIL IS SET TO UNITY.
C
C     WORKSPACE PARAMETERS.
C
C     WORK  - AN (M*M+4*M) ELEMENT VECTOR.
C             ON RETURN THE FIRST M ELEMENTS OF WORK WILL
C             CONTAIN THE SINGULAR VALUES OF A ARRANGED IN
C             DESCENDING ORDER. WORK(M+1) WILL CONTAIN THE
C             TOTAL NUMBER OF ITERATIONS TAKEN BY THE
C             QR-ALGORITHM.
C
C     LWORK - THE LENGTH OF THE VECTOR WORK.
C             LWORK MUST BE AT LEAST M*M+4*M.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04JDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SIGMA, TOL
      INTEGER           IFAIL, IRANK, LWORK, M, N, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,N), B(N), WORK(LWORK)
C     .. Local Scalars ..
      INTEGER           IERR, K, MP1, MP2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           F02WDY, P01ABF
      EXTERNAL          F02WDY, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02WBX, F04JAZ
C     .. Executable Statements ..
      IERR = IFAIL
      IF (IERR.EQ.0) IFAIL = 1
C
      IF (NRA.LT.M .OR. N.LT.M .OR. M.LT.1 .OR. LWORK.LT.M*(M+4))
     *    GO TO 20
C
      K = LWORK - M
      MP1 = M + 1
      MP2 = MP1 + 1
C
      CALL F02WBX(M,N,A,NRA,.TRUE.,B,WORK,WORK(MP1),K,IFAIL)
C
      IF (IFAIL.NE.0) GO TO 20
C
      IRANK = F02WDY(M,WORK,TOL)
C
      CALL F04JAZ(M,N,IRANK,WORK,M,B,A,NRA,B,SIGMA,WORK(MP2))
C
      RETURN
C
   20 IFAIL = P01ABF(IERR,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
