      DOUBLE PRECISION FUNCTION F02WDZ(N,A,NRA,WORK)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (CONDTN)
C
C     F02WDZ RETURNS THE CONDITION NUMBER OF THE N*N UPPER
C     TRIANGULAR MATRIX A GIVEN BY
C
C     F02WDZ = NORM(A)*NORM(A**(-1)) ,
C
C     WHERE NORM DENOTES THE EUCLIDEAN NORM.
C
C     NRA MUST BE THE ACTUAL ROW DIMENSION OF A AS DECLARED IN
C     THE CALLING PROGRAM AND MUST BE AT LEAST N.
C
C     WORK IS A WORK ARRAY WHOSE LENGTH MUST BE AT LEAST N.
C
C     ONLY THE UPPER TRIANGULAR PART OF A IS REFERENCED.
C
C     IF THE CONDITION NUMBER WOULD OVERFLOW THEN F02WDZ IS
C     RETURNED AS 1.0/SMALL, WHERE SMALL IS THE SMALL POSITIVE
C     REAL NUMBER RETURNED FROM ROUTINE X02AMF.
C
C     .. Scalar Arguments ..
      INTEGER                          N, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(NRA,N), WORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 AINVNM, ANORM, BIG, SMALL
C     .. External Functions ..
      DOUBLE PRECISION                 F04JGR, F04JGS, X02AMF
      EXTERNAL                         F04JGR, F04JGS, X02AMF
C     .. Executable Statements ..
      SMALL = X02AMF()
      BIG = 1.0D0/SMALL
C
      F02WDZ = BIG
C
      AINVNM = F04JGR(N,A,NRA,WORK)
C
      IF (AINVNM.EQ.BIG) RETURN
C
      ANORM = F04JGS(N,A,NRA)
C
      IF (AINVNM.LT.1.0D0) GO TO 20
      IF (ANORM.GE.BIG/AINVNM) RETURN
C
   20 F02WDZ = AINVNM*ANORM
C
      RETURN
      END