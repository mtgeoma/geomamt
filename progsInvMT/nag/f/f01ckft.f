      SUBROUTINE F01CKF(A,B,C,N,P,M,Z,IZ,OPT,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     RETURNS WITH THE RESULT OF THE MATRIX MULTIPLICATION OF B AND
C     C
C     IN THE MATRIX A, WITH THE OPTION TO OVERWRITE B OR C.
C     **************************************************************
C     ****
C     MATRIX MULTIPLICATION
C     A=B*C          IF OPT= 1 A,B,C ASSUMED DISTINCT
C     IF OPT=2 B IS OVERWRITTEN
C     IF OPT=3 C IS OVERWRITTEN
C     **************************************************************
C     ****
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01CKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IZ, M, N, OPT, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N,P), B(N,M), C(M,P), Z(IZ)
C     .. Local Scalars ..
      INTEGER           I, IB1, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DGEMV
C     .. Executable Statements ..
      IF (M.GE.1 .AND. N.GE.1 .AND. P.GE.1) GO TO 20
C     *****ERROR CHECK 1 FAILS:- NEGATIVE DIMENSIONS******
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 IF (OPT.EQ.1) GO TO 60
      IF (IZ.GE.M) GO TO 40
      IFAIL = P01ABF(IFAIL,4,SRNAME,0,P01REC)
      RETURN
   40 IF (OPT.EQ.2) GO TO 80
      IF (OPT.EQ.3) GO TO 160
   60 IFAIL = 0
      CALL DGEMM('No transpose','No Transpose',N,P,M,1.0D0,B,N,C,M,
     *           0.0D0,A,N)
      RETURN
C     ******OVERWRITE B******
   80 IF (M.EQ.P) GO TO 100
C     ******ERROR CHECK 2 FAILS:- OVERWRITING INCOMPATIBLE
C     ARRAYS******
      IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      RETURN
  100 IFAIL = 0
      IB1 = N
      IF (M.EQ.1) IB1 = 1
      DO 140 I = 1, N
         CALL DGEMV('T',M,P,1.0D0,C,M,B(I,1),IB1,0.0D0,Z,1)
         DO 120 J = 1, M
            B(I,J) = Z(J)
  120    CONTINUE
  140 CONTINUE
      RETURN
C     ******OVERWRITE C******
  160 IF (N.EQ.M) GO TO 180
C     ******ERROR CHECK 3 FAILS:- OVERWRITING INCOMPATIBLE
C     ARRAYS******
      IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
      RETURN
  180 IFAIL = 0
      DO 220 J = 1, P
         CALL DGEMV('N',M,M,1.0D0,B,N,C(1,J),1,0.0D0,Z,1)
         DO 200 I = 1, M
            C(I,J) = Z(I)
  200    CONTINUE
  220 CONTINUE
      RETURN
      END
