      SUBROUTINE F01ACZ(N,EPS,A,IA,B,IB,Z,L,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993
C
C     This routine originally called F01ACF.
C     ACCINVERSE
C     THE UPPER TRIANGLE OF A POSITIVE DEFINITE SYMMETRIC MATRIX,
C     A, IS STORED IN THE UPPER TRIANGLE OF AN (N+1)*N ARRAY
C     A(I,J), I=1,N+1, J=1,N. X, THE INVERSE OF A, IS FORMED IN
C     THE REMAINDER OF THE ARRAY A(I,J) BY THE SUBROUTINE F01ADF
C     CHOLINVERSION 1. THE INVERSE IS IMPROVED BY CALCULATING X=X+Z
C     UNTIL THE CORRECTION, Z, IS SUCH THAT MAXIMUM ABS(Z(I,J)) IS
C     LESS THAN 2 EPS TIMES MAXIMUM ABS(X(I,J)), WHERE Z=XB AND
C     B=I-AX. B IS AN N*N ARRAY AND Z IS A 1*N ARRAY, X BEING
C     OVERWRITTEN A ROW AT A TIME. EXITS WITH IFAIL = 1 IF A IS
C     NOT POSITIVE DEFINITE AND WITH IFAIL = 2 IF THE MAXIMUM
C     CORRECTION AT ANY STAGE IS NOT LESS THAN HALF THAT AT THE
C     PREVIOUS STAGE. L IS THE NUMBER OF CORRECTIONS APPLIED.
C     ADDITIONAL PRECISION INNERPRODUCTS ARE ABSOLUTELY NECESSARY.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01ACZ')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IA, IB, IFAIL, L, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, C1, C2, D, D1, D2, E, XMAX, ZMAX
      INTEGER           I, IFAIL1, ISAVE, J, J1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01ADF, X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL1 = 0
      E = 1.0D0
      L = 0
      IFAIL = 1
      CALL F01ADF(N,A,IA,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
   20 DO 100 I = 1, N
         DO 80 J = 1, N
            J1 = J + 1
            C1 = 0.0D0
            IF (J.LT.I) GO TO 40
            IF (I.EQ.J) C1 = -1.0D0
            CALL X03AAF(A(1,I),(N-I+1)*IA,A(J1,1)
     *                  ,N*IA-J1+1,I,1,IA,C1,0.0D0,D1,D2,.TRUE.,IFAIL1)
            IF (I.EQ.N) GO TO 60
            C1 = D1
            C2 = D2
            CALL X03AAF(A(I,I+1),(N-I)*IA-I+1,A(J1,I+1),(N-I)
     *                  *IA-J1+1,J-I,IA,IA,C1,C2,D1,D2,.TRUE.,IFAIL1)
            IF (J1.GT.N) GO TO 60
            C1 = D1
            C2 = D2
            CALL X03AAF(A(I,J1),(N-J)*IA-I+1,A(J1+1,J),(N-J+1)
     *                  *IA-J1,N-J,IA,1,C1,C2,D1,D2,.TRUE.,IFAIL1)
            GO TO 60
   40       CALL X03AAF(A(1,I),(N-I+1)*IA,A(J1,1)
     *                  ,N*IA-J1+1,J,1,IA,C1,0.0D0,D1,D2,.TRUE.,IFAIL1)
            C1 = D1
            C2 = D2
            CALL X03AAF(A(J1,I),(N-I+1)*IA-J1+1,A(J1+1,J),(N-J+1)
     *                  *IA-J1,I-J1+1,1,1,C1,C2,D1,D2,.TRUE.,IFAIL1)
            IF (I.EQ.N) GO TO 60
            C1 = D1
            C2 = D2
            CALL X03AAF(A(I,I+1),(N-I)*IA-I+1,A(I+2,J),(N-J+1)
     *                  *IA-I-1,N-I,IA,1,C1,C2,D1,D2,.TRUE.,IFAIL1)
   60       B(I,J) = -D1
   80    CONTINUE
  100 CONTINUE
      XMAX = 0.0D0
      ZMAX = 0.0D0
      DO 180 I = 1, N
         DO 140 J = 1, I
            CALL X03AAF(A(I+1,1),N*IA-I,B(1,J),(N-J+1)
     *                  *IB,I,IA,1,0.0D0,0.0D0,D1,D2,.TRUE.,IFAIL1)
            IF (I.EQ.N) GO TO 120
            C1 = D1
            C2 = D2
            CALL X03AAF(A(I+2,I),(N-I+1)*IA-I-1,B(I+1,J),(N-J+1)
     *                  *IB-I,N-I,1,1,C1,C2,D1,D2,.TRUE.,IFAIL1)
  120       Z(J) = D1
  140    CONTINUE
         DO 160 J = 1, I
            C = ABS(A(I+1,J))
            D = ABS(Z(J))
            IF (C.GT.XMAX) XMAX = C
            IF (D.GT.ZMAX) ZMAX = D
            A(I+1,J) = A(I+1,J) + Z(J)
  160    CONTINUE
  180 CONTINUE
      L = L + 1
      D = ZMAX/XMAX
      IF (D.GT.E/2.0D0) GO TO 200
      E = D
      IF (D.GT.2.0D0*EPS) GO TO 20
      IFAIL = 0
      RETURN
  200 IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
      END
