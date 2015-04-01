      SUBROUTINE F05AAF(A,IA,M,N1,N2,S,CC,ICOL,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 10A REVISED. IER-390 (OCT 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     PERFORMS SCHMIDT ORTHONORMALISATION ON VECTORS N1
C     TO N2 OF ARRAY A.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F05AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CC
      INTEGER           IA, ICOL, IFAIL, M, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N2), S(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D1, D2, SCALE
      INTEGER           I, IFAIL1, ISAVE, J, JM1, K, N, N11
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          X03AAF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 0
      IFAIL1 = 0
      IF (N1.GT.N2) GO TO 240
      CC = 0.0D0
      ICOL = N1
C     NORMALISE VECTORS
      DO 60 J = N1, N2
         CALL X03AAF(A(1,J),M,A(1,J),M,M,1,1,0.0D0,0.0D0,SCALE,D2,
     *               .FALSE.,IFAIL1)
         IF (SCALE.GT.0.0D0) GO TO 20
         CC = 1.0D0
         ICOL = J
         GO TO 60
   20    SCALE = 1.0D0/SQRT(SCALE)
         DO 40 I = 1, M
            A(I,J) = SCALE*A(I,J)
   40    CONTINUE
   60 CONTINUE
C     IF ONLY ONE VECTOR SKIP ORTHOGONALISATION
      IF (N1.EQ.N2) RETURN
C     FOR EACH VECTOR ORTHOGONALISE WITH RESPECT TO PREVIOUS
C     VECTORS AND THEN RENORMALISE
      N11 = N1 + 1
      DO 220 J = N11, N2
         JM1 = J - 1
   80    DO 100 K = N1, JM1
            CALL X03AAF(A(1,J),M,A(1,K),M,M,1,1,0.0D0,0.0D0,D1,D2,
     *                  .TRUE.,IFAIL1)
            S(K) = D1
  100    CONTINUE
         N = J - N1
C        ORTHOGONALISATION
         DO 120 I = 1, M
            C = -A(I,J)
            CALL X03AAF(A(I,N1),(N-1)*IA+1,S(N1)
     *                  ,N,N,IA,1,C,0.0D0,D1,D2,.TRUE.,IFAIL1)
            A(I,J) = -D1
  120    CONTINUE
C        RE-NORMALISATION
         CALL X03AAF(A(1,J),M,A(1,J),M,M,1,1,0.0D0,0.0D0,SCALE,D2,
     *               .FALSE.,IFAIL1)
         IF (SCALE.GT.0.0D0) GO TO 140
         CC = 1.0D0
         ICOL = J
         GO TO 180
  140    SCALE = 1.0D0/SQRT(SCALE)
         DO 160 I = 1, M
            A(I,J) = SCALE*A(I,J)
  160    CONTINUE
  180    CALL X03AAF(S(N1),N,S(N1),N,N,1,1,0.0D0,0.0D0,C,D2,.FALSE.,
     *               IFAIL1)
         IF (C.LE.CC) GO TO 200
         CC = C
         ICOL = J
C        REPEAT ORTHONORMALISATION IF C IS TOO LARGE
  200    IF (C.GE.0.5D0) GO TO 80
  220 CONTINUE
      RETURN
  240 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
