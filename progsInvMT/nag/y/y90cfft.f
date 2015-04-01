      SUBROUTINE Y90CFF(TTYPE,UPLO,M,N,A,IA,DIAG,ODIAG,DETMAN,DETEXP,
     *                  DIST,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ============================================
C         *  Y90CFF :  Generate Triangular Matrices  *
C         ============================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      DOUBLE PRECISION  ZERO, HALF
      PARAMETER         (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0),
     *                  ZERO=0.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      COMPLEX*16        DETMAN
      INTEGER           DETEXP, DIST, IA, M, N, TTYPE
      CHARACTER*1       UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), DIAG(*), ODIAG(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, X1, X2, Y, Y1, Y2
      INTEGER           I, J, J1, J2, ND
C     .. External Functions ..
      DOUBLE PRECISION  Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          Y90DMF
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, ANINT, DCMPLX, MIN, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Random Triangular Matrices
C
C-----------------------------------------------------------------------
      DO 40 J = 1, N
         DO 20 I = 1, M
            A(I,J) = CZERO
   20    CONTINUE
   40 CONTINUE
C
      ND = MIN(M,N)
      X1 = DBLE(DIAG(1))
      X2 = DBLE(DIAG(2))
      Y1 = DIMAG(DIAG(1))
      Y2 = DIMAG(DIAG(2))
C
      DO 60 I = 1, ND
         IF (DIST.LE.1) THEN
            X = X1 + (X2-X1)*Y90TBF(DIST,SEED)
            Y = Y1 + (Y2-Y1)*Y90TBF(DIST,SEED)
         ELSE IF (DIST.EQ.2) THEN
            X = HALF*(X1+X2) + (X2-X1)*Y90TBF(DIST,SEED)
            Y = HALF*(Y1+Y2) + (Y2-Y1)*Y90TBF(DIST,SEED)
         ELSE
            X = X1 + X2*Y90TBF(DIST,SEED)
            Y = Y1 + Y2*Y90TBF(DIST,SEED)
         END IF
         IF (TTYPE.LE.1) THEN
            A(I,I) = DCMPLX(X,Y)
         ELSE IF (TTYPE.EQ.2) THEN
            A(I,I) = DCMPLX(ANINT(X),ANINT(Y))
         ELSE IF (TTYPE.EQ.3) THEN
            A(I,I) = DCMPLX(X,ZERO)
         ELSE IF (TTYPE.EQ.4) THEN
            A(I,I) = DCMPLX(ANINT(X),ZERO)
         END IF
   60 CONTINUE
C
      DO 100 I = 1, ND
         IF (Y90WAF(UPLO,'L')) THEN
            J1 = I + 1
            J2 = N
         ELSE
            J1 = 1
            J2 = I - 1
         END IF
         X1 = DBLE(ODIAG(1))
         X2 = DBLE(ODIAG(2))
         Y1 = DIMAG(ODIAG(1))
         Y2 = DIMAG(ODIAG(2))
C
         DO 80 J = J1, J2
            IF (DIST.LE.1) THEN
               X = X1 + (X2-X1)*Y90TBF(DIST,SEED)
               Y = Y1 + (Y2-Y1)*Y90TBF(DIST,SEED)
            ELSE IF (DIST.EQ.2) THEN
               X = HALF*(X1+X2) + (X2-X1)*Y90TBF(DIST,SEED)
               Y = HALF*(Y1+Y2) + (Y2-Y1)*Y90TBF(DIST,SEED)
            ELSE
               X = X1 + X2*Y90TBF(DIST,SEED)
               Y = Y1 + Y2*Y90TBF(DIST,SEED)
            END IF
            IF (TTYPE.LE.1 .OR. TTYPE.EQ.3) THEN
               A(J,I) = DCMPLX(X,Y)
            ELSE
               A(J,I) = DCMPLX(ANINT(X),ANINT(Y))
            END IF
   80    CONTINUE
  100 CONTINUE
C
C     Calculate the determinant
C
      IF (M.EQ.N) THEN
         DETMAN = CONE
         DETEXP = 0
         DO 120 I = 1, M
            DETMAN = DETMAN*A(I,I)
            CALL Y90DMF(DETMAN,DETEXP,4)
  120    CONTINUE
      ELSE
         DETMAN = CZERO
         DETEXP = 0
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90CFF
C
C-----------------------------------------------------------------------
      RETURN
      END
