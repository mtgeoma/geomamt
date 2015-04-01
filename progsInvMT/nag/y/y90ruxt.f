      SUBROUTINE Y90RUX(N,P,A,IROW,ICOL,IDIMA,ISTR,ISTC,AP)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
*-----------------------------------------------------------------------
*
*         ====================================
*         *  Y90RUX :  Auxiliary for Y90RUF  *
*         ====================================
*
*     Purpose
*     =======
*     Maps row P of a real symmetric matrix stored in linked-list format
*     into dense format.
*
*-----------------------------------------------------------------------
*
*     Initialize AP.
*
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IDIMA, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDIMA), AP(N)
      INTEGER           ICOL(IDIMA), IROW(IDIMA), ISTC(N), ISTR(N)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Executable Statements ..
      CALL F06FBF(N,ZERO,AP,1)
*
*     Copy lower triangular non-zero elements of row P into AP.
*
      I = ISTR(P)
   20 CONTINUE
      J = I
   40 CONTINUE
      J = ICOL(J)
      IF (J.GT.0) GO TO 40
      J = -J
      AP(J) = A(I)
      I = IROW(I)
      IF (I.GT.0) GO TO 20
*
*     Copy upper triangular non-zero elements of row P into AP.
*
      I = ISTC(P)
   60 CONTINUE
      J = I
   80 CONTINUE
      J = IROW(J)
      IF (J.GT.0) GO TO 80
      J = -J
      AP(J) = A(I)
      I = ICOL(I)
      IF (I.GT.0) GO TO 60
*
*     End of subroutine Y90RUX
*
      RETURN
      END
