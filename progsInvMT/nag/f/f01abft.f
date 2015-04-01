      SUBROUTINE F01ABF(A,IA,N,B,IB,Z,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 16 REVISED. IER-1001 (JUN 1993).
C
C     ACCURATE INVERSE OF A REAL SYMMETRIC POSITIVE DEFINITE
C     MATRIX USING CHOLESKYS METHOD, SIMPLIFIED PARAMETERS.
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01ABF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, XXXX
      INTEGER           I, ISAVE, J, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01ACZ
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      EPS = X02AJF()
      CALL F01ACZ(N,EPS,A,IA,B,IB,Z,L,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 DO 60 I = 1, N
         DO 40 J = I, N
            B(J,I) = A(J+1,I)
   40    CONTINUE
   60 CONTINUE
      RETURN
      END
