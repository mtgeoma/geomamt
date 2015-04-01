      SUBROUTINE F02ADF(A,IA,B,IB,N,R,DE,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-736 (DEC 1989).
C
C     EIGENVALUES OF A-LAMBDA*B
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02ADF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), DE(N), R(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XXXX
      INTEGER           ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01AEF, F01AGF, F02AVF
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      CALL F01AEF(N,A,IA,B,IB,DE,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
      CALL F01AGF(N,XXXX,A,IA,R,DE,DE)
      IFAIL = 1
      CALL F02AVF(N,X02AJF(),R,DE,IFAIL)
      IF (IFAIL.EQ.0) GO TO 40
      IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
   40 RETURN
      END
