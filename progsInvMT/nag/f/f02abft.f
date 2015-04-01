      SUBROUTINE F02ABF(A,IA,N,R,V,IV,E,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-735 (DEC 1989).
C
C     EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIX MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02ABF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, IV, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), E(N), R(N), V(IV,N)
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
      EXTERNAL          F01AJF, F02AMF
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      CALL F01AJF(N,XXXX,A,IA,R,E,V,IV)
      CALL F02AMF(N,X02AJF(),R,E,V,IV,IFAIL)
      IF (IFAIL.NE.0) IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
