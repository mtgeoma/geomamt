      SUBROUTINE F02BBF(A,IA,N,ALB,UB,M,MM,R,V,IV,D,E,E2,X,G,C,ICOUNT,
     *                  IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-740 (DEC 1989).
C     WRITTEN BY W.PHILLIPS     1ST OCTOBER 1975
C     OXFORD UNIVERSITY COMPUTING LABORATORY.
C     THIS ROUTINE REPLACES F02ACF.
C
C     SELECTED EIGENVALUES AND EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALB, UB
      INTEGER           IA, IFAIL, IV, M, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N), E(N), E2(N), G(N), R(M), V(IV,M),
     *                  X(N,7)
      INTEGER           ICOUNT(M)
      LOGICAL           C(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, XXXX, ZERO
      INTEGER           I, ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01AGF, F01AHF, F02BEF
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      EPS = ZERO
      CALL F01AGF(N,XXXX,A,IA,D,E,E2)
      DO 20 I = 1, N
         G(I) = E(I)
   20 CONTINUE
      CALL F02BEF(N,D,ALB,UB,X02AJF(),EPS,E,E2,M,MM,R,V,IV,ICOUNT,X,C,
     *            IFAIL)
      IF (IFAIL.EQ.0) GO TO 40
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   40 IF (MM.EQ.0) RETURN
      CALL F01AHF(N,1,MM,A,IA,G,V,IV)
      RETURN
      END
