      INTEGER FUNCTION F02BEZ(IP,IQ,D,F,N,AMBDA,ACHEPS)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY W.PHILLIPS     1ST OCTOBER 1975
C     OXFORD UNIVERSITY COMPUTING LABORATORY.
C     THIS ROUTINE REPLACES F02ASZ.
C
C     STURMCNT
C     AUXILIARY ROUTINE CALLED BY F02BEF.
C     CALCULATES THE NUMBER OF EIGENVALUES LESS THAN AMBDA
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION        ACHEPS, AMBDA
      INTEGER                 IP, IQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION        D(N), F(N)
C     .. Local Scalars ..
      DOUBLE PRECISION        ONE, X, ZERO
      INTEGER                 I, ICOUNT
C     .. Intrinsic Functions ..
      INTRINSIC               SQRT
C     .. Data statements ..
      DATA                    ONE/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      ICOUNT = 0
      X = ONE
      DO 60 I = IP, IQ
         IF (X.EQ.ZERO) GO TO 20
         X = D(I) - AMBDA - (F(I)/X)
         GO TO 40
   20    X = D(I) - AMBDA - (SQRT(F(I)))/ACHEPS
   40    IF (X.GE.ZERO) GO TO 60
         ICOUNT = ICOUNT + 1
   60 CONTINUE
      F02BEZ = ICOUNT
      RETURN
      END
