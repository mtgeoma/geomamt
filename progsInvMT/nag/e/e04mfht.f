      SUBROUTINE E04MFH(N,NCLIN,LDA,ISTATE,BIGBND,NUMINF,SUMINF,BL,BU,A,
     *                  FEATOL,CVEC,X,WTINF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFH  finds the number and weighted sum of infeasibilities for
C     the bounds and linear constraints.   An appropriate gradient
C     is returned in cvec.
C
C     Positive values of  istate(j)  will not be altered.  These mean
C     the following...
C
C               1             2           3
C           a'x = bl      a'x = bu     bl = bu
C
C     Other values of  istate(j)  will be reset as follows...
C           a'x lt bl     a'x gt bu     a'x free
C              - 2           - 1           0
C
C     Original version written 31-October-1984.
C     This version of E04MFH dated  1-January-1987.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, SUMINF
      INTEGER           LDA, N, NCLIN, NUMINF
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(*), BU(*), CVEC(N), FEATOL(*),
     *                  WTINF(*), X(N)
      INTEGER           ISTATE(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CTX, FEASJ, S, WEIGHT
      INTEGER           J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DAXPY, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      BIGUPP = BIGBND
      BIGLOW = -BIGBND
C
      NUMINF = 0
      SUMINF = ZERO
      CALL F06FBF(N,(ZERO),CVEC,1)
C
      DO 40 J = 1, N + NCLIN
         IF (ISTATE(J).LE.0) THEN
            FEASJ = FEATOL(J)
            IF (J.LE.N) THEN
               CTX = X(J)
            ELSE
               K = J - N
               CTX = DDOT(N,A(K,1),LDA,X,1)
            END IF
            ISTATE(J) = 0
C
C           See if the lower bound is violated.
C
            IF (BL(J).GT.BIGLOW) THEN
               S = BL(J) - CTX
               IF (S.GT.FEASJ) THEN
                  ISTATE(J) = -2
                  WEIGHT = -WTINF(J)
                  GO TO 20
               END IF
            END IF
C
C           See if the upper bound is violated.
C
            IF (BU(J).GE.BIGUPP) GO TO 40
            S = CTX - BU(J)
            IF (S.LE.FEASJ) GO TO 40
            ISTATE(J) = -1
            WEIGHT = WTINF(J)
C
C           Add the infeasibility.
C
   20       NUMINF = NUMINF + 1
            SUMINF = SUMINF + ABS(WEIGHT)*S
            IF (J.LE.N) THEN
               CVEC(J) = WEIGHT
            ELSE
               CALL DAXPY(N,WEIGHT,A(K,1),LDA,CVEC,1)
            END IF
         END IF
   40 CONTINUE
      RETURN
C
C     End of  E04MFH.  (CMSINF)
C
      END
