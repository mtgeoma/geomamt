      SUBROUTINE E04MFQ(N,NCLIN,ISTATE,BIGBND,NVIOL,JMAX,ERRMAX,AX,BL,
     *                  BU,FEATOL,X)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1561 (JUN 1995).
C
C     ******************************************************************
C     E04MFQ  checks the residuals of the constraints that are believed
C     to be feasible.  The number of constraints violated by more than
C     featol is computed, along with the maximum constraint violation.
C
C     Original version written by PEG,   April    1984.
C     This version of  E04MFQ  dated  30-Jun-1988.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, ERRMAX
      INTEGER           JMAX, N, NCLIN, NVIOL
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATOL(N+NCLIN), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CON, FEASJ, RES
      INTEGER           IS, J
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Compute the number of constraints (NVIOL) violated by more than
C     FEATOL and  the maximum constraint violation (ERRMAX).
C     (The residual of a constraint in the working set is treated as if
C     it were an equality constraint fixed at that bound.)
C     ==================================================================
      NVIOL = 0
      JMAX = 0
      ERRMAX = ZERO
C
      DO 40 J = 1, N + NCLIN
         IS = ISTATE(J)
C
         IF (IS.GE.0) THEN
            FEASJ = FEATOL(J)
C
            IF (J.LE.N) THEN
               CON = X(J)
            ELSE
               CON = AX(J-N)
            END IF
C
C           Check for constraint violations.
C
            IF (BL(J).GT.BIGLOW) THEN
               RES = BL(J) - CON
               IF (RES.GT.FEASJ) THEN
                  NVIOL = NVIOL + 1
                  GO TO 20
               END IF
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               RES = BU(J) - CON
               IF (RES.LT.(-FEASJ)) THEN
                  NVIOL = NVIOL + 1
                  RES = -RES
                  GO TO 20
               END IF
            END IF
C
C           this constraint is satisfied,  but count a large residual
C           as a violation if the constraint is in the working set.
C
            RES = ZERO
C
            IF (IS.EQ.1) THEN
               RES = ABS(BL(J)-CON)
C
            ELSE IF (IS.EQ.2) THEN
               RES = ABS(BU(J)-CON)
C
            ELSE IF (IS.EQ.3) THEN
               RES = ABS(BU(J)-CON)
            END IF
C
            IF (RES.GT.FEASJ) NVIOL = NVIOL + 1
C
   20       IF (RES.GT.ERRMAX) THEN
               JMAX = J
               ERRMAX = RES
            END IF
         END IF
   40 CONTINUE
C
      RETURN
C
C     End of  E04MFQ.  (CMFEAS)
C
      END
