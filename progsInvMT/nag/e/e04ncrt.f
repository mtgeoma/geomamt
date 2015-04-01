      SUBROUTINE E04NCR(N,NCLIN,ISTATE,BIGBND,CVNORM,ERRMAX,JMAX,NVIOL,
     *                  AX,BL,BU,FEATOL,X,WORK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-1069 (JUL 1993).
C     MARK 17 REVISED. IER-1581 (JUN 1995).
C
C     ******************************************************************
C     E04NCR  computes the following...
C     (1)  The number of constraints that are violated by more
C          than  FEATOL  and the 2-norm of the constraint violations.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version      April    1984.
C     This version of  E04NCR  dated  17-October-1985.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CVNORM, ERRMAX
      INTEGER           JMAX, N, NCLIN, NVIOL
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATOL(N+NCLIN), WORK(N+NCLIN), X(N)
      INTEGER           ISTATE(N+NCLIN)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CON, FEASJ, RES, TOLJ
      INTEGER           I, IS, J
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DNRM2, IDAMAX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ==================================================================
C     Compute NVIOL,  the number of constraints violated by more than
C     FEATOL,  and CVNORM,  the 2-norm of the constraint violations and
C     residuals of the constraints in the working set.
C     ==================================================================
      NVIOL = 0
C
      DO 40 J = 1, N + NCLIN
         FEASJ = FEATOL(J)
         IS = ISTATE(J)
         RES = ZERO
C
         IF (IS.GE.0 .AND. IS.LT.4) THEN
            IF (J.LE.N) THEN
               CON = X(J)
            ELSE
               I = J - N
               CON = AX(I)
            END IF
C
            TOLJ = FEASJ
C
C           Check for constraint violations.
C
            IF (BL(J).GT.BIGLOW) THEN
               RES = BL(J) - CON
               IF (RES.GT.FEASJ) NVIOL = NVIOL + 1
               IF (RES.GT.TOLJ) GO TO 20
            END IF
C
            IF (BU(J).LT.BIGUPP) THEN
               RES = BU(J) - CON
               IF (RES.LT.(-FEASJ)) NVIOL = NVIOL + 1
               IF (RES.LT.(-TOLJ)) GO TO 20
            END IF
C
C           This constraint is satisfied,  but count the residual as a
C           violation if the constraint is in the working set.
C
            IF (IS.LE.0) RES = ZERO
            IF (IS.EQ.1) RES = BL(J) - CON
            IF (IS.GE.2) RES = BU(J) - CON
            IF (ABS(RES).GT.FEASJ) NVIOL = NVIOL + 1
         END IF
   20    WORK(J) = RES
   40 CONTINUE
C
      JMAX = IDAMAX(N+NCLIN,WORK,1)
      ERRMAX = ABS(WORK(JMAX))
C
      CVNORM = DNRM2(N+NCLIN,WORK,1)
C
      RETURN
C
C
C     End of  E04NCR. (LSFEAS)
C
      END
