      SUBROUTINE E04UCW(N,NCLIN,NCNLN,ISTATE,BIGBND,CVNORM,ERRMAX,JMAX,
     *                  NVIOL,AX,BL,BU,C,FEATOL,X,WORK)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1090 (JUL 1993).
C     MARK 17 REVISED. IER-1610 (JUN 1995).
C
C     ******************************************************************
C     E04UCW  computes the following...
C     (1)  The number of constraints that are violated by more
C          than  FEATOL  and the 2-norm of the constraint violations.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version      April    1984.
C     This version of  E04UCW  dated  16-October-1985.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, CVNORM, ERRMAX
      INTEGER           JMAX, N, NCLIN, NCNLN, NVIOL
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(*), BL(N+NCLIN+NCNLN), BU(N+NCLIN+NCNLN),
     *                  C(*), FEATOL(N+NCLIN+NCNLN),
     *                  WORK(N+NCLIN+NCNLN), X(N)
      INTEGER           ISTATE(N+NCLIN+NCNLN)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CON, FEASJ, RES, TOLJ
      INTEGER           IS, J
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
C     Compute NVIOL, the number of constraints violated by more than
C     FEATOL,  and CVNORM,  the 2-norm of the constraint
C     violations and residuals of the constraints in the QP working set.
C     ==================================================================
      NVIOL = 0
C
      DO 40 J = 1, N + NCLIN + NCNLN
         FEASJ = FEATOL(J)
         RES = ZERO
C
         IF (J.LE.N+NCLIN) THEN
C
C           Bound or general linear constraint.
C
            IF (J.LE.N) THEN
               CON = X(J)
            ELSE
               CON = AX(J-N)
            END IF
C
            TOLJ = FEASJ
         ELSE
C
C           Nonlinear constraint.
C
            CON = C(J-N-NCLIN)
            TOLJ = ZERO
         END IF
C
C        Check for constraint violations.
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
C        This constraint is satisfied,  but count the residual as a
C        violation if the constraint is in the working set.
C
         IS = ISTATE(J)
C
         IF (IS.EQ.0) THEN
            RES = ZERO
         ELSE IF (IS.EQ.1 .OR. IS.LE.-2) THEN
            RES = BL(J) - CON
         ELSE IF (IS.GE.2 .OR. IS.EQ.-1) THEN
            RES = BU(J) - CON
         END IF
C
         IF (ABS(RES).GT.FEASJ) NVIOL = NVIOL + 1
C
C        Set the array of violations.
C
   20    WORK(J) = RES
   40 CONTINUE
C
      JMAX = IDAMAX(N+NCLIN+NCNLN,WORK,1)
      ERRMAX = ABS(WORK(JMAX))
C
      CVNORM = DNRM2(N+NCLIN+NCNLN,WORK,1)
C
      RETURN
C
C
C     End of  E04UCW. (NPFEAS)
C
      END
