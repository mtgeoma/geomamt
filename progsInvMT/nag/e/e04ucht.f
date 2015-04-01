      SUBROUTINE E04UCH(FIRSTV,NEGSTP,BIGALF,BIGBND,PNORM,JADD1,JADD2,
     *                  PALFA1,PALFA2,ISTATE,N,NCTOTL,ANORM,AP,AX,BL,BU,
     *                  FEATOL,P,X)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1079 (JUL 1993).
C     MARK 17 REVISED. IER-1601 (JUN 1995).
C
C     ******************************************************************
C     E04UCH  finds steps PALFA1, PALFA2 such that
C        X + PALFA1*P  reaches a linear constraint that is currently not
C                      in the working set but is satisfied.
C        X + PALFA2*P  reaches a linear constraint that is currently not
C                      in the working set but is violated.
C     The constraints are perturbed by an amount FEATOL, so that PALFA1
C     is slightly larger than it should be,  and PALFA2 is slightly
C     smaller than it should be.  This gives some leeway later when the
C     exact steps are computed by E04UCG.
C
C     Constraints in the working set are ignored  (ISTATE(j) .GE. 1).
C
C     If NEGSTP is true, the search direction will be taken to be  - P.
C
C
C     Values of ISTATE(j)....
C
C        - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     The values  -2  and  -1  do not occur once a feasible point has
C     been found.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original Fortran 66 version written  May 1980.
C     This version of E04UCH dated 26-June-1986.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGALF, BIGBND, PALFA1, PALFA2, PNORM
      INTEGER           JADD1, JADD2, N, NCTOTL
      LOGICAL           FIRSTV, NEGSTP
C     .. Array Arguments ..
      DOUBLE PRECISION  ANORM(*), AP(*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  FEATOL(NCTOTL), P(N), X(N)
      INTEGER           ISTATE(NCTOTL)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSATP, ATP, ATX, RES, ROWNRM
      INTEGER           I, J, JS
      LOGICAL           LASTV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /CE04NB/EPSPT3, EPSPT5, EPSPT8, EPSPT9
C     .. Executable Statements ..
C
      LASTV = .NOT. FIRSTV
      JADD1 = 0
      JADD2 = 0
      PALFA1 = BIGALF
C
      PALFA2 = ZERO
      IF (FIRSTV) PALFA2 = BIGALF
C
      DO 20 J = 1, NCTOTL
         JS = ISTATE(J)
         IF (JS.LE.0) THEN
            IF (J.LE.N) THEN
               ATX = X(J)
               ATP = P(J)
               ROWNRM = ONE
            ELSE
               I = J - N
               ATX = AX(I)
               ATP = AP(I)
               ROWNRM = ONE + ANORM(I)
            END IF
            IF (NEGSTP) ATP = -ATP
C
            IF (ABS(ATP).LE.EPSPT9*ROWNRM*PNORM) THEN
C
C              This constraint appears to be constant along P.  It is
C              not used to compute the step.  Give the residual a value
C              that can be spotted in the debug output.
C
               RES = -ONE
            ELSE IF (ATP.LE.ZERO .AND. JS.NE.-2) THEN
C              ---------------------------------------------------------
C              a'x  is decreasing and the lower bound is not violated.
C              ---------------------------------------------------------
C              First test for smaller PALFA1.
C
               ABSATP = -ATP
               IF (BL(J).GT.(-BIGBND)) THEN
                  RES = ATX - BL(J) + FEATOL(J)
                  IF (BIGALF*ABSATP.GT.ABS(RES)) THEN
                     IF (PALFA1*ABSATP.GT.RES) THEN
                        PALFA1 = RES/ABSATP
                        JADD1 = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-1) THEN
C
C                 The upper bound is violated.  Test for either larger
C                 or smaller PALFA2, depending on the value of FIRSTV.
C
                  RES = ATX - BU(J) - FEATOL(J)
                  IF (BIGALF*ABSATP.GT.ABS(RES)) THEN
                     IF (FIRSTV .AND. PALFA2*ABSATP.GT.RES .OR.
     *                   LASTV .AND. PALFA2*ABSATP.LT.RES) THEN
                        PALFA2 = RES/ABSATP
                        JADD2 = J
                     END IF
                  END IF
               END IF
            ELSE IF (ATP.GT.ZERO .AND. JS.NE.-1) THEN
C              ---------------------------------------------------------
C              a'x  is increasing and the upper bound is not violated.
C              ---------------------------------------------------------
C              Test for smaller PALFA1.
C
               IF (BU(J).LT.BIGBND) THEN
                  RES = BU(J) - ATX + FEATOL(J)
                  IF (BIGALF*ATP.GT.ABS(RES)) THEN
                     IF (PALFA1*ATP.GT.RES) THEN
                        PALFA1 = RES/ATP
                        JADD1 = J
                     END IF
                  END IF
               END IF
C
               IF (JS.EQ.-2) THEN
C
C                 The lower bound is violated.  Test for a new PALFA2.
C
                  RES = BL(J) - ATX - FEATOL(J)
                  IF (BIGALF*ATP.GT.ABS(RES)) THEN
                     IF (FIRSTV .AND. PALFA2*ATP.GT.RES .OR. LASTV .AND.
     *                   PALFA2*ATP.LT.RES) THEN
                        PALFA2 = RES/ATP
                        JADD2 = J
                     END IF
                  END IF
               END IF
            END IF
C
         END IF
   20 CONTINUE
C
      RETURN
C
C
C     End of  E04UCH. (CMALF1)
C
      END
