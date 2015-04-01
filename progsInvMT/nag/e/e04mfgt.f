      SUBROUTINE E04MFG(NFREE,LDA,N,NCLIN,NCTOTL,NACTIV,ISTATE,KACTIV,
     *                  KX,A,BL,BU,C,CLAMDA,FEATOL,R,RLAMDA,X)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     ==================================================================
C     E04MFG  creates the expanded Lagrange multiplier vector CLAMDA
C     and resets ISTATE for the printed solution.
C
C     KACTIV holds the general constraint indices in the order in which
C     they were added to the working set. The reverse ordering is used
C     for T since new rows are added at the front of T.
C
C     E04MFG is called by...
C        E04MFF and E04NFF just before exiting.
C
C     Original Fortran 77 version written  05-May-93.
C     This version of  E04MFG  dated  05-May-93.
C     ==================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           LDA, N, NACTIV, NCLIN, NCTOTL, NFREE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(NCTOTL), BU(NCTOTL), C(*),
     *                  CLAMDA(NCTOTL), FEATOL(NCTOTL), R(NCTOTL),
     *                  RLAMDA(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, RJ, RLAM, SLK1, SLK2, TOL
      INTEGER           I, IS, J, K, NFIXED, NPLIN, NZ
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      NPLIN = N + NCLIN
      NZ = NFREE - NACTIV
C
C     Expand multipliers for bounds, linear and nonlinear constraints
C     into the  CLAMDA  array.
C
      CALL F06FBF(NCTOTL,ZERO,CLAMDA,1)
      DO 20 K = 1, NACTIV + NFIXED
         IF (K.LE.NACTIV) THEN
            J = KACTIV(K) + N
            RLAM = RLAMDA(NACTIV-K+1)
         ELSE
            J = KX(NZ+K)
            RLAM = RLAMDA(K)
         END IF
         CLAMDA(J) = RLAM
   20 CONTINUE
C
C     Reset ISTATE if necessary.
C
      DO 40 J = 1, NCTOTL
         B1 = BL(J)
         B2 = BU(J)
C
         IF (J.LE.N) THEN
            RJ = X(J)
         ELSE IF (J.LE.NPLIN) THEN
            I = J - N
            RJ = DDOT(N,A(I,1),LDA,X,1)
         ELSE
            I = J - NPLIN
            RJ = C(I)
         END IF
C
         IS = ISTATE(J)
         SLK1 = RJ - B1
         SLK2 = B2 - RJ
         TOL = FEATOL(J)
         IF (SLK1.LT.-TOL) IS = -2
         IF (SLK2.LT.-TOL) IS = -1
         IF (IS.EQ.1 .AND. SLK1.GT.TOL) IS = 0
         IF (IS.EQ.2 .AND. SLK2.GT.TOL) IS = 0
         ISTATE(J) = IS
         R(J) = RJ
   40 CONTINUE
C
C     End of E04MFG.  (CMWRAP in CMOPTSUBS.F)
C
      RETURN
      END
