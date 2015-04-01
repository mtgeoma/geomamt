      SUBROUTINE E04UCS(COLD,N,NCLIN,NCNLN,NCTOTL,NACTIV,NFREE,NZ,
     *                  ISTATE,KACTIV,BIGBND,TOLACT,BL,BU,C)
C     MARK 13 RE-ISSUE.  NAG COPYRIGHT 1988.
C     MARK 16 REVISED. IER-1087 (JUL 1993).
C     MARK 17 REVISED. IER-1607 (JUN 1995).
C
C     ******************************************************************
C     E04UCS  adds indices of nonlinear constraints to the initial
C     working set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   14-February 1985.
C     This version of  E04UCS  dated 14-November-1985.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, TOLACT
      INTEGER           N, NACTIV, NCLIN, NCNLN, NCTOTL, NFREE, NZ
      LOGICAL           COLD
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(NCTOTL), BU(NCTOTL), C(*)
      INTEGER           ISTATE(NCTOTL), KACTIV(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, BIGLOW, BIGUPP, CMIN, RES, RESL, RESU,
     *                  TOOBIG
      INTEGER           I, IMIN, IS, J, NFIXED, NPLIN
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
C
      NFIXED = N - NFREE
      NPLIN = N + NCLIN
C
C     If a cold start is being made, initialize the status of the QP
C     working set.  First,  if  BL(j) = BU(j),  set ISTATE(j)=3.
C
      IF (COLD) THEN
         DO 20 J = NPLIN + 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J).EQ.BU(J)) ISTATE(J) = 3
   20    CONTINUE
      END IF
C
C     Increment NACTIV and KACTIV.
C     Ensure that the number of bounds and general constraints in the
C     QP  working set does not exceed N.
C
      DO 40 J = NPLIN + 1, NCTOTL
         IF (NFIXED+NACTIV.EQ.N) ISTATE(J) = 0
         IF (ISTATE(J).GT.0) THEN
            NACTIV = NACTIV + 1
            KACTIV(NACTIV) = J - N
         END IF
   40 CONTINUE
C
      IF (COLD) THEN
C
C        ---------------------------------------------------------------
C        If a cold start is required, an attempt is made to add as many
C        nonlinear constraints as possible to the working set.
C        ---------------------------------------------------------------
C        The following loop finds the most violated constraint.  If
C        there is room in KACTIV, it will be added to the working set
C        and the process will be repeated.
C
C
         IS = 1
         BIGLOW = -BIGBND
         BIGUPP = BIGBND
         TOOBIG = TOLACT + TOLACT
C
C        while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
   60    IF (IS.GT.0 .AND. NFIXED+NACTIV.LT.N) THEN
            IS = 0
            CMIN = TOLACT
C
            DO 80 I = 1, NCNLN
               J = NPLIN + I
               IF (ISTATE(J).EQ.0) THEN
                  B1 = BL(J)
                  B2 = BU(J)
                  RESL = TOOBIG
                  RESU = TOOBIG
                  IF (B1.GT.BIGLOW) RESL = ABS(C(I)-B1)/(ONE+ABS(B1))
                  IF (B2.LT.BIGUPP) RESU = ABS(C(I)-B2)/(ONE+ABS(B2))
                  RES = MIN(RESL,RESU)
                  IF (RES.LT.CMIN) THEN
                     CMIN = RES
                     IMIN = I
                     IS = 1
                     IF (RESL.GT.RESU) IS = 2
                  END IF
               END IF
   80       CONTINUE
C
            IF (IS.GT.0) THEN
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = NCLIN + IMIN
               J = NPLIN + IMIN
               ISTATE(J) = IS
            END IF
            GO TO 60
C           end while
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     An initial working set has now been selected.
C     ------------------------------------------------------------------
      NZ = NFREE - NACTIV
C
      RETURN
C
C
C     End of  E04UCS. (NPCRSH)
C
      END
