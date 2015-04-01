      SUBROUTINE E04NCU(COLD,VERTEX,NCLIN,NCTOTL,NACTIV,NARTIF,NFREE,N,
     *                  LDA,ISTATE,KACTIV,BIGBND,TOLACT,A,AX,BL,BU,X,WX,
     *                  WORK)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-720 (DEC 1989).
C     MARK 14C REVISED. IER-890 (NOV 1990).
C     MARK 16A REVISED. IER-996 (JUN 1993).
C     MARK 17 REVISED. IER-1584 (JUN 1995).
C
C     ******************************************************************
C     E04NCU  computes the quantities  ISTATE (optionally), KACTIV,
C     NACTIV, nz and NFREE  associated with the working set at X.
C     The computation depends upon the value of the input parameter
C     COLD,  as follows...
C
C     COLD = TRUE.  An initial working set will be selected. First,
C                   nearly-satisfied or violated bounds are added.
C                   Next,  general linear constraints are added that
C                   have small residuals.
C
C     COLD = FALSE. The quantities KACTIV, NACTIV, nz and NFREE are
C                   computed from ISTATE,  specified by the user.
C
C     Values of ISTATE(j)....
C
C        - 2         - 1         0           1          2         3
C     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCU dated 14-May-1992.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, TOLACT
      INTEGER           LDA, N, NACTIV, NARTIF, NCLIN, NCTOTL, NFREE
      LOGICAL           COLD, VERTEX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL),
     *                  WORK(N), WX(N), X(N)
      INTEGER           ISTATE(NCTOTL), KACTIV(N)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  B1, B2, BIGLOW, BIGUPP, COLMIN, COLSIZ, FLMAX,
     *                  RESIDL, RESL, RESMIN, RESU, TOOBIG
      INTEGER           I, IMIN, IS, J, JMIN, K, NFIXED
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Executable Statements ..
C
      FLMAX = WMACH(7)
      BIGLOW = -BIGBND
      BIGUPP = BIGBND
C
C     ------------------------------------------------------------------
C     Move the variables inside their bounds.
C     ------------------------------------------------------------------
C
      DO 20 J = 1, N
         B1 = BL(J)
         B2 = BU(J)
C
         IF (B1.GT.BIGLOW) THEN
            IF (X(J).LT.B1) X(J) = B1
         END IF
C
         IF (B2.LT.BIGUPP) THEN
            IF (X(J).GT.B2) X(J) = B2
         END IF
   20 CONTINUE
C
      CALL DCOPY(N,X,1,WX,1)
C
      NFIXED = 0
      NACTIV = 0
      NARTIF = 0
C
C     If a cold start is being made, initialize  ISTATE.
C     If  BL(j) = BU(j),  set  ISTATE(j)=3  for all variables and linear
C     constraints.
C
      IF (COLD) THEN
         DO 40 J = 1, NCTOTL
            ISTATE(J) = 0
            IF (BL(J).EQ.BU(J)) ISTATE(J) = 3
   40    CONTINUE
      ELSE
         DO 60 J = 1, NCTOTL
            IF (BL(J).EQ.BU(J)) THEN
               ISTATE(J) = 3
            ELSE IF (ISTATE(J).GE.3 .OR. ISTATE(J).LT.0) THEN
               ISTATE(J) = 0
            END IF
            IF (BL(J).LE.BIGLOW .AND. BU(J).GE.BIGUPP) ISTATE(J) = 0
            IF (BL(J).LE.BIGLOW .AND. ISTATE(J).EQ.1) ISTATE(J) = 0
            IF (BU(J).GE.BIGUPP .AND. ISTATE(J).EQ.2) ISTATE(J) = 0
   60    CONTINUE
      END IF
C
C     Initialize NFIXED, NFREE and KACTIV.
C     Ensure that the number of bounds and general constraints in the
C     working set does not exceed N.
C
      DO 80 J = 1, NCTOTL
         IF (NFIXED+NACTIV.EQ.N) ISTATE(J) = 0
         IF (ISTATE(J).GT.0) THEN
            IF (J.LE.N) THEN
               NFIXED = NFIXED + 1
               IF (ISTATE(J).EQ.1) WX(J) = BL(J)
               IF (ISTATE(J).GE.2) WX(J) = BU(J)
            ELSE
               NACTIV = NACTIV + 1
               KACTIV(NACTIV) = J - N
            END IF
         END IF
   80 CONTINUE
C
C     ------------------------------------------------------------------
C     If a cold start is required,  attempt to add as many
C     constraints as possible to the working set.
C     ------------------------------------------------------------------
      IF (COLD) THEN
C
C        See if any bounds are violated or nearly satisfied.
C        If so,  add these bounds to the working set and set the
C        variables exactly on their bounds.
C
         J = N
C        +       WHILE (J .GE. 1  .AND.  NFIXED + NACTIV .LT. N) DO
  100    IF (J.GE.1 .AND. NFIXED+NACTIV.LT.N) THEN
            IF (ISTATE(J).EQ.0) THEN
               B1 = BL(J)
               B2 = BU(J)
               IS = 0
               IF (B1.GT.BIGLOW) THEN
                  IF (WX(J)-B1.LE.(ONE+ABS(B1))*TOLACT) IS = 1
               END IF
               IF (B2.LT.BIGUPP) THEN
                  IF (B2-WX(J).LE.(ONE+ABS(B2))*TOLACT) IS = 2
               END IF
               IF (IS.GT.0) THEN
                  ISTATE(J) = IS
                  IF (IS.EQ.1) WX(J) = B1
                  IF (IS.EQ.2) WX(J) = B2
                  NFIXED = NFIXED + 1
               END IF
            END IF
            J = J - 1
            GO TO 100
C           +       END WHILE
         END IF
C
C        ---------------------------------------------------------------
C        The following loop finds the linear constraint (if any) with
C        smallest residual less than or equal to TOLACT  and adds it
C        to the working set.  This is repeated until the working set
C        is complete or all the remaining residuals are too large.
C        ---------------------------------------------------------------
C        First, compute the residuals for all the constraints not in the
C        working set.
C
         IF (NCLIN.GT.0 .AND. NACTIV+NFIXED.LT.N) THEN
            DO 120 I = 1, NCLIN
               IF (ISTATE(N+I).LE.0) AX(I) = DDOT(N,A(I,1),LDA,WX,1)
  120       CONTINUE
C
            IS = 1
            TOOBIG = TOLACT + TOLACT
C
C           + WHILE (IS .GT. 0  .AND.  NFIXED + NACTIV .LT. N) DO
  140       IF (IS.GT.0 .AND. NFIXED+NACTIV.LT.N) THEN
               IS = 0
               RESMIN = TOLACT
C
               DO 160 I = 1, NCLIN
                  J = N + I
                  IF (ISTATE(J).EQ.0) THEN
                     B1 = BL(J)
                     B2 = BU(J)
                     RESL = TOOBIG
                     RESU = TOOBIG
                     IF (B1.GT.BIGLOW) RESL = ABS(AX(I)-B1)/(ONE+ABS(B1)
     *                                        )
                     IF (B2.LT.BIGUPP) RESU = ABS(AX(I)-B2)/(ONE+ABS(B2)
     *                                        )
                     RESIDL = MIN(RESL,RESU)
                     IF (RESIDL.LT.RESMIN) THEN
                        RESMIN = RESIDL
                        IMIN = I
                        IS = 1
                        IF (RESL.GT.RESU) IS = 2
                     END IF
                  END IF
  160          CONTINUE
C
               IF (IS.GT.0) THEN
                  NACTIV = NACTIV + 1
                  KACTIV(NACTIV) = IMIN
                  J = N + IMIN
                  ISTATE(J) = IS
               END IF
               GO TO 140
C              +          END WHILE
            END IF
         END IF
      END IF
C
      IF (VERTEX .AND. NACTIV+NFIXED.LT.N) THEN
C        ---------------------------------------------------------------
C        Find an initial vertex by temporarily fixing some variables.
C        ---------------------------------------------------------------
C        Compute lengths of columns of selected linear constraints
C        (just the ones corresponding to variables eligible to be
C        temporarily fixed).
C
         DO 200 J = 1, N
            IF (ISTATE(J).EQ.0) THEN
               COLSIZ = ZERO
               DO 180 K = 1, NCLIN
                  IF (ISTATE(N+K).GT.0) COLSIZ = COLSIZ + ABS(A(K,J))
  180          CONTINUE
               WORK(J) = COLSIZ
            END IF
  200    CONTINUE
C
C        Find the  NARTIF  smallest such columns.
C        This is an expensive loop.  Later we can replace it by a
C        4-pass process (say), accepting the first col that is within
C        t  of  COLMIN, where  t = 0.0, 0.001, 0.01, 0.1 (say).
C        (This comment written in 1980).
C
C        +       WHILE (NFIXED + NACTIV .LT. N) DO
  220    IF (NFIXED+NACTIV.LT.N) THEN
            COLMIN = FLMAX
            DO 240 J = 1, N
               IF (ISTATE(J).EQ.0) THEN
                  IF (NCLIN.EQ.0) GO TO 260
                  COLSIZ = WORK(J)
                  IF (COLMIN.GT.COLSIZ) THEN
                     COLMIN = COLSIZ
                     JMIN = J
                  END IF
               END IF
  240       CONTINUE
            J = JMIN
  260       ISTATE(J) = 4
            NARTIF = NARTIF + 1
            NFIXED = NFIXED + 1
            GO TO 220
C           +       END WHILE
         END IF
      END IF
C
      NFREE = N - NFIXED
C
      RETURN
C
C
C     End of  E04NCU. (LSCRSH)
C
      END
