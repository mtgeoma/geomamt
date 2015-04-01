      SUBROUTINE E04MFJ(ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NZ,N,LDQ,LDA,
     *                  LDT,ISTATE,KACTIV,KX,JMAX,ERRMAX,XNORM,A,AX,BL,
     *                  BU,FEATOL,T,X,Q,P,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1557 (JUN 1995).
C
C     ******************************************************************
C     E04MFJ  computes the point on a working set that is closest in the
C     least-squares sense to the input vector X.
C
C     If the computed point gives a row error of more than the
C     feasibility tolerance, an extra step of iterative refinement is
C     used.  If  X  is still infeasible,  the logical variable ROWERR
C     is set.
C
C     Original version derived from LSSETX January-1987.
C     This version of  E04MFJ  dated   5-Jul-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      INTEGER           NTRY
      PARAMETER         (NTRY=5)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERRMAX, XNORM
      INTEGER           JMAX, LDA, LDQ, LDT, N, NACTIV, NCLIN, NFREE, NZ
      LOGICAL           ROWERR, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(N+NCLIN), BU(N+NCLIN),
     *                  FEATOL(N+NCLIN), P(N), Q(LDQ,*), T(LDT,*),
     *                  WORK(N), X(N)
      INTEGER           ISTATE(N+NCLIN), KACTIV(N), KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
      INTEGER           I, IS, J, K, KTRY
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DDOT, DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRSV, E04NBW, F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
C     ------------------------------------------------------------------
C     Move  x  onto the simple bounds in the working set.
C     ------------------------------------------------------------------
      DO 20 K = NFREE + 1, N
         J = KX(K)
         IS = ISTATE(J)
         BND = BL(J)
         IF (IS.GE.2) BND = BU(J)
         IF (IS.NE.4) X(J) = BND
   20 CONTINUE
C
C     ------------------------------------------------------------------
C     Move  x  onto the general constraints in the working set.
C     ntry  attempts are made to get acceptable row errors.
C     ------------------------------------------------------------------
      KTRY = 1
      JMAX = 1
      ERRMAX = ZERO
C
C     repeat
   40 IF (NACTIV.GT.0) THEN
C
C        Set work = residuals for constraints in the working set.
C        Solve for P, the smallest correction to x that gives a point
C        on the constraints in the working set.  Define  P = Y*(py),
C        where  py  solves the triangular system  T*(py) = residuals.
C
         DO 60 I = 1, NACTIV
            K = KACTIV(I)
            J = N + K
            BND = BL(J)
            IF (ISTATE(J).EQ.2) BND = BU(J)
            WORK(NACTIV-I+1) = BND - DDOT(N,A(K,1),LDA,X,1)
   60    CONTINUE
C
         CALL DTRSV('U','N','N',NACTIV,T(1,NZ+1),LDT,WORK,1)
         CALL F06FBF(N,ZERO,P,1)
         CALL DCOPY(NACTIV,WORK,1,P(NZ+1),1)
C
         CALL E04NBW(2,N,NZ,NFREE,LDQ,UNITQ,KX,P,Q,WORK)
         CALL DAXPY(N,ONE,P,1,X,1)
      END IF
C
C     ---------------------------------------------------------------
C     Compute the 2-norm of  x.
C     Initialize  Ax  for all the general constraints.
C     ---------------------------------------------------------------
      XNORM = DNRM2(N,X,1)
      IF (NCLIN.GT.0) CALL DGEMV('N',NCLIN,N,ONE,A,LDA,X,1,ZERO,AX,1)
C
C     ---------------------------------------------------------------
C     Check the row residuals.
C     ---------------------------------------------------------------
      IF (NACTIV.GT.0) THEN
         DO 80 K = 1, NACTIV
            I = KACTIV(K)
            J = N + I
            IS = ISTATE(J)
            IF (IS.EQ.1) WORK(K) = BL(J) - AX(I)
            IF (IS.GE.2) WORK(K) = BU(J) - AX(I)
   80    CONTINUE
C
         JMAX = IDAMAX(NACTIV,WORK,1)
         ERRMAX = ABS(WORK(JMAX))
      END IF
C
      KTRY = KTRY + 1
C     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      IF ( .NOT. (ERRMAX.LE.FEATOL(JMAX) .OR. KTRY.GT.NTRY)) GO TO 40
C
      ROWERR = ERRMAX .GT. FEATOL(JMAX)
C
      RETURN
C
C     End of  E04MFJ.  (CMSETX)
C
      END
