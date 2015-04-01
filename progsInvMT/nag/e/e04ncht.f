      SUBROUTINE E04NCH(LINOBJ,ROWERR,UNITQ,NCLIN,NACTIV,NFREE,NRANK,NZ,
     *                  N,NCTOTL,LDZY,LDA,LDR,LDT,ISTATE,KACTIV,KX,JMAX,
     *                  ERRMAX,CTX,XNORM,A,AX,BL,BU,CQ,RES,RES0,FEATOL,
     *                  R,T,X,ZY,P,WORK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1061 (JUL 1993).
C     MARK 17 REVISED. IER-1573 (JUN 1995).
C
C     ******************************************************************
C     E04NCH  computes the point on a working set that is closest to the
C     input vector  x  (in the least-squares sense).  The norm of  x,
C     the transformed residual vector  Pr - RQ'x,  and the constraint
C     values Ax  are also initialized.
C
C     If the computed point gives a row error of more than the
C     feasibility tolerance, an extra step of iterative refinement is
C     used.  If  x  is still infeasible,  the logical variable  ROWERR
C     is set.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version dated 29-December-1985.
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           NTRY
      PARAMETER         (NTRY=5)
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CTX, ERRMAX, XNORM
      INTEGER           JMAX, LDA, LDR, LDT, LDZY, N, NACTIV, NCLIN,
     *                  NCTOTL, NFREE, NRANK, NZ
      LOGICAL           LINOBJ, ROWERR, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AX(*), BL(NCTOTL), BU(NCTOTL), CQ(*),
     *                  FEATOL(NCTOTL), P(N), R(LDR,*), RES(*), RES0(*),
     *                  T(LDT,*), WORK(NCTOTL), X(N), ZY(LDZY,*)
      INTEGER           ISTATE(NCTOTL), KACTIV(N), KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BND
      INTEGER           I, IS, J, K, KTRY
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           IDAMAX
      EXTERNAL          DDOT, DNRM2, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRMV, E04NBT, E04NBW,
     *                  F06FBF
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
C     We shall make  NTRY  tries at getting acceptable row errors.
C     ------------------------------------------------------------------
      KTRY = 1
      JMAX = 1
      ERRMAX = ZERO
C
C     REPEAT
   40 IF (NACTIV.GT.0) THEN
C
C        Set  work = residuals for constraints in the working set.
C        Solve for p, the smallest correction to x that gives a point
C        on the constraints in the working set.  Define  p = Y*(py),
C        where  py  solves the triangular system  T*(py) = residuals.
C
         DO 60 I = 1, NACTIV
            K = KACTIV(I)
            J = N + K
            BND = BL(J)
            IF (ISTATE(J).EQ.2) BND = BU(J)
            WORK(I) = BND - DDOT(N,A(K,1),LDA,X,1)
   60    CONTINUE
C
         CALL E04NBT(1,LDT,NACTIV,T(1,NZ+1),WORK)
         CALL F06FBF(N,ZERO,P,1)
         CALL DCOPY(NACTIV,WORK,1,P(NZ+1),1)
C
         CALL E04NBW(2,N,NZ,NFREE,LDZY,UNITQ,KX,P,ZY,WORK)
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
C     UNTIL    (ERRMAX .LE. FEATOL(JMAX) .OR. KTRY .GT. NTRY
      IF ( .NOT. (ERRMAX.LE.FEATOL(JMAX) .OR. KTRY.GT.NTRY)) GO TO 40
C
      ROWERR = ERRMAX .GT. FEATOL(JMAX)
C
C     ==================================================================
C     Compute the linear objective value  c'x  and the transformed
C     residual  Pr  -  RQ'x = RES0  -  RQ'x.
C     ==================================================================
      IF (NRANK.GT.0 .OR. LINOBJ) THEN
         CALL DCOPY(N,X,1,P,1)
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,P,ZY,WORK)
      END IF
C
      CTX = ZERO
      IF (LINOBJ) CTX = DDOT(N,CQ,1,P,1)
C
      IF (NRANK.GT.0) THEN
         CALL DTRMV('U','N','N',NRANK,R,LDR,P,1)
         IF (NRANK.LT.N) CALL DGEMV('N',NRANK,N-NRANK,ONE,R(1,NRANK+1),
     *                              LDR,P(NRANK+1),1,ONE,P,1)
         CALL DCOPY(NRANK,RES0,1,RES,1)
         CALL DAXPY(NRANK,-ONE,P,1,RES,1)
      END IF
C
      RETURN
C
C
C     End of  E04NCH. (LSSETX)
C
      END
