      SUBROUTINE E04NCP(PRBTYP,LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,
     *                  LDA,LDZY,LDR,NRANK,NZ,NRZ,ISTATE,KX,BIGBND,
     *                  TOLRNK,NUMINF,SUMINF,BL,BU,A,RES,FEATOL,GQ,CQ,R,
     *                  X,WTINF,ZY,WRK)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 16 REVISED. IER-1067 (JUL 1993).
C     MARK 17 REVISED. IER-1579 (JUN 1995).
C
C     ******************************************************************
C     E04NCP  finds the number and weighted sum of infeasibilities for
C     the bounds and linear constraints.   An appropriate transformed
C     gradient vector is returned in  GQ.
C
C     Positive values of  ISTATE(j)  will not be altered.  These mean
C     the following...
C
C               1             2           3
C           a'x = bl      a'x = bu     bl = bu
C
C     Other values of  ISTATE(j)  will be reset as follows...
C           a'x lt bl     a'x gt bu     a'x free
C              - 2           - 1           0
C
C     If  x  is feasible,  E04NCP computes the vector Q(free)'g(free),
C     where  g  is the gradient of the the sum of squares plus the
C     linear term.  The matrix Q is of the form
C                    ( Q(free)  0       ),
C                    (   0      I(fixed))
C     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
C     the matrix of constraints in the working set.  The transformed
C     gradients are stored in GQ.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     This version of E04NCP dated 14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGBND, SUMINF, TOLRNK
      INTEGER           LDA, LDR, LDZY, N, NCLIN, NFREE, NRANK, NRZ,
     *                  NUMINF, NZ
      LOGICAL           LINOBJ, SINGLR, UNITGZ, UNITQ
      CHARACTER*2       PRBTYP
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), BL(*), BU(*), CQ(*), FEATOL(*), GQ(N),
     *                  R(LDR,*), RES(*), WRK(N), WTINF(*), X(N),
     *                  ZY(LDZY,*)
      INTEGER           ISTATE(*), KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGLOW, BIGUPP, CTX, FEASJ, ROWNRM, S, WEIGHT
      INTEGER           J, K
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      INTEGER           F06KLF
      EXTERNAL          DDOT, DNRM2, F06KLF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DSCAL, DTRMV, E04NBW,
     *                  F06FBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      BIGUPP = BIGBND
      BIGLOW = -BIGBND
C
      NUMINF = 0
      SUMINF = ZERO
      CALL F06FBF(N,ZERO,GQ,1)
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
               GQ(J) = WEIGHT
            ELSE
               CALL DAXPY(N,WEIGHT,A(K,1),LDA,GQ,1)
            END IF
         END IF
   40 CONTINUE
C
C     ------------------------------------------------------------------
C     Install  GQ,  the transformed gradient.
C     ------------------------------------------------------------------
      SINGLR = .FALSE.
      UNITGZ = .TRUE.
C
      IF (NUMINF.GT.0) THEN
         CALL E04NBW(6,N,NZ,NFREE,LDZY,UNITQ,KX,GQ,ZY,WRK)
      ELSE IF (NUMINF.EQ.0 .AND. PRBTYP.EQ.'FP') THEN
         CALL F06FBF(N,ZERO,GQ,1)
      ELSE
C
C        Ready for the optimality phase.
C        Set NRZ so that Rz1 is nonsingular.
C
         IF (NRANK.EQ.0) THEN
            IF (LINOBJ) THEN
               CALL DCOPY(N,CQ,1,GQ,1)
            ELSE
               CALL F06FBF(N,ZERO,GQ,1)
            END IF
            NRZ = 0
         ELSE
C
C           Compute GQ = - R' * (transformed residual)
C
            CALL DCOPY(NRANK,RES,1,GQ,1)
            CALL DSCAL(NRANK,(-ONE),GQ,1)
            CALL DTRMV('U','T','N',NRANK,R,LDR,GQ,1)
            IF (NRANK.LT.N) CALL DGEMV('T',NRANK,N-NRANK,-ONE,
     *                                 R(1,NRANK+1),LDR,RES,1,ZERO,
     *                                 GQ(NRANK+1),1)
            IF (LINOBJ) CALL DAXPY(N,ONE,CQ,1,GQ,1)
            UNITGZ = .FALSE.
            ROWNRM = DNRM2(N,R(1,1),LDR)
            IF (ROWNRM.LE.TOLRNK .OR. ABS(R(1,1)).LE.ROWNRM*TOLRNK) THEN
               NRZ = 0
            ELSE
               NRZ = F06KLF(MIN(NRANK,NZ),R,LDR+1,TOLRNK)
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of  E04NCP. (LSGSET)
C
      END
