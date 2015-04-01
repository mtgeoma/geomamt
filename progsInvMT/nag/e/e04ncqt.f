      SUBROUTINE E04NCQ(LINOBJ,SINGLR,UNITGZ,UNITQ,N,NCLIN,NFREE,LDA,
     *                  LDZY,LDR,NRANK,NUMINF,NRZ,KX,CTP,PNORM,A,AP,RES,
     *                  HZ,P,GQ,CQ,R,ZY,WORK)
C     MARK 16 REVISED. IER-1068 (JUL 1993).
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1580 (JUN 1995).
C
C     ******************************************************************
C     E04NCQ  computes the following quantities for  E04NCZ.
C     (1) The vector  (hz1) = (Rz1)(pz1).
C         If X is not yet feasible,  the product is computed directly.
C         If  Rz1 is singular,  hz1  is zero.  Otherwise  hz1  satisfies
C         the equations
C                        Rz1'hz1 = -gz1,
C         where  g  is the total gradient.  If there is no linear term
C         in the objective,  hz1  is set to  dz1  directly.
C     (2) The search direction P (and its 2-norm).  The vector P is
C         defined as  Z*(pz1), where  (pz1)  depends upon whether or
C         not X is feasible and the nonsingularity of  (Rz1).
C         If  NUMINF .GT. 0,  (pz1)  is the steepest-descent direction.
C         Otherwise,  x  is the solution of the  NRZ*NRZ  triangular
C         system   (Rz1)*(pz1) = (hz1).
C     (3) The vector Ap,  where A is the matrix of linear constraints.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level 2 BLAS added 11-June-1986.
C     This version of E04NCQ dated 28-July-1987.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CTP, PNORM
      INTEGER           LDA, LDR, LDZY, N, NCLIN, NFREE, NRANK, NRZ,
     *                  NUMINF
      LOGICAL           LINOBJ, SINGLR, UNITGZ, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AP(*), CQ(*), GQ(N), HZ(*), P(N),
     *                  R(LDR,*), RES(*), WORK(N), ZY(LDZY,*)
      INTEGER           KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  GTP
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, DTRSV, E04NBW, F06FBF
C     .. Executable Statements ..
      IF (SINGLR) THEN
C        ---------------------------------------------------------------
C        The triangular factor for the current objective function is
C        singular,  i.e., the objective is linear along the last column
C        of Z1.  This can only occur when UNITGZ is TRUE.
C        ---------------------------------------------------------------
         IF (NRZ.GT.1) THEN
            CALL DCOPY(NRZ-1,R(1,NRZ),1,P,1)
            CALL DTRSV('U','N','N',NRZ-1,R,LDR,P,1)
         END IF
         P(NRZ) = -ONE
C
         GTP = DDOT(NRZ,GQ,1,P,1)
         IF (GTP.GT.ZERO) CALL DSCAL(NRZ,(-ONE),P,1)
C
         IF (NRZ.LE.NRANK) THEN
            IF (NUMINF.EQ.0) THEN
               IF (UNITGZ) THEN
                  HZ(NRZ) = R(NRZ,NRZ)*P(NRZ)
               ELSE
                  CALL F06FBF(NRZ,ZERO,HZ,1)
               END IF
            ELSE
               HZ(1) = R(1,1)*P(1)
            END IF
         END IF
      ELSE
C        ---------------------------------------------------------------
C        The objective is quadratic in the space spanned by Z1.
C        ---------------------------------------------------------------
         IF (LINOBJ) THEN
            IF (UNITGZ) THEN
               IF (NRZ.GT.1) CALL F06FBF(NRZ-1,ZERO,HZ,1)
               HZ(NRZ) = -GQ(NRZ)/R(NRZ,NRZ)
            ELSE
               CALL DCOPY(NRZ,GQ,1,HZ,1)
               CALL DSCAL(NRZ,(-ONE),HZ,1)
               CALL DTRSV('U','T','N',NRZ,R,LDR,HZ,1)
            END IF
         ELSE
            CALL DCOPY(NRZ,RES,1,HZ,1)
         END IF
C
C        Solve  Rz1*pz1 = hz1.
C
         CALL DCOPY(NRZ,HZ,1,P,1)
         CALL DTRSV('U','N','N',NRZ,R,LDR,P,1)
C
      END IF
C
C     Compute  p = Z1*pz1  and its norm.
C
      IF (LINOBJ) CTP = DDOT(NRZ,CQ,1,P,1)
      PNORM = DNRM2(NRZ,P,1)
C
      CALL E04NBW(1,N,NRZ,NFREE,LDZY,UNITQ,KX,P,ZY,WORK)
C
C     Compute  Ap.
C
      IF (NCLIN.GT.0) THEN
         CALL DGEMV('No transpose',NCLIN,N,ONE,A,LDA,P,1,ZERO,AP,1)
      END IF
C
      RETURN
C
C
C     End of  E04NCQ. (LSGETP)
C
      END
