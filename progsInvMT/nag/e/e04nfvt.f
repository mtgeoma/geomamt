      SUBROUTINE E04NFV(DELREG,POSDEF,STATPT,UNITGZ,UNITQ,N,NCLIN,NFREE,
     *                  LDA,LDQ,LDR,NRZ,ISSAVE,JDSAVE,KX,DNORM,GZDZ,A,
     *                  AD,D,GQ,R,Q,V)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1594 (JUN 1995).
C
C     ******************************************************************
C     E04NFV  computes the following quantities for  E04NFZ.
C
C     (1) The search direction d (and its 2-norm).  The vector d is
C         defined as  Z*(dz), where  (dz)  is defined as follows.
C         If Hz is positive definite, (dz) is the solution of the
C         (NRZ x NRZ)  triangular system  (Rz)'(Rz)*(dz) = - (gz).
C         Otherwise  (dz) is the solution of the triangular system
C         (Rz)(dz) =  gamma ez,  where ez is the NRZ-th unit vector and
C         gamma = -sgn(GQ(NRZ)).
C
C     (2) The vector Ad,  where A is the matrix of linear constraints.
C
C     Original version written 31-December-1986.
C     This version of  E04NFV  dated 21-Dec-1990.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DNORM, GZDZ
      INTEGER           ISSAVE, JDSAVE, LDA, LDQ, LDR, N, NCLIN, NFREE,
     *                  NRZ
      LOGICAL           DELREG, POSDEF, STATPT, UNITGZ, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), AD(*), D(N), GQ(N), Q(LDQ,*),
     *                  R(LDR,*), V(N)
      INTEGER           KX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ATD
      LOGICAL           DELLOW, REVERS
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2
      EXTERNAL          DDOT, DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DSCAL, DTRSV, E04NBW, F06FBF
C     .. Executable Statements ..
C
      IF (POSDEF) THEN
         IF (UNITGZ) THEN
            IF (NRZ.GT.1) CALL F06FBF(NRZ-1,(ZERO),V,1)
            V(NRZ) = -GQ(NRZ)/R(NRZ,NRZ)
         ELSE
            CALL DCOPY(NRZ,GQ,1,V,1)
            CALL DSCAL(NRZ,(-ONE),V,1)
            CALL DTRSV('U','T','N',NRZ,R,LDR,V,1)
         END IF
      ELSE
         IF (NRZ.GT.1) CALL F06FBF(NRZ-1,(ZERO),V,1)
         IF (GQ(NRZ).GT.ZERO) THEN
            V(NRZ) = -ONE
         ELSE
            V(NRZ) = ONE
         END IF
      END IF
C
C     Solve  (Rz)*(dz) =  V.
C
      CALL DCOPY(NRZ,V,1,D,1)
      CALL DTRSV('U','N','N',NRZ,R,LDR,D,1)
C
C     Compute  d = Z*(dz)  and its norm.  Find  gz'dz
C
      DNORM = DNRM2(NRZ,D,1)
      GZDZ = DDOT(NRZ,D,1,GQ,1)
C
      CALL E04NBW(1,N,NRZ,NFREE,LDQ,UNITQ,KX,D,Q,V)
C
C     Compute  Ad.
C
      IF (NCLIN.GT.0) CALL DGEMV('No transpose',NCLIN,N,ONE,A,LDA,D,1,
     *                           ZERO,AD,1)
C
      IF (DELREG .AND. (GZDZ.GT.ZERO .OR. STATPT)) THEN
C        ---------------------------------------------------------------
C        The reduced-gradient norm is small enough that we need to worry
C        about the sign of d.  Make  d  point away from the last deleted
C        constraint.
C        ---------------------------------------------------------------
C        JDSAVE  is the index of the last deleted regular constraint.
C
         IF (JDSAVE.LE.N) THEN
            ATD = D(JDSAVE)
         ELSE
            ATD = AD(JDSAVE-N)
         END IF
C
         DELLOW = ISSAVE .EQ. 1
         IF (DELLOW) THEN
            REVERS = ATD .LT. ZERO
         ELSE
            REVERS = ATD .GT. ZERO
         END IF
C
         IF (REVERS) THEN
            CALL DSCAL(N,(-ONE),D,1)
            IF (NCLIN.GT.0) CALL DSCAL(NCLIN,(-ONE),AD,1)
            GZDZ = -GZDZ
         END IF
      END IF
C
      RETURN
C
C
C     End of  E04NFV.  (QPGETD)
C
      END
