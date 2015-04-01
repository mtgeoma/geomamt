      SUBROUTINE E04NFY(SINGLR,POSDEF,RENEWR,UNITQ,N,NRZ,NFREE,LDQ,LDH,
     *                  LDR,KX,HSIZE,DRZZ,TOLRNK,QPHESS,H,R,Q,HZ,WRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1597 (JUN 1995).
C
C     ******************************************************************
C     E04NFY  is used to compute the last column of the (NRZ x NRZ)
C     triangular factor Rz such that
C                    Hz  =  (Rz)'D(Rz),
C     where Hz is the reduced Hessian Z'HZ, and D is a diagonal
C     matrix.  If  Hz  is positive definite, Rz is the Cholesky factor
C     of  Hz  and  D  is the identity matrix;  otherwise, D(NRZ) is
C     negative or small and the last diagonal of Rz is one.
C
C     The element D(NRZ) is stored in DRZZ.  DRZZ is equal to one if
C     POSDEF is true.
C
C     Original f66 version written by PEG,   March-1982.
C     This version of  E04NFY  dated 27-Jun-1989.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TEN
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,TEN=10.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DRZZ, HSIZE, TOLRNK
      INTEGER           LDH, LDQ, LDR, N, NFREE, NRZ
      LOGICAL           POSDEF, RENEWR, SINGLR, UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LDH,*), HZ(N), Q(LDQ,*), R(LDR,*), WRK(N)
      INTEGER           KX(N)
C     .. Subroutine Arguments ..
      EXTERNAL          QPHESS
C     .. Local Scalars ..
      DOUBLE PRECISION  DRZMAX, DRZMIN, DRZNEW, RDSMIN, RZNORM, RZZ,
     *                  ZTHZ
      INTEGER           J, JTHCOL, K, NRZ1
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2
      EXTERNAL          DNRM2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DTRSV, E04NBW, F06FBF, F06FLF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
C
      IF (NRZ.EQ.0) THEN
         POSDEF = .TRUE.
         RENEWR = .FALSE.
         SINGLR = .FALSE.
         DRZZ = ONE
      ELSE
         IF (RENEWR) THEN
C           ------------------------------------------------------------
C           Compute the first NRZ-1 elements of the last column of Rz
C           and DRZNEW, the square of the last diagonal element.
C           ------------------------------------------------------------
            CALL F06FBF(N,ZERO,WRK,1)
            IF (UNITQ) THEN
C
C              Only bounds are in the working set.  The NRZ-th column of
C              Z is just a column of the identity matrix.
C
               JTHCOL = KX(NRZ)
               WRK(JTHCOL) = ONE
            ELSE
C
C              Expand the new column of  Z  into an n-vector.
C
               DO 20 K = 1, NFREE
                  J = KX(K)
                  WRK(J) = Q(K,NRZ)
   20          CONTINUE
               JTHCOL = 0
            END IF
C
C           Compute the NRZ-th column of Z'HZ.
C
            CALL QPHESS(N,JTHCOL,H,LDH,WRK,HZ)
            CALL E04NBW(4,N,NRZ,NFREE,LDQ,UNITQ,KX,HZ,Q,WRK)
            CALL DCOPY(NRZ,HZ,1,R(1,NRZ),1)
C
            NRZ1 = NRZ - 1
            ZTHZ = R(NRZ,NRZ)
            DRZNEW = ZTHZ
            IF (NRZ1.GT.0) THEN
               CALL DTRSV('U','T','N',NRZ1,R,LDR,R(1,NRZ),1)
               RZNORM = DNRM2(NRZ1,R(1,NRZ),1)
               DRZNEW = ZTHZ - RZNORM*RZNORM
            END IF
C
            R(NRZ,NRZ) = ONE
            DRZZ = DRZNEW
C
C           Update the estimate of the norm of the Hessian.
C
            HSIZE = MAX(HSIZE,ABS(ZTHZ))
C
         END IF
C
         DRZNEW = DRZZ*R(NRZ,NRZ)**2
C
C        ---------------------------------------------------------------
C        Attempt to compute RZZ, the square root of  DRZNEW.  The last
C        diagonal of Rz.  The variables POSDEF and SINGLR are set here.
C        They are used to indicate if the new Z'HZ is positive definite
C        or singular.  If the required diagonal modification is large
C        the last row and column of Rz are marked for recomputation next
C        iteration.
C        ---------------------------------------------------------------
C        RDSMIN is the square of the smallest allowable diagonal element
C        for a positive-definite Cholesky factor.  Note that the test
C        for positive definiteness is unavoidably scale dependent.
C
         IF (NRZ.EQ.1) THEN
            RDSMIN = TOLRNK*HSIZE
         ELSE
            CALL F06FLF(NRZ,R,LDR+1,DRZMAX,DRZMIN)
            RDSMIN = (TOLRNK*DRZMAX)*DRZMAX
         END IF
C
         POSDEF = DRZNEW .GT. RDSMIN
C
         IF (POSDEF) THEN
            DRZZ = ONE
            RZZ = SQRT(DRZNEW)
            RENEWR = .FALSE.
            SINGLR = .FALSE.
         ELSE
            DRZZ = DRZNEW
            RZZ = ONE
            SINGLR = DRZNEW .GE. -RDSMIN
            RENEWR = DRZNEW .LT. -TEN*HSIZE
         END IF
C
         R(NRZ,NRZ) = RZZ
C
      END IF
C
      RETURN
C
C
C     End of  E04NFY.  (QPCOLR)
C
      END
