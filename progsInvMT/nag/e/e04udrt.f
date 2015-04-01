      SUBROUTINE E04UDR(UNITQ,N,NFREE,NZ,NQ,NROWR,IPERM,KX,GQ,R,ZY,WORK,
     *                  QRWORK)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-793 (DEC 1989).
C     MARK 16 REVISED. IER-1094 (JUL 1993).
C     MARK 17 REVISED. IER-1614 (JUN 1995).
C
C     ******************************************************************
C     E04UDR  bounds the condition estimator of the transformed Hessian.
C     On exit, R is of the form
C                  ( DRz   0     )
C                  (  0  sigma*I )
C     where D is a diagonal matrix such that DRz has a bounded condition
C     number,  I is the identity matrix and sigma  is the geometric mean
C     of the largest and smallest elements of DRz. The QR factorization
C     with interchanges is used to give diagonals of DRz that are
C     decreasing in modulus.
C
C     Systems Optimization Laboratory, Stanford University.
C
C     Original version of E04UDR dated  4-August-1986.
C     This version dated  14-Sep-92.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF, ONE
      PARAMETER         (ZERO=0.0D+0,HALF=0.5D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           N, NFREE, NQ, NROWR, NZ
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  GQ(N), QRWORK(2*N), R(NROWR,*), WORK(N),
     *                  ZY(NQ,*)
      INTEGER           IPERM(N), KX(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DRMAX, DRMIN, RCNDBD, RFROBN
C     .. Local Scalars ..
      DOUBLE PRECISION  DRGM, DRGS, GJMAX, SCLE, SUMSQ
      INTEGER           INFO, J, JMAX, JSAVE, NRANK
C     .. External Functions ..
      DOUBLE PRECISION  F06BMF
      INTEGER           F06KLF
      EXTERNAL          F06BMF, F06KLF
C     .. External Subroutines ..
      EXTERNAL          DSWAP, F01QFF, F06FBF, F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Common blocks ..
      COMMON            /EE04NB/RCNDBD, RFROBN, DRMAX, DRMIN
C     .. Executable Statements ..
C
C     ==================================================================
C     Bound the condition estimator of Q'HQ.
C     The scheme used here reduces the modulus of the larger
C     diagonals while increasing the modulus of the smaller one.
C     ==================================================================
      IF (NZ.GT.1) THEN
C        ---------------------------------------------------------------
C        Refactorize Rz.  Interchanges are used to give diagonals
C        of decreasing magnitude.
C        ---------------------------------------------------------------
         DO 20 J = 1, NZ - 1
            CALL F06FBF(NZ-J,ZERO,R(J+1,J),1)
   20    CONTINUE
C
         CALL F01QFF('Column iterchanges',NZ,NZ,R,NROWR,WORK,IPERM,
     *               QRWORK,INFO)
C
         DO 40 J = 1, NZ
            JMAX = IPERM(J)
            IF (JMAX.GT.J) THEN
               IF (UNITQ) THEN
                  JSAVE = KX(JMAX)
                  KX(JMAX) = KX(J)
                  KX(J) = JSAVE
               ELSE
                  CALL DSWAP(NFREE,ZY(1,JMAX),1,ZY(1,J),1)
               END IF
C
               GJMAX = GQ(JMAX)
               GQ(JMAX) = GQ(J)
               GQ(J) = GJMAX
            END IF
   40    CONTINUE
      END IF
C
      DRGM = ONE
C
      IF (NZ.GT.0) THEN
         NRANK = F06KLF(NZ,R,NROWR+1,ONE/RCNDBD)
         DRGM = HALF*SQRT(ABS(R(1,1)*R(NRANK,NRANK)))
         DRGS = ABS(R(1,1))/RCNDBD
C
         IF (NZ.GT.NRANK) THEN
            DO 60 J = NRANK + 1, NZ
               CALL F06FBF(J-1,ZERO,R(1,J),1)
   60       CONTINUE
            CALL F06FBF(NZ-NRANK,DRGS,R(NRANK+1,NRANK+1),NROWR+1)
         END IF
      END IF
C
C     ------------------------------------------------------------------
C     Reset the range-space partition of the Hessian.
C     ------------------------------------------------------------------
      IF (NZ.LT.N) THEN
         DO 80 J = NZ + 1, N
            CALL F06FBF(J,ZERO,R(1,J),1)
   80    CONTINUE
         CALL F06FBF(N-NZ,DRGM,R(NZ+1,NZ+1),NROWR+1)
      END IF
C
C     Recompute the Frobenius norm of R.
C
      SCLE = SQRT(DBLE(N-NZ))*DRGM
      SUMSQ = ONE
      DO 100 J = 1, NZ
         CALL F06FJF(J,R(1,J),1,SCLE,SUMSQ)
  100 CONTINUE
      RFROBN = F06BMF(SCLE,SUMSQ)
C
      RETURN
C
C     End of  E04UDR. (NPRSET)
C
      END
