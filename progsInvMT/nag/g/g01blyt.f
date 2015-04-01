      SUBROUTINE G01BLY(K,N,L,M,IK)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes I(K,N,L,M) = Sum from j=K+1 to j=min(M,N) of
C                            (K!(M-K)!(L-K)!(N-M-L+K)!)/
C                            (J!(M-J)!(L-J)!(N-M-L+J)!)
C
C     .. Parameters ..
      DOUBLE PRECISION  EPS
      PARAMETER         (EPS=1.0D-9)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  IK
      INTEGER           K, L, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, B1, B2, FAC
      INTEGER           I, NSTEP
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     Determine NSTEP, the number of steps in the backward recursion
C
      NSTEP = 0
      FAC = 1.0D0
      A1 = DBLE(L-K)
      B1 = DBLE(M-K)
      A2 = DBLE(K+1)
      B2 = DBLE(N-L-M+K+1)
   20 CONTINUE
      FAC = FAC*A1*B1/(A2*B2)
      A1 = A1 - 1.0D0
      B1 = B1 - 1.0D0
      A2 = A2 + 1.0D0
      B2 = B2 + 1.0D0
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Backward recursion with NSTEP iteration steps
C
      IK = 0.0D0
      DO 40 I = 1, NSTEP
         A1 = A1 + 1.0D0
         B1 = B1 + 1.0D0
         A2 = A2 - 1.0D0
         B2 = B2 - 1.0D0
         IK = (1.0D0+IK)*A1*B1/(A2*B2)
   40 CONTINUE
      RETURN
      END
