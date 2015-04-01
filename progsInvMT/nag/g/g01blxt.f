      SUBROUTINE G01BLX(K,N,L,M,JK)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Computes J(K,N,L,M) = Sum from j=max(0,L+M-N) to j=K of
C                            (K!(M-K)!(L-K)!(N-M-L+K)!)/
C                            (J!(M-J)!(L-J)!(N-M-L+J)!)
C
C     .. Parameters ..
      DOUBLE PRECISION  EPS
      PARAMETER         (EPS=1.0D-9)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  JK
      INTEGER           K, L, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, B1, B2, FAC
      INTEGER           I, NSTEP
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     Determine NSTEP, the number of steps in the forward recursion
C
      NSTEP = 0
      FAC = 1.0D0
      A1 = DBLE(K+1)
      B1 = DBLE(N-L-M+K+1)
      A2 = DBLE(L-K)
      B2 = DBLE(M-K)
   20 CONTINUE
      A1 = A1 - 1.0D0
      B1 = B1 - 1.0D0
      A2 = A2 + 1.0D0
      B2 = B2 + 1.0D0
      FAC = FAC*A1*B1/(A2*B2)
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Forward recursion with NSTEP iteration steps
C
      JK = 1.0D0
      DO 40 I = 2, NSTEP
         A1 = A1 + 1.0D0
         B1 = B1 + 1.0D0
         A2 = A2 - 1.0D0
         B2 = B2 - 1.0D0
         JK = 1.0D0 + JK*A1*B1/(A2*B2)
   40 CONTINUE
      RETURN
      END
