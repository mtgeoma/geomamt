      SUBROUTINE G13DJX(Z,K,N,LMAX,V,IK,ND,MEAN,PAR,NPAR,PHISTA,IP,IQ,
     *                  IPSTAR,PREDZ,B)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 15B REVISED. IER-958 (NOV 1991).
C     MARK 17 REVISED. IER-1689 (JUN 1995).
C     .. Scalar Arguments ..
      INTEGER           IK, IP, IPSTAR, IQ, K, LMAX, N, ND, NPAR
      CHARACTER         MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  B(K), PAR(NPAR), PHISTA(K,*),
     *                  PREDZ(IK,LMAX), V(IK,ND), Z(K,N+LMAX)
C     .. Local Scalars ..
      INTEGER           I, J, K3, L, L2, M2
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06FBF, DGEMV
C     .. Executable Statements ..
C
      K3 = K*K
      M2 = IP*K*K
C
      IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
C
C        set THETA(0) = PHI(B) * MU
C
         L2 = (IP+IQ)*K3
         DO 15 I = 1, K
            B(I) = PAR(L2+I)
            DO 10 J = 1, K
               DO 5 L = 1, IP
                  B(I) = B(I) - PAR((L-1)*K3+(I-1)*K+J)*PAR(L2+J)
    5          CONTINUE
   10       CONTINUE
   15    CONTINUE
C
      ELSE
         CALL F06FBF(K,0.0D0,B,1)
      END IF
C
C     predict Z(n+1), ..., Z(n+lmax)
C
      DO 80 L = 1, LMAX
         CALL DCOPY(K,B,1,Z(1,N+L),1)
C
         DO 40 I = 1, IPSTAR
            CALL DGEMV('N',K,K,1.0D0,PHISTA(1,(I-1)*K+1),K,Z(1,N+L-I),
     *                  1,1.0D0,Z(1,N+L),1)
   40    CONTINUE
C
         DO 60 I = L, IQ
            CALL DGEMV('T',K,K,-1.0D0,PAR(M2+(I-1)*K3+1),K,V(1,ND+L-I),
     *                  1,1.0D0,Z(1,N+L),1)
   60    CONTINUE
C
         CALL DCOPY(K,Z(1,N+L),1,PREDZ(1,L),1)
C
   80 CONTINUE
C
      RETURN
      END
