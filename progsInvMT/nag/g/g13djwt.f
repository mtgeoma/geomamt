      SUBROUTINE G13DJW(PHISTA,K,IPSTAR,IP,IQ,PAR,NPAR,SEFZ,IK,LMAX,QQ,
     *                  PSI,A,C,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1688 (JUN 1995).
C     .. Scalar Arguments ..
      INTEGER           IFAULT, IK, IP, IPSTAR, IQ, K, LMAX, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(K,K), C(K,K), PAR(NPAR), PHISTA(K,*),
     *                  PSI(K,*), QQ(IK,K), SEFZ(IK,LMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I2, J, J2, K3, L, L2, M2
C     .. External Subroutines ..
      EXTERNAL          F06ECF, DCOPY, F06FBF, F06QFF, DGEMM
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      IFAULT = 0
      K3 = K*K
      M2 = IP*K*K
      CALL F06QFF('G',K,K,QQ,IK,A,K)
C
      DO 20 I = 1, K
         IF (QQ(I,I).GT.0.0D0) THEN
            SEFZ(I,1) = QQ(I,I)
         ELSE
            IFAULT = 2
            RETURN
         END IF
   20 CONTINUE
C
      DO 220 L = 1, LMAX
C
         IF (L.LE.LMAX-1) THEN
C
C           calculate psi(l)
C
            IF (L.LE.IPSTAR) THEN
               DO 40 J = 1, K
                  CALL DCOPY(K,PHISTA(1,(L-1)*K+J),1,PSI(1,(L-1)*K+J),
     *                        1)
   40          CONTINUE
            ELSE
               DO 60 J = 1, K
                  CALL F06FBF(K,0.0D0,PSI(1,(L-1)*K+J),1)
   60          CONTINUE
            END IF
C
            DO 80 L2 = 1, MIN(L-1,IPSTAR)
               CALL DGEMM('N','N',K,K,K,1.0D0,PHISTA(1,(L2-1)*K+1),K,
     *                     PSI(1,(L-L2-1)*K+1),K,1.0D0,PSI(1,(L-1)*K+1),
     *                     K)
   80       CONTINUE
C
            IF (L.LE.IQ) THEN
               DO 100 J = 1, K
                  CALL F06ECF(K,-1.0D0,PAR(M2+(L-1)*K3+J),K,PSI(1,(L-1)
     *                        *K+J),1)
  100          CONTINUE
            END IF
C
         END IF
C
         IF (L.GE.2) THEN
C
C           compute V(l)
C
            DO 180 I = 1, K
               DO 160 J = 1, I
                  SUM = A(I,J)
                  DO 140 I2 = 1, K
                     DO 120 J2 = 1, K
                        SUM = SUM + PSI(I,(L-2)*K+I2)*QQ(I2,J2)*PSI(J,
     *                        (L-2)*K+J2)
  120                CONTINUE
  140             CONTINUE
                  C(I,J) = SUM
  160          CONTINUE
  180       CONTINUE
C
            DO 200 I = 1, K
               IF (C(I,I).GT.0.0D0) THEN
                  SEFZ(I,L) = C(I,I)
               ELSE
                  IFAULT = 2
                  RETURN
               END IF
  200       CONTINUE
C
C           copy C onto A
C
            CALL F06QFF('L',K,K,C,K,A,K)
         END IF
C
  220 CONTINUE
C
      RETURN
      END
