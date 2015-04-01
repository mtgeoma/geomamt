      SUBROUTINE G13DBX(W,WB,NSM,NS,K,WA,NWA,LW)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 15 REVISED. IER-924 (APR 1991).
C
C        G13DBX UPDATES 1ST K-1 PREDICTION AND BACKWARD
C        PREDICTION COEFFICIENTS AFTER THE KTH COEFFICIENT
C        HAS BEEN CALCULATED.
C
C     .. Scalar Arguments ..
      INTEGER           K, LW, NS, NSM, NWA
C     .. Array Arguments ..
      DOUBLE PRECISION  W(LW), WA(NWA), WB(LW)
C     .. Local Scalars ..
      INTEGER           I, IERR, IFAIL1, K1, K2, K3, K4, NSMQ, NSQ, NSQ1
C     .. External Subroutines ..
      EXTERNAL          F01CTF, G13DBT
C     .. Executable Statements ..
      K1 = K - 1
      IF (K1.LT.1) GO TO 40
      NSMQ = NSM*NSM
      NSQ = NS*NS
      NSQ1 = NSQ + 1
      K4 = K1*NSMQ + 1
      K2 = 1
      K3 = K4 - NSMQ
      DO 20 I = 1, K1
         CALL G13DBT(WA(1),NS,W(K4),NSM,WB(K3),NSM,NS,IERR)
         CALL G13DBT(WA(NSQ1),NS,WB(K4),NSM,W(K2),NSM,NS,IERR)
         IFAIL1 = 1
         CALL F01CTF('N','N',NS,NS,1.0D0,W(K2),NSM,-1.0D0,WA(1),NS,W(K2)
     *               ,NSM,IFAIL1)
         IFAIL1 = 1
         CALL F01CTF('N','N',NS,NS,1.0D0,WB(K3),NSM,-1.0D0,WA(NSQ1),NS,
     *               WB(K3),NSM,IFAIL1)
C
         K2 = K2 + NSMQ
         K3 = K3 - NSMQ
   20 CONTINUE
   40 CONTINUE
      RETURN
      END
