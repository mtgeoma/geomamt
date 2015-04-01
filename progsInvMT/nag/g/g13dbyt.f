      SUBROUTINE G13DBY(C,W,G,NSM,NS,K,WA,NWA,LC,LG)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-797 (DEC 1989).
C     MARK 15 REVISED. IER-925 (APR 1991).
C
C        G13DBY CALCULATES G AT STEP K
C           1 .GE. K .LE. NK
C
C        G =C -W C   - .......... -W   C
C         K  K  1 K-1               K-1 1
C
C     .. Scalar Arguments ..
      INTEGER           K, LC, LG, NS, NSM, NWA
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LC), G(LG), W(LC), WA(NWA)
C     .. Local Scalars ..
      INTEGER           IERR, IFAIL1, K1, K2, NSMQ
C     .. External Subroutines ..
      EXTERNAL          F01CTF, F06QFF, G13DBT
C     .. Executable Statements ..
      NSMQ = NSM*NSM
      K1 = (K-1)*NSMQ + 1
      K2 = 1
      CALL F06QFF('General',NS,NS,C(K1),NSM,G,NS)
      K1 = K1 - NSMQ
      IF (K1.LE.0) GO TO 40
   20 CONTINUE
      CALL G13DBT(WA,NS,W(K2),NSM,C(K1),NSM,NS,IERR)
      IFAIL1 = 1
      CALL F01CTF('N','N',NS,NS,1.0D0,G,NS,-1.0D0,WA,NS,G,NS,IFAIL1)
      K2 = K2 + NSMQ
      K1 = K1 - NSMQ
      IF (K1.GT.0) GO TO 20
   40 CONTINUE
      RETURN
      END
