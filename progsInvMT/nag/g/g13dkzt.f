      SUBROUTINE G13DKZ(K,LMAX,ZOLD,PREDZ,SEFZ,IK,M,PSI,Z,MS,TR,ZT,
     *                  MLAST,V,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           IFAULT, IK, K, LMAX, M, MLAST
C     .. Array Arguments ..
      DOUBLE PRECISION  MS(K,LMAX), PREDZ(IK,LMAX), PSI(K,*),
     *                  SEFZ(IK,LMAX), TR(K), V(IK,M), Z(K,LMAX),
     *                  ZOLD(IK,M), ZT(K,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA
      INTEGER           I, L2, LLL, T
      CHARACTER*1       TRAN
C     .. External Subroutines ..
      EXTERNAL          F06PAF, F06QFF, F06QHF, G13DJV, G13DLZ
C     .. Executable Statements ..
C
C     First transform Z(n+1), Z(n+2), ..., Z(n+l) and store in the array
C     ZT
C
      DO 20 I = 1, K
C
C        Test that transformations in REF not corrupted.
C
         AA = TR(I)
         IF (AA.EQ.100.0D0) THEN
            TRAN = 'N'
         ELSE IF (AA.EQ.200.0D0) THEN
            TRAN = 'L'
         ELSE IF (AA.EQ.300.0D0) THEN
            TRAN = 'S'
         ELSE
            IFAULT = 2
            RETURN
         END IF
         CALL G13DLZ(I,K,M,ZOLD,IK,TRAN,ZT,IFAULT)
         IF (IFAULT.EQ.1) THEN
            IFAULT = 3
            RETURN
         END IF
   20 CONTINUE
C
C     update predictors
C
C     Z array contains transformed predictors of Z(n+1), Z(n+2), ...,
C     Z(n+lmax)
C
      DO 80 T = 1, M
C
C        calculate a(n+t)
C
         DO 40 I = 1, K
            V(I,T) = ZT(I,T) - Z(I,MLAST+T)
   40    CONTINUE
C
C        update predictors given additional observation Z(n+t)
C
         DO 60 L2 = 1, LMAX - MLAST - T
            CALL F06PAF('N',K,K,1.0D0,PSI(1,(L2-1)*K+1),K,V(1,T),1,
     *                  1.0D0,Z(1,MLAST+T+L2),1)
   60    CONTINUE
C
   80 CONTINUE
C
C     update PREDZ and SEFZ
C
      CALL F06QHF('G',K,M,0.0D0,0.0D0,SEFZ(1,MLAST+1),IK)
      CALL F06QFF('G',K,M,ZOLD,IK,PREDZ(1,MLAST+1),IK)
      CALL F06QFF('G',K,LMAX-MLAST-M,MS,K,SEFZ(1,MLAST+M+1),IK)
      CALL F06QFF('G',K,LMAX-MLAST-M,Z(1,MLAST+M+1),K,PREDZ(1,MLAST+M+1)
     *            ,IK)
C
C     back transform the predictors and their standard errors
C
      LLL = MLAST + M + 1
      DO 100 I = 1, K
         AA = TR(I)
         IF (AA.EQ.100.0D0) THEN
            TRAN = 'N'
         ELSE IF (AA.EQ.200.0D0) THEN
            TRAN = 'L'
         ELSE IF (AA.EQ.300.0D0) THEN
            TRAN = 'S'
         END IF
         CALL G13DJV(I,LMAX,IK,TRAN,PREDZ,SEFZ,LLL,IFAULT)
         IF (IFAULT.EQ.1) IFAULT = 4
         IF (IFAULT.NE.0) RETURN
  100 CONTINUE
C
      RETURN
      END
