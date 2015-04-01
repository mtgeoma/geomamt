      SUBROUTINE G13DJZ(K,N,IP,IQ,MEAN,PAR,NPAR,QQ,IK,ZOLD,TR,ID,DELTA,
     *                  IDMAX,V,LMAX,ND,PREDZ,SEFZ,Z,A,C,PHISTA,PSI,REF,
     *                  LR,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 15B REVISED. IER-959 (NOV 1991).
C     .. Scalar Arguments ..
      INTEGER           IDMAX, IFAULT, IK, IP, IQ, K, LMAX, LR, N, ND,
     *                  NPAR
      CHARACTER         MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(K,K), C(K,K), DELTA(IK,*), PAR(NPAR),
     *                  PHISTA(K,*), PREDZ(IK,LMAX), PSI(K,*), QQ(IK,K),
     *                  REF(LR), SEFZ(IK,LMAX), V(IK,ND), Z(K,N+LMAX),
     *                  ZOLD(IK,N)
      INTEGER           ID(K)
      CHARACTER         TR(K)
C     .. Local Scalars ..
      INTEGER           I, IPSTAR, LLL, MK
C     .. External Subroutines ..
      EXTERNAL          F06QFF, G13DJV, G13DJW, G13DJX, G13DJY, G13DLZ
C     .. Executable Statements ..
C
C     transform each series (if necessary)
C
      DO 20 I = 1, K
         CALL G13DLZ(I,K,N,ZOLD,IK,TR(I),Z,IFAULT)
         IF (IFAULT.EQ.1) RETURN
   20 CONTINUE
C
C     set   PHISTA(B) = PHI(B) * DELTA(B)
C
      IPSTAR = IP + IDMAX
      IF (IPSTAR.GT.0) CALL G13DJY(K,IP,PAR,NPAR,ID,DELTA,IDMAX,IK,
     *                             PHISTA,IPSTAR)
C
C     calculate the forecasts for l = 1, 2, ..., lmax
C
      CALL G13DJX(Z,K,N,LMAX,V,IK,ND,MEAN,PAR,NPAR,PHISTA,IP,IQ,IPSTAR,
     *            PREDZ,C)
C
C     copy predictors onto REF array
C
      CALL F06QFF('G',K,LMAX,PREDZ,IK,REF,K)
C
C     calculate the psi weights up to lag LMAX-1 and at the same time
C     the standard deviations of the forecast errors
C
      CALL G13DJW(PHISTA,K,IPSTAR,IP,IQ,PAR,NPAR,SEFZ,IK,LMAX,QQ,PSI,A,
     *            C,IFAULT)
      IF (IFAULT.EQ.2) RETURN
C
C     copy MSE's onto REF array
C
      MK = LMAX*K
      CALL F06QFF('G',K,LMAX,SEFZ,IK,REF(MK+1),K)
C
C     apply corrections to forecasts and their MSE's
C
      LLL = 1
      DO 40 I = 1, K
         CALL G13DJV(I,LMAX,IK,TR(I),PREDZ,SEFZ,LLL,IFAULT)
         IF (IFAULT.EQ.1) THEN
            IFAULT = 3
            RETURN
         END IF
C
C        Note that IFAULT = 2 in G13DJV cannot occur here.
C
   40 CONTINUE
C
      RETURN
      END
