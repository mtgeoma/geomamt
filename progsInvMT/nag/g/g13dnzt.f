      SUBROUTINE G13DNZ(K,N,M,IK,PARLAG,GAMMA,VU,VV,VVU,ALPHAT,BETAT,B,
     *                  WK,X,PVALUE,MAXLAG,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1692 (JUN 1995).
C
C     For details of this algorithm see
C       WEI, W.W.S., Time Series Analysis: Univariate and multivariate
C                    methods; Chpater 14. (In particular page 361).
C     Note that in this subroutine only the transposes of the alphas
C     and betas are stored. They arise naturally from solving the
C     system of equations with F04ABF and all the recursive updating
C     may be written in terms of the transposes.
C     Note this means that we are storing Vu and Vv transpose but they
C     are symmetric so this makes no difference.
C     GAMMA stores the Gamma(0), Gamma(1), etc in sequential order.
C     Each Gamma(i) is stored row by row, that is GAMMA(2) contains
C     the (1,2) element of Gamma(0).
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, K, M, MAXLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ALPHAT(K+K,M*K), B(K,K), BETAT(K+K,M*K),
     *                  GAMMA((M+1)*K*K), PARLAG(IK,IK,M), PVALUE(M),
     *                  VU(K,K), VV(K,K), VVU(K,K), WK(K*(K+1)), X(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DF, SUM
      INTEGER           I, IFAIL2, IS, IS1, J, K1, K2, KK, S
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      EXTERNAL          G01ECF
C     .. External Subroutines ..
      EXTERNAL          F04ABF, DCOPY, DGEMV, F06QFF, DGEMM
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE, SQRT
C     .. Executable Statements ..
C
C     Calculate ALPHAT(1,1) and BETAT(1,1)
C
C     First copy GAMMA(0) onto upper triangle of VU
C
      K2 = K*K
      KK = K + K
      DO 20 I = 1, K
         CALL DCOPY(I,GAMMA((I-1)*K+1),1,VU(1,I),1)
   20 CONTINUE
C
C     Set B = GAMMA(1)
C
      DO 40 I = 1, K
         CALL DCOPY(K,GAMMA(K2+(I-1)*K+1),1,B(I,1),K)
   40 CONTINUE
C
      IFAIL = 1
      CALL F04ABF(VU,K,B,K,K,K,ALPHAT,KK,WK,WK(K+1),K,IFAIL)
      IF (IFAIL.NE.0) RETURN
C
C     Set B = GAMMA(1)'
C
      CALL F06QFF('G',K,K,GAMMA(K2+1),K,B(1,1),K)
C
      IFAIL = 1
      CALL F04ABF(VU,K,B,K,K,K,BETAT,KK,WK,WK(K+1),K,IFAIL)
      IF (IFAIL.NE.0) RETURN
C
C     Perform recursions
C
      DF = DBLE(K2)
      DO 220 S = 2, M
C
C        Update upper triangles of Vu(s) and Vv(s)
C
         IS = MOD(S,2)
         IS1 = MOD(S-1,2)
         CALL F06QFF('U',K,K,GAMMA,K,VU,K)
         CALL F06QFF('U',K,K,GAMMA,K,VV,K)
         DO 80 K1 = 1, S - 1
            DO 60 J = 1, K
               CALL DGEMV('N',J,K,-1.0D0,GAMMA(K1*K2+1),K,
     *                     ALPHAT(IS*K+1,(K1-1)*K+J),1,1.0D0,VU(1,J),1)
               CALL DGEMV('T',K,J,-1.0D0,GAMMA(K1*K2+1),K,
     *                     BETAT(IS*K+1,(K1-1)*K+J),1,1.0D0,VV(1,J),1)
   60       CONTINUE
   80    CONTINUE
C
C        Update Vvu(s)
C
         DO 100 I = 1, K
            CALL DCOPY(K,GAMMA(S*K2+I),K,VVU(1,I),1)
  100    CONTINUE
         DO 120 K1 = 1, S - 1
            CALL DGEMM('T','N',K,K,K,-1.0D0,GAMMA((S-K1)*K2+1),K,
     *                  ALPHAT(IS*K+1,(K1-1)*K+1),KK,1.0D0,VVU,K)
  120    CONTINUE
C
C        Calculate PARLAG(s) and X(s) and PVALUE(s)
C
         SUM = 0.0D0
         DO 160 I = 1, K
            DO 140 J = 1, K
               IF (VV(I,I).GT.0.0D0 .AND. VU(J,J).GT.0.0D0) THEN
                  PARLAG(I,J,S) = VVU(I,J)/SQRT(VV(I,I)*VU(J,J))
                  SUM = SUM + PARLAG(I,J,S)*PARLAG(I,J,S)
               ELSE
                  IFAIL = 1
                  RETURN
               END IF
  140       CONTINUE
  160    CONTINUE
         X(S) = DBLE(N)*SUM
         IFAIL2 = 0
C
C        Note that G01ECF cannot fail
C
         PVALUE(S) = G01ECF('Upper',X(S),DF,IFAIL2)
         MAXLAG = S
C
         IF (S.LT.M) THEN
C
C           Update ALPHAT(s,s)
C
            IFAIL = 1
            CALL F04ABF(VV,K,VVU,K,K,K,ALPHAT(IS1*K+1,(S-1)*K+1),KK,WK,
     *                  WK(K+1),K,IFAIL)
            IF (IFAIL.NE.0) RETURN
C
C           Update BETAT(s,s)
C
C           Set B = Vvu(s)'
C
            DO 180 I = 1, K
               CALL DCOPY(K,VVU(1,I),1,B(I,1),K)
  180       CONTINUE
            IFAIL = 1
            CALL F04ABF(VU,K,B,K,K,K,BETAT(IS1*K+1,(S-1)*K+1),KK,WK,
     *                  WK(K+1),K,IFAIL)
            IF (IFAIL.NE.0) RETURN
C
C           Calculate ALPHAT(s,k) for k = 1,2,...,s-1
C           and BETAT(s,k) for k = 1,2,...,s-1
C
            DO 200 K1 = 1, S - 1
C
               CALL F06QFF('G',K,K,ALPHAT(IS*K+1,(K1-1)*K+1),KK,
     *                     ALPHAT(IS1*K+1,(K1-1)*K+1),KK)
C
               CALL F06QFF('G',K,K,BETAT(IS*K+1,(K1-1)*K+1),KK,
     *                     BETAT(IS1*K+1,(K1-1)*K+1),KK)
C
               CALL DGEMM('N','N',K,K,K,-1.0D0,BETAT(IS*K+1,(S-K1-1)
     *                     *K+1),KK,ALPHAT(IS1*K+1,(S-1)*K+1),KK,1.0D0,
     *                     ALPHAT(IS1*K+1,(K1-1)*K+1),KK)
C
               CALL DGEMM('N','N',K,K,K,-1.0D0,ALPHAT(IS*K+1,(S-K1-1)
     *                     *K+1),KK,BETAT(IS1*K+1,(S-1)*K+1),KK,1.0D0,
     *                     BETAT(IS1*K+1,(K1-1)*K+1),KK)
  200       CONTINUE
C
         END IF
C
  220 CONTINUE
C
      RETURN
      END
