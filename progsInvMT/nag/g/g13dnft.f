      SUBROUTINE G13DNF(K,N,M,IK,R0,R,MAXLAG,PARLAG,X,PVALUE,WORK,LWORK,
     *                  IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     This subroutine calculates the sample partial lag correlation
C     matrices of a multivariate time series using the recursive
C     algorithm of Heyse and Wei.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DNF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, K, LWORK, M, MAXLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION  PARLAG(IK,IK,M), PVALUE(M), R(IK,IK,M),
     *                  R0(IK,K), WORK(LWORK), X(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DF, SUM
      INTEGER           I, IERROR, IFAIL2, J, K2, L, L2, LW1, LW2, LW3,
     *                  LW4, LW5, LW6, LW7, LW8, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13DNZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     First test for errors in the input arguments
C
      IERROR = 1
      NREC = 1
      IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99999) K
      ELSE IF (N.LT.2) THEN
         WRITE (P01REC,FMT=99998) N
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99997) IK, K
      ELSE IF (M.LT.1) THEN
         WRITE (P01REC,FMT=99996) M
      ELSE IF (M.GE.N) THEN
         WRITE (P01REC,FMT=99995) M, N
      ELSE
         IERROR = 0
         NREC = 0
         MAXLAG = 0
C
C        Test whether the workspace array is big enough
C
         K2 = K*K
         J = (5*M+6)*K2 + K
         IF (LWORK.LT.J) THEN
            IERROR = 1
            NREC = 2
            WRITE (P01REC,FMT=99994) LWORK, J
            GO TO 160
         END IF
C
C        Recover the sample cross-covariance matrices and store in the
C        workspace array
C
         L2 = 0
         DO 40 I = 1, K
            DO 20 J = 1, K
               L2 = L2 + 1
               IF (I.EQ.J) THEN
                  WORK(L2) = R0(I,I)*R0(I,I)
               ELSE
                  WORK(L2) = R0(I,J)*R0(I,I)*R0(J,J)
               END IF
   20       CONTINUE
   40    CONTINUE
         DO 100 L = 1, M
            DO 80 I = 1, K
               DO 60 J = 1, K
                  L2 = L2 + 1
                  WORK(L2) = R(I,J,L)*R0(I,I)*R0(J,J)
   60          CONTINUE
   80       CONTINUE
  100    CONTINUE
C
C        Set PARLAG(1) = rho(1) and calculate X(1) and PVALUE(1)
C
         SUM = 0.0D0
         DO 140 J = 1, K
            DO 120 I = 1, K
               PARLAG(I,J,1) = R(I,J,1)
               SUM = SUM + PARLAG(I,J,1)*PARLAG(I,J,1)
  120       CONTINUE
  140    CONTINUE
         X(1) = DBLE(N)*SUM
         DF = DBLE(K2)
         IFAIL2 = 0
C        Note that G01ECF cannot fail
         PVALUE(1) = G01ECF('Upper',X(1),DF,IFAIL2)
         MAXLAG = 1
C        Partition workspace array
         LW1 = 1
         LW2 = LW1 + (M+1)*K2
         LW3 = LW2 + K2
         LW4 = LW3 + K2
         LW5 = LW4 + K2
         LW6 = LW5 + (K+K)*M*K
         LW7 = LW6 + (K+K)*M*K
         LW8 = LW7 + K2
         CALL G13DNZ(K,N,M,IK,PARLAG,WORK(LW1),WORK(LW2),WORK(LW3),
     *               WORK(LW4),WORK(LW5),WORK(LW6),WORK(LW7),WORK(LW8),
     *               X,PVALUE,MAXLAG,IERROR)
         IF (IERROR.NE.0) THEN
            IERROR = 2
            NREC = 2
            WRITE (P01REC,FMT=99993) MAXLAG
         END IF
C
      END IF
  160 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT ('  ** On entry, K.lt.1 : K = ',I16)
99998 FORMAT ('  ** On entry, N.lt.2 : N = ',I16)
99997 FORMAT ('  ** On entry, IK.lt.K : IK = ',I16,'  and K = ',I16)
99996 FORMAT ('  ** On entry, M.lt.1 : M = ',I16)
99995 FORMAT ('  ** On entry, M.ge.N : M = ',I16,'  and N = ',I16)
99994 FORMAT ('  ** On entry, LWORK is too small : LWORK = ',I16,' and',
     *       ' must be at',/'     least ',I16)
99993 FORMAT ('  ** The recursive equations used to compute the partia',
     *       'l lag correlation',/'     matrices are ill-conditioned (',
     *       'they have been computed up to lag ',I4,')')
      END
