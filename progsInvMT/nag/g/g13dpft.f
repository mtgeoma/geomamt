      SUBROUTINE G13DPF(K,N,Z,IK,M,MAXLAG,PARLAG,SE,QQ,X,PVALUE,LOGLHD,
     *                  WORK,LWORK,IWORK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DPF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, K, LWORK, M, MAXLAG, N
C     .. Array Arguments ..
      DOUBLE PRECISION  LOGLHD(M), PARLAG(IK,IK,M), PVALUE(M),
     *                  QQ(IK,IK,M), SE(IK,IK,M), WORK(LWORK), X(M),
     *                  Z(IK,N)
      INTEGER           IWORK(K*M)
C     .. Local Scalars ..
      INTEGER           IERROR, J, L, LW1, LW2, LW3, LW4, LW5, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13DPZ
C     .. Executable Statements ..
C
C     ******************************************************************
C
C     This subroutine calculates the sample partial autoregression
C     matrices of a multivariate time series using the algorithm
C     described by Wei (see also Box and Tiao (1981))
C
C     ******************************************************************
C
C     First test for errors in the input arguments
C
      MAXLAG = 0
      IERROR = 1
      NREC = 1
      L = M*K + 1
      IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99999) K
      ELSE IF (N.LT.4) THEN
         WRITE (P01REC,FMT=99998) N
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99997) IK, K
      ELSE IF (M.LT.1) THEN
         WRITE (P01REC,FMT=99996) M
      ELSE IF (N-M-L.LT.K) THEN
         NREC = 2
         WRITE (P01REC,FMT=99995) K, M, N
      ELSE
         IERROR = 0
C
C        Test whether the workspace array is big enough
C
         J = (K+1)*K + L*(4+K) + 2*L*L
         IF (LWORK.LT.J) THEN
            IERROR = 1
            NREC = 2
            WRITE (P01REC,FMT=99994) LWORK, J
            GO TO 20
         END IF
         NREC = 0
C
C        Partition workspace array
C
         LW1 = 1
         LW2 = LW1 + L
         LW3 = LW2 + L*(L+K)
         LW4 = LW3 + L*(L+1)/2
         LW5 = LW4 + K*K
C
         CALL G13DPZ(K,N,Z,IK,M,PARLAG,SE,QQ,X,PVALUE,LOGLHD,MAXLAG,L,
     *               WORK(LW1),WORK(LW2),WORK(LW3),WORK(LW4),WORK(LW5),
     *               IWORK,IERROR)
         IF (IERROR.NE.0) THEN
            NREC = 3
            WRITE (P01REC,FMT=99993) MAXLAG
         END IF
C
      END IF
   20 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, K.lt.1: K = ',I16)
99998 FORMAT (' ** On entry, N.lt.4: N = ',I16)
99997 FORMAT (' ** On entry, IK.lt.K: IK = ',I16,' and K = ',I16)
99996 FORMAT (' ** On entry, M.lt.1: M = ',I16)
99995 FORMAT (' ** On entry, N-M-(K*M+1).lt.K: K = ',I16,' M = ',I16,
     *       /'    and N = ',I16)
99994 FORMAT (' ** On entry, LWORK is too small: LWORK = ',I16,/'    I',
     *       't should be at least',I16)
99993 FORMAT (' ** The recursive equations used to compute the partial',
     *       ' autoregression',/'    matrices are ill-conditioned. The',
     *       'y have been computed up to',/'    lag ',I16)
      END
