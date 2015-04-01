      SUBROUTINE G01BLF(N,L,M,K,PLEK,PGTK,PEQK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Returns the lower tail, upper tail and point probabilities
C     associated with a hypergeometric distribution.
C
C     Let X denote a random variable having a hypergeometric
C     distribution with parameters N, L and M. The routine computes for
C     given N, L, M and K:
C
C     PLEK = Prob (X .LE. K)
C     PGTK = Prob (X .GT. K)
C     PEQK = Prob (X .EQ. K)
C
C     Reference: L. Knusel, Computation of the Chisquare and Poisson
C     distribution, SIAM J. Sci. Stat. Comput. 7, pp 1022-1036, 1986.
C
C     .. Parameters ..
      DOUBLE PRECISION  VARMAX
      PARAMETER         (VARMAX=1.0D6)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01BLF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PEQK, PGTK, PLEK
      INTEGER           IFAIL, K, L, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION  VAR
      INTEGER           IER, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01BLZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IF (N.LE.1) THEN
         VAR = 0.0D0
      ELSE
         VAR = (DBLE(L)*DBLE(M)/DBLE(N-1))*(DBLE(N-L)/DBLE(N))
     *         *(DBLE(N-M)/DBLE(N))
      END IF
C
      IF (N.LT.0) THEN
         IER = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (L.LT.0 .OR. L.GT.N) THEN
         IER = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) L, N
      ELSE IF (M.LT.0 .OR. M.GT.N) THEN
         IER = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) M, N
      ELSE IF (K.LT.0 .OR. K.GT.L .OR. K.GT.M .OR. K.LT.L+M-N) THEN
         IER = 4
         NREC = 2
         WRITE (P01REC,FMT=99996) K, L, M, L + M - N
      ELSE IF (DBLE(N)*X02AJF().GT.1.0D0) THEN
         IER = 5
         NREC = 1
         WRITE (P01REC,FMT=99995) N
      ELSE IF (VAR.GT.VARMAX) THEN
         IER = 6
         NREC = 2
         WRITE (P01REC,FMT=99994) VARMAX, VAR, N, L, M
      ELSE
         IER = 0
         IF (L.EQ.0 .OR. M.EQ.0 .OR. L.EQ.N .OR. M.EQ.N) THEN
            PLEK = 1.0D0
            PGTK = 0.0D0
            PEQK = 1.0D0
         ELSE
            CALL G01BLZ(K,N,L,M,PLEK,PGTK,PEQK)
         END IF
         IFAIL = 0
      END IF
      IF (IER.NE.0) IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.0: N =',I16)
99998 FORMAT (' ** On entry, L.lt.0 or L.gt.N: L =',I16,'  N =',I16)
99997 FORMAT (' ** On entry, M.lt.0 or M.gt.N: M =',I16,'  N =',I16)
99996 FORMAT (' ** On entry, K.lt.0 or K.gt.L or K.gt.M or K.lt.L+M-N:',
     *       '  K =',I16,/'    L =',I16,'  M =',I16,'  L+M-N =',I16)
99995 FORMAT (' ** On entry, N is too large: N =',I16)
99994 FORMAT (' ** On entry, the variance exceeds',1P,D10.1,': varianc',
     *       'e =',D14.6,/'    N =',I16,'  L =',I16,'  M =',I16)
      END
