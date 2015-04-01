      SUBROUTINE G01BJF(N,P,K,PLEK,PGTK,PEQK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Returns the lower tail, upper tail and point probabilities
C     associated with a Binomial distribution.
C
C     Let X denote a random variable having a binomial distribution with
C     parameters N and P. The routine computes for given N, P, K:
C
C     PLEK = Prob (X .LE. K)
C     PGTK = Prob (X .GT. K)
C     PEQK = Prob (X .EQ. K)
C
C     Reference: L. Knusel, Computation of the Chi-square and Poisson
C     distribution, SIAM J. Sci. Statist. Comput. 7, pp 1022-1036, 1986.
C
C     .. Parameters ..
      DOUBLE PRECISION  VARMAX
      PARAMETER         (VARMAX=1.0D6)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01BJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, PEQK, PGTK, PLEK
      INTEGER           IFAIL, K, N
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, PC, RLAM
      INTEGER           IER, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01BJZ
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, DBLE, SQRT
C     .. Executable Statements ..
C
      IF (N.LT.0) THEN
         IER = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (P.LE.0.0D0 .OR. P.GE.1.0D0) THEN
         IER = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) P
      ELSE IF (K.LT.0 .OR. K.GT.N) THEN
         IER = 3
         NREC = 1
         WRITE (P01REC,FMT=99997) K, N
      ELSE IF (DBLE(N)*X02AJF().GT.1.0D0) THEN
         IER = 4
         NREC = 1
         WRITE (P01REC,FMT=99996) N
      ELSE IF (N*P*(1.0D0-P).GT.VARMAX) THEN
         IER = 5
         NREC = 2
         WRITE (P01REC,FMT=99995) VARMAX, N*P*(1.0D0-P), N, P
      ELSE
         IER = 0
C
         RLAM = N*P
         IF (RLAM.LT.SQRT(X02AMF())) THEN
C           This case is treated here in order to reduce underflow
C           problems
            PLEK = 1.0D0
            IF (K.EQ.0) THEN
               PEQK = 1.0D0
               PGTK = RLAM
            ELSE IF (K.EQ.1) THEN
               PEQK = RLAM
               PGTK = 0.0D0
            ELSE
               PEQK = 0.0D0
               PGTK = 0.0D0
            END IF
         ELSE IF (K.EQ.N) THEN
            PLEK = 1.0D0
            PGTK = 0.0D0
            PEQK = N*LOG(P)
            IF (PEQK.LT.LOG(X02AMF())) THEN
               PEQK = 0.0D0
            ELSE
               PEQK = EXP(PEQK)
            END IF
         ELSE
            A = DBLE(K+1)
            B = DBLE(N-K)
            PC = 1.0D0 - P
            IF (P.LE.0.5D0) THEN
               CALL G01BJZ(P,A,B,PGTK,PLEK,PEQK)
            ELSE
               CALL G01BJZ(PC,B,A,PLEK,PGTK,PEQK)
            END IF
            PEQK = PEQK*PC/B
         END IF
         IFAIL = 0
      END IF
      IF (IER.NE.0) IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.1: N =',I16)
99998 FORMAT (' ** On entry, P.le.0.0 or P.ge.1.0: P =',1P,D13.5)
99997 FORMAT (' ** On entry, K.lt.0 or K.gt.N: K =',I16,'  N =',I16)
99996 FORMAT (' ** On entry, N is too large: N =',I16)
99995 FORMAT (' ** On entry, the variance exceeds',1P,D10.1,': varianc',
     *       'e =',D13.6,/'              N =',I16,'  P =',D13.5)
      END
