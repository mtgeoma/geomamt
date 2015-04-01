      SUBROUTINE G08ALF(N,K,X,LDX,Q,PROB,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Performs the Cochran Q test.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08ALF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PROB, Q
      INTEGER           IFAIL, K, LDX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(LDX,K)
C     .. Local Scalars ..
      DOUBLE PRECISION  CSUM2, DF, DIV, RK, RSUM2, SUM, TSUM, TX
      INTEGER           I, IERROR, IF2, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF
      INTEGER           P01ABF
      EXTERNAL          G01ECF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      NREC = 0
      IERROR = 0
      IF (N.LT.2) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (K.LT.2) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99998) K
      ELSE IF (LDX.LT.N) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99997) LDX, N
      ELSE
         SUM = 0.0D0
         CSUM2 = 0.0D0
         DO 40 J = 1, K
            TSUM = 0.0D0
            DO 20 I = 1, N
               TX = X(I,J)
               IF (TX.EQ.0.0D0 .OR. TX.EQ.1.0D0) THEN
                  TSUM = TSUM + TX
               ELSE
                  IERROR = 2
                  NREC = 1
                  WRITE (P01REC,FMT=99996)
                  GO TO 100
               END IF
   20       CONTINUE
            SUM = SUM + TSUM
            CSUM2 = CSUM2 + TSUM*TSUM
   40    CONTINUE
         RSUM2 = 0.0D0
         DO 80 I = 1, N
            TSUM = 0.0D0
            DO 60 J = 1, K
               TSUM = TSUM + X(I,J)
   60       CONTINUE
            RSUM2 = RSUM2 + TSUM*TSUM
   80    CONTINUE
         RK = DBLE(K)
         DF = RK - 1.0D0
         DIV = RK*SUM - RSUM2
         IF (DIV.GT.0.0D0) THEN
            Q = DF*(RK*CSUM2-SUM*SUM)/DIV
            IF2 = 1
            PROB = G01ECF('UPPER',Q,DF,IF2)
            IF (IF2.EQ.4) THEN
               IERROR = 3
               NREC = 2
               WRITE (P01REC,FMT=99995)
            END IF
         ELSE
            Q = 0.0D0
            PROB = 1.0D0
         END IF
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (1X,'** On entry, K.lt.2 : K = ',I16)
99997 FORMAT (1X,'** On entry, LDX.lt.N : LDX = ',I16,' and N = ',I16)
99996 FORMAT (1X,'** An element of X is not equal to 0 or 1.')
99995 FORMAT (1X,'** The solution has failed to converge while calcula',
     *       'ting the tail probability.',/'    The result may be a re',
     *       'asonable approximation.')
      END
