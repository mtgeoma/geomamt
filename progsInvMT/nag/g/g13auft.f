      SUBROUTINE G13AUF(N,Z,M,K,RS,Y,MEAN,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13AUF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, K, M, N
      CHARACTER*1       RS
C     .. Array Arguments ..
      DOUBLE PRECISION  MEAN(K), Y(K), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  MEANI, SUM, TEMP, ZMAX, ZMIN
      INTEGER           I, IERROR, IM, J, NREC
      LOGICAL           RANGE
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, DBLE, SQRT
C     .. Executable Statements ..
C
C     This routine divides a time series into consecutive groups
C     of M (or more) observations and calculates the range or standard
C     deviation and mean for each group.  These values may then be used
C     to produce a range-mean plot or a standard deviation-mean plot.
C
C     Test for errors in the input data
C
      NREC = 1
      IERROR = 0
      IF (M.LT.2) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) M
      ELSE IF (N.LT.M) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) N, M
      ELSE IF (K.NE.N/M) THEN
         IERROR = 1
         NREC = 2
         WRITE (P01REC,FMT=99997) N, M, K
      ELSE IF (RS.NE.'R' .AND. RS.NE.'r' .AND. RS.NE.'S' .AND. RS.NE.
     *         's') THEN
         IERROR = 2
         WRITE (P01REC,FMT=99996) RS
      ELSE
         NREC = 0
         RANGE = RS .EQ. 'R' .OR. RS .EQ. 'r'
         IF (RANGE) THEN
            DO 40 I = 1, K
               IM = N - M*(I-1)
               SUM = Z(IM)
               ZMIN = SUM
               ZMAX = SUM
               DO 20 J = 2, M
                  TEMP = Z(IM-J+1)
                  SUM = SUM + TEMP
                  ZMIN = MIN(ZMIN,TEMP)
                  ZMAX = MAX(ZMAX,TEMP)
   20          CONTINUE
               MEAN(K-I+1) = SUM/DBLE(M)
               Y(K-I+1) = ZMAX - ZMIN
   40       CONTINUE
         ELSE
            DO 100 I = 1, K
               IM = N - M*(I-1)
               SUM = 0.0D0
               DO 60 J = 1, M
                  SUM = SUM + Z(IM-J+1)
   60          CONTINUE
               MEANI = SUM/DBLE(M)
               MEAN(K-I+1) = MEANI
               SUM = 0.0D0
               DO 80 J = 1, M
                  TEMP = Z(IM-J+1) - MEANI
                  SUM = SUM + TEMP*TEMP
   80          CONTINUE
               Y(K-I+1) = SQRT(SUM/DBLE(M-1))
  100       CONTINUE
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, M.lt.2 : M = ',I16)
99998 FORMAT (1X,'** On entry, N.lt.M : N = ',I16,' and M = ',I16)
99997 FORMAT (1X,'** On entry, NGRPS does not equal N/M :',/5X,'N = ',
     *       I16,' M = ',I16,' and NGRPS = ',I16)
99996 FORMAT (1X,'** On entry, RS is not valid : RS = ',A1)
      END
