      SUBROUTINE G10ZAF(WEIGHT,N,X,Y,WT,NORD,XORD,YORD,WWT,RSS,IWRK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C      Subroutine G10ZAF is an auxiliary to spline fitting routines,
C     given an unordered data set, weighted or non-weighted, it
C     will order the observations and reweight them if necessary.
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G10ZAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS
      INTEGER           IFAIL, N, NORD
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  WT(*), WWT(N), X(N), XORD(N), Y(N), YORD(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, WTSUM, WTSUM1, YK, YSUM
      INTEGER           I, IERROR, IFAIL2, IJ, J, K, M
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01DAF, M01ZAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IFAIL2 = 0
      IERROR = 0
C
      IF (N.LE.0) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) N
      ELSE IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
C
C        Sort X and corresponding Y and WT
C        Remove observations with zero weight
C
         CALL M01DAF(X,1,N,'A',IWRK,IFAIL2)
         CALL M01ZAF(IWRK,1,N,IFAIL2)
         J = 0
         DO 20 I = 1, N
            IJ = IWRK(I)
            IF (WT(IJ).GT.0.0D0) THEN
               J = J + 1
               XORD(J) = X(IJ)
               YORD(J) = Y(IJ)
               WWT(J) = WT(IJ)
            ELSE IF (WT(IJ).LT.0.0D0) THEN
               IERROR = 2
               WRITE (REC,FMT=99997)
               GO TO 100
            END IF
   20    CONTINUE
         M = J
C
C        Calculate  average Y's and weight vector
C
         IF (M.EQ.0) THEN
            IERROR = 2
            WRITE (REC,FMT=99996)
            I = 0
         ELSE
            I = 1
            RSS = 0.0D0
            YSUM = 0.0D0
            WTSUM = 0.0D0
            DO 40 K = 1, M - 1
               IF (WTSUM.EQ.0.0D0) THEN
                  SUM = YORD(K)
                  WTSUM = WWT(K)
               ELSE
                  YK = YORD(K)
                  WTSUM1 = WTSUM
                  WTSUM = WTSUM + WWT(K)
                  RSS = RSS + (YK-SUM)*(YK-SUM)*WWT(K)*WTSUM1/WTSUM
                  SUM = SUM + WWT(K)*(YK-SUM)/WTSUM
               END IF
               IF (XORD(K+1).NE.XORD(K)) THEN
                  WWT(I) = WTSUM
                  XORD(I) = XORD(K)
                  YORD(I) = SUM
                  I = I + 1
                  YSUM = 0.0D0
                  WTSUM = 0.0D0
               END IF
   40       CONTINUE
C
C           Set end values
C
            XORD(I) = XORD(M)
            IF (WTSUM.EQ.0) THEN
               YORD(I) = YORD(M)
               WWT(I) = WWT(M)
            ELSE
               YK = YORD(M)
               WTSUM1 = WTSUM
               WWT(I) = WTSUM + WWT(M)
               RSS = RSS + (YK-SUM)*(YK-SUM)*WWT(M)*WTSUM1/WWT(I)
               YORD(I) = SUM + WWT(M)*(YK-SUM)/WWT(I)
            END IF
         END IF
C
      ELSE IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
C
C        Sort X and corresponding Y
C
         CALL M01DAF(X,1,N,'A',IWRK,IFAIL2)
         DO 60 I = 1, N
            XORD(IWRK(I)) = X(I)
            YORD(IWRK(I)) = Y(I)
   60    CONTINUE
C
         I = 1
         J = 0
         SUM = 0.0D0
         RSS = 0.0D0
         DO 80 K = 1, N - 1
            IF (J.EQ.0) THEN
               SUM = YORD(K)
               J = 1
            ELSE
               YK = YORD(K)
               RSS = RSS + (YK-SUM)*(YK-SUM)*DBLE(J)/DBLE(J+1)
               J = J + 1
               SUM = SUM + (YK-SUM)/DBLE(J)
            END IF
            IF (XORD(K+1).NE.XORD(K)) THEN
               WWT(I) = DBLE(J)
               XORD(I) = XORD(K)
               YORD(I) = SUM
               I = I + 1
               J = 0
            END IF
   80    CONTINUE
C
C        Set end values
C
         WWT(I) = DBLE(J+1)
         XORD(I) = XORD(N)
         IF (J.EQ.0) THEN
            YORD(I) = YORD(N)
         ELSE
            YK = YORD(N)
            RSS = RSS + (YK-SUM)*(YK-SUM)*DBLE(J)/DBLE(J+1)
            J = J + 1
            YORD(I) = SUM + (YK-SUM)/DBLE(J)
         END IF
      ELSE
         IERROR = 1
         WRITE (REC,FMT=99998) WEIGHT
         I = 0
      END IF
      NORD = I
  100 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, N.le.0: N = ',I16)
99998 FORMAT (1X,'** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99997 FORMAT (1X,'** On entry, at least one weight is negative.')
99996 FORMAT (1X,'** On entry, all weights are zero.')
      END
