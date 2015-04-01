      SUBROUTINE G02BTF(MEAN,M,WT,X,INCX,SW,XBAR,C,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G02BTF RETURNS THE SUM OF MEAN SQUARED DEVIATIONS
C     AND THE SUM OF SQUARES OF A DATA SET.
C     THESE VALUES ARE CALCULATED USING A STABLE UPDATING METHOD
C     WHICH REVISES THE EXISTING VALUES IN A SINGLE PASS.
C
C     BASED ON ALGORITHM AS (WV2) COMM. ACM., (1979), VOL.22, NO.9
C
C     USES NAG LIBRARY ROUTINE P01ABF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      DOUBLE PRECISION  ZERO
      PARAMETER         (SRNAME='G02BTF',ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SW, WT
      INTEGER           IFAIL, INCX, M
      CHARACTER         MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  C((M*M+M)/2), X(M*INCX), XBAR(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  RW, RWW, SW1
      INTEGER           I, IERROR, II, IJ, J
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IF (INCX.LT.1) THEN
         IERROR = 1
         WRITE (REC,FMT=99996) INCX
      ELSE IF (M.LT.1) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) M
      ELSE IF (SW.LT.ZERO) THEN
         IERROR = 2
         WRITE (REC,FMT=99998) 'SW', 'SW', SW
      ELSE IF ((SW+WT).LT.ZERO) THEN
         IERROR = 3
         WRITE (REC,FMT=99998) '(SW+WT)', '(SW+WT)', (SW+WT)
      ELSE
         IERROR = 0
         IF (INCX.GT.0) THEN
            II = 1 - INCX
         ELSE
            II = 1 - M*INCX
         END IF
         IF (WT.NE.ZERO) THEN
            IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
C
C              Calculate the Sum of Mean Squared Deviations
C
               IF (SW.EQ.0) THEN
                  SW = WT
                  DO 20 I = 1, M
                     XBAR(I) = X(I*INCX+II)
                     C(I) = ZERO
   20             CONTINUE
                  DO 40 I = M + 1, (M*M+M)/2
                     C(I) = ZERO
   40             CONTINUE
               ELSE
                  SW1 = SW + WT
                  IF (SW1.EQ.0) THEN
                     DO 60 I = 1, M
                        XBAR(I) = ZERO
                        C(I) = ZERO
   60                CONTINUE
                     DO 80 I = M + 1, (M*M+M)/2
                        C(I) = ZERO
   80                CONTINUE
                  ELSE
                     RW = WT/SW1
                     RWW = RW*SW
                     SW = SW1
                     IJ = 0
                     DO 120 I = 1, M
                        DO 100 J = 1, I
                           IJ = IJ + 1
                           C(IJ) = (XBAR(I)-X(I*INCX+II))*RWW*(XBAR(J)
     *                             -X(J*INCX+II)) + C(IJ)
  100                   CONTINUE
  120                CONTINUE
                     DO 140 I = 1, M
                        XBAR(I) = (X(I*INCX+II)-XBAR(I))*RW + XBAR(I)
  140                CONTINUE
                  END IF
               END IF
            ELSE IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
C
C              Calculate the Sum of Squares
C
               IF (SW.EQ.0) THEN
                  SW = WT
                  IJ = 0
                  DO 180 I = 1, M
                     XBAR(I) = X(I*INCX+II)
                     DO 160 J = 1, I
                        IJ = IJ + 1
                        C(IJ) = X(I*INCX+II)*WT*X(J*INCX+II)
  160                CONTINUE
  180             CONTINUE
               ELSE
                  SW = SW + WT
                  IF (SW.EQ.0) THEN
                     DO 200 I = 1, M
                        XBAR(I) = ZERO
                        C(I) = ZERO
  200                CONTINUE
                     DO 220 I = M + 1, (M*M+M)/2
                        C(I) = ZERO
  220                CONTINUE
                  ELSE
                     RW = WT/SW
                     IJ = 0
                     DO 260 I = 1, M
                        DO 240 J = 1, I
                           IJ = IJ + 1
                           C(IJ) = X(I*INCX+II)*WT*X(J*INCX+II) + C(IJ)
  240                   CONTINUE
  260                CONTINUE
                     DO 280 I = 1, M
                        XBAR(I) = (X(I*INCX+II)-XBAR(I))*RW + XBAR(I)
  280                CONTINUE
                  END IF
               END IF
            ELSE
               IERROR = 4
               WRITE (REC,FMT=99997) MEAN
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
99999 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99998 FORMAT (' ** On entry, ',A,'.lt.0.0 : ',A,' = ',D13.5)
99997 FORMAT (' ** On entry, MEAN is not valid : MEAN = ',A1)
99996 FORMAT (' ** On entry, INCX.lt.1 : INCX = ',I16)
      END
