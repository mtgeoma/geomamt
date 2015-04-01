      SUBROUTINE F01LHZ(BLKSTR,NBLOKS,N,LENARQ,PASNAM,REPORT,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IERR, LENARQ, N, NBLOKS
      LOGICAL           REPORT
      CHARACTER*6       PASNAM
C     .. Array Arguments ..
      INTEGER           BLKSTR(3,NBLOKS)
C     .. Local Scalars ..
      INTEGER           IDEV, JOVER, K, KCOLS, KOVER, KROWS, NCOLS,
     *                  NCOLS1, NROWS
      CHARACTER*80      REC
C     .. External Subroutines ..
      EXTERNAL          X04AAF, X04BAF
C     .. Executable Statements ..
C
      CALL X04AAF(0,IDEV)
C
      IF (NBLOKS.LT.1) THEN
         IERR = -3
         IF (REPORT) THEN
            WRITE (REC,FMT=99999) PASNAM, NBLOKS
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      IF (N.LT.1) THEN
         IERR = -2
         IF (REPORT) THEN
            WRITE (REC,FMT=99998) PASNAM, N
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      IF (IERR.LT.0) GO TO 60
      IF (N.LT.NBLOKS) THEN
         IERR = -1
         IF (REPORT) THEN
            WRITE (REC,FMT=99997) PASNAM, N, NBLOKS
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      DO 20 K = 1, NBLOKS - 1
         KROWS = BLKSTR(1,K)
         KCOLS = BLKSTR(2,K)
         KOVER = BLKSTR(3,K)
         IF (KROWS.LT.1) THEN
            IERR = -4
            IF (REPORT) THEN
               WRITE (REC,FMT=99996) PASNAM, K, KROWS
               CALL X04BAF(IDEV,REC)
            END IF
         END IF
         IF (KCOLS.LT.1) THEN
            IERR = -4
            IF (REPORT) THEN
               WRITE (REC,FMT=99995) PASNAM, K, KCOLS
               CALL X04BAF(IDEV,REC)
            END IF
         END IF
         IF (KOVER.LT.0) THEN
            IERR = -4
            IF (REPORT) THEN
               WRITE (REC,FMT=99994) PASNAM, K, KOVER
               CALL X04BAF(IDEV,REC)
            END IF
         END IF
   20 CONTINUE
      KROWS = BLKSTR(1,NBLOKS)
      KCOLS = BLKSTR(2,NBLOKS)
      IF (KROWS.LT.1) THEN
         IERR = -4
         IF (REPORT) THEN
            WRITE (REC,FMT=99996) PASNAM, NBLOKS, KROWS
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      IF (KCOLS.LT.1) THEN
         IERR = -4
         IF (REPORT) THEN
            WRITE (REC,FMT=99995) PASNAM, NBLOKS, KCOLS
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      IF (IERR.LT.0) GO TO 60
      KROWS = BLKSTR(1,1)
      KCOLS = BLKSTR(2,1)
      KOVER = BLKSTR(3,1)
      IF (KCOLS.LT.KROWS) THEN
         IERR = -5
         IF (REPORT) THEN
            WRITE (REC,FMT=99988) PASNAM, KCOLS, KROWS
            CALL X04BAF(IDEV,REC)
         END IF
         GO TO 60
      ELSE IF (KCOLS-KOVER.GT.KROWS) THEN
         IERR = -5
         IF (REPORT) THEN
            WRITE (REC,FMT=99987) PASNAM, KCOLS - KOVER, KROWS
            CALL X04BAF(IDEV,REC)
         END IF
         GO TO 60
      END IF
      LENARQ = KROWS*KCOLS
      NROWS = KROWS
      NCOLS = KCOLS
      JOVER = KOVER
      DO 40 K = 2, NBLOKS - 1
         KROWS = BLKSTR(1,K)
         KCOLS = BLKSTR(2,K)
         KOVER = BLKSTR(3,K)
         LENARQ = LENARQ + KROWS*KCOLS
         NROWS = NROWS + KROWS
         NCOLS = NCOLS + KCOLS - JOVER
         NCOLS1 = NCOLS - KOVER
         IF ( .NOT. (NCOLS1.LE.NROWS .AND. NROWS.LE.NCOLS)) THEN
            IERR = -5
            IF (REPORT) THEN
               WRITE (REC,FMT=99986) PASNAM, K
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99985)
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99984)
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99983)
               CALL X04BAF(IDEV,REC)
            END IF
            GO TO 60
         ELSE IF (KOVER+JOVER.GT.KCOLS) THEN
            IERR = -5
            IF (REPORT) THEN
               WRITE (REC,FMT=99993) PASNAM, K - 1, K, JOVER + KOVER
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99992) K, KCOLS
               CALL X04BAF(IDEV,REC)
            END IF
            GO TO 60
         END IF
         JOVER = KOVER
   40 CONTINUE
      IF (NBLOKS.GT.1) THEN
         KROWS = BLKSTR(1,NBLOKS)
         KCOLS = BLKSTR(2,NBLOKS)
         NROWS = NROWS + KROWS
         NCOLS = NCOLS + KCOLS - JOVER
         LENARQ = LENARQ + KROWS*KCOLS
      END IF
      IF (NROWS.NE.N) THEN
         IERR = -5
         IF (REPORT) THEN
            WRITE (REC,FMT=99991) PASNAM
            CALL X04BAF(IDEV,REC)
            WRITE (REC,FMT=99990)
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
      IF (NCOLS.NE.N) THEN
         IERR = -5
         IF (REPORT) THEN
            WRITE (REC,FMT=99991) PASNAM
            CALL X04BAF(IDEV,REC)
            WRITE (REC,FMT=99989)
            CALL X04BAF(IDEV,REC)
         END IF
      END IF
   60 CONTINUE
      RETURN
C
99999 FORMAT (' ** ',A6,' - NBLOKS(=',I16,') .le. 0 **')
99998 FORMAT (' ** ',A6,' - N(=',I16,') .le. 0 **')
99997 FORMAT (' ** ',A6,' - N(=',I16,') .lt. NBLOKS(=',I16,') **')
99996 FORMAT (' ** ',A6,' - BLKSTR(1,',I6,') (=',I16,') .lt. 1 **')
99995 FORMAT (' ** ',A6,' - BLKSTR(2,',I6,') (=',I16,') .lt. 1 **')
99994 FORMAT (' ** ',A6,' - BLKSTR(3,',I6,') (=',I16,') .lt. 0 **')
99993 FORMAT (' ** ',A6,' - BLKSTR(3,',I6,')+BLKSTR(3,',I6,')(=',I16,
     *       ')')
99992 FORMAT (13X,'.lt. BLKSTR(2,',I6,') (=',I16,') **')
99991 FORMAT (' ** ',A6,' - the following equality does not hold')
99990 FORMAT ('     sum( BLKSTR(1,k) :k=1,NBLOKS) .eq. N')
99989 FORMAT ('     BLKSTR(2,1) + sum(BLKSTR(2,k)-BLKSTR(3,k-1) :k=2,N',
     *       'BLOKS) .eq. N **')
99988 FORMAT (' ** ',A6,' - BLKSTR(2,1) (=',I6,') .lt. BLKSTR(1,1) (=',
     *       I6,') **')
99987 FORMAT (' ** ',A6,' - BLKSTR(2,1)-BLKSTR(3,1) (=',I6,') .gt. BLK',
     *       'STR(1,1) (=',I6,') **')
99986 FORMAT (' ** ',A6,' - the following inequality was not satisfied',
     *       ' for j (=',I6,')')
99985 FORMAT ('     sum(BLKSTR(2,k)-BLKSTR(3,k):k=1,J)  .le.  ')
99984 FORMAT ('     sum(BLKSTR(1,k):k=1,J)  .le.  ')
99983 FORMAT ('     BLKSTR(2,1)+sum(BLKSTR(2,k)-BLKSTR(3,k-1):k=2,j) **'
     *       )
      END
