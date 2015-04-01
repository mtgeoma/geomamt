      SUBROUTINE G02DCF(UPDATE,MEAN,WEIGHT,M,ISX,Q,LDQ,IP,X,IX,Y,WT,RSS,
     *                  WK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-919 (APR 1991).
C
C     ADD (UPDATE='A') OR DROPS (UPDATE='D') AN OBSERVATION
C     FROM A LINEAR REGRESSION RETURNING THE UPDATED RSS AND NEW R
C     MATRIX AND Q'Y MATRIX STORED IN Q (SEE G02DAF)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02DCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RSS, WT, Y
      INTEGER           IFAIL, IP, IX, LDQ, M
      CHARACTER         MEAN, UPDATE, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*), WK(3*IP), X(*)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  QYI, QYN, SQRWT, ZETA
      INTEGER           I, IERROR, K, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06QQF, G02DCZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) M
      ELSE IF (LDQ.LT.IP) THEN
         WRITE (P01REC(1),FMT=99998) LDQ, IP
      ELSE IF (IX.LT.1) THEN
         WRITE (P01REC(1),FMT=99989) IX
      ELSE IF (RSS.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99993) RSS
      ELSE IF (MEAN.NE.'M' .AND. MEAN.NE.'m' .AND. MEAN.NE.'Z' .AND.
     *         MEAN.NE.'z') THEN
         WRITE (P01REC(1),FMT=99991) MEAN
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         WRITE (P01REC(1),FMT=99997) WEIGHT
      ELSE IF (UPDATE.NE.'A' .AND. UPDATE.NE.'a' .AND. UPDATE.NE.
     *         'D' .AND. UPDATE.NE.'d') THEN
         WRITE (P01REC(1),FMT=99996) UPDATE
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            K = 1
         ELSE
            K = 0
         END IF
         DO 20 I = 1, M
            IF (ISX(I).GT.0) THEN
               K = K + 1
            END IF
   20    CONTINUE
         IF (K.NE.IP) THEN
            IERROR = 1
            NREC = 2
            IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
               WRITE (P01REC,FMT=99990) K, IP
            ELSE
               WRITE (P01REC,FMT=99988) K - 1, IP - 1
            END IF
            GO TO 160
         END IF
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            IF (WT.LT.0.0D0) THEN
               IERROR = 2
               WRITE (P01REC(1),FMT=99992) WT
               GO TO 160
            ELSE IF (WT.EQ.0.0D0) THEN
               GO TO 160
            ELSE
               SQRWT = SQRT(WT)
               QYN = SQRWT*Y
               IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
                  K = 1
                  WK(2*IP+1) = SQRWT
               ELSE
                  K = 0
               END IF
               DO 40 I = 1, M
                  IF (ISX(I).GT.0) THEN
                     K = K + 1
                     WK(2*IP+K) = SQRWT*X(IX*(I-1)+1)
                  END IF
   40          CONTINUE
            END IF
         ELSE
            QYN = Y
            IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
               K = 1
               WK(2*IP+1) = 1.0D0
            ELSE
               K = 0
            END IF
            DO 60 I = 1, M
               IF (ISX(I).GT.0) THEN
                  K = K + 1
                  WK(2*IP+K) = X(IX*(I-1)+1)
               END IF
   60       CONTINUE
         END IF
         IF (UPDATE.EQ.'A' .OR. UPDATE.EQ.'a') THEN
C
C           UPDATE R
C
            CALL F06QQF(IP,1.0D0,WK(2*IP+1),1,Q(1,2),LDQ,WK(1),WK(IP+1))
C
C           UPDATE Q'Y
C
            DO 80 I = 1, IP
               QYI = Q(I,1)
               Q(I,1) = WK(IP+I)*QYN + WK(I)*QYI
               QYN = WK(I)*QYN - WK(IP+I)*QYI
   80       CONTINUE
C
C           UPDATE RSS
C
            RSS = QYN*QYN + RSS
         ELSE
            DO 100 I = 1, IP
               IF (Q(I,I+1).EQ.0.0D0) GO TO 140
  100       CONTINUE
C
C           DOWNDATE R
C
            CALL G02DCZ(IP,1.0D0,WK(2*IP+1),1,Q(1,2),LDQ,WK(1),WK(IP+1),
     *                  ZETA)
            IF (ZETA.GE.0.0D0) THEN
C
C              DOWNDATE Q'Y
C
               ZETA = QYN
               DO 120 I = 1, IP
                  Q(I,1) = (Q(I,1)-WK(IP+I)*ZETA)/WK(I)
                  ZETA = WK(I)*ZETA - WK(IP+I)*Q(I,1)
  120          CONTINUE
C
C              DOWNDATE RSS
C
               ZETA = ZETA*ZETA
               IF (ZETA.GT.RSS) THEN
                  IERROR = 4
                  WRITE (P01REC(1),FMT=99994)
                  RSS = 0.0D0
               ELSE
                  RSS = RSS - ZETA
               END IF
               GO TO 160
            END IF
  140       IERROR = 3
            WRITE (P01REC(1),FMT=99995)
         END IF
      END IF
  160 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry IP .lt. 1 : IP = ',I16)
99998 FORMAT (' ** On entry LDQ .lt. IP : LDQ = ',I16,' IP = ',I16)
99997 FORMAT (' ** On entry WEIGHT does not have a valid value: WEIGHT',
     *       ' = ',A1)
99996 FORMAT (' ** On entry UPDATE does not have a valid value: UPDATE',
     *       ' = ',A1)
99995 FORMAT (' ** The R matrix could not be updated')
99994 FORMAT (' ** The RSS could not be updated')
99993 FORMAT (' ** On entry RSS .lt. 0.0: RSS = ',D13.5)
99992 FORMAT (' ** On entry WT .lt. 0.0: WT = ',D13.5)
99991 FORMAT (' ** On entry MEAN does not have a valid value: MEAN = ',
     *       A1)
99990 FORMAT (' ** On entry ',I16,' elements of ISX .gt. 0',/'        ',
     *       '     instead of IP = ',I16)
99989 FORMAT (' ** On entry IX .lt. 1 : IX = ',I16)
99988 FORMAT (' ** On entry ',I16,' elements of ISX .gt. 0',/'        ',
     *       '     instead of IP-1 (for mean) = ',I16)
      END
