      SUBROUTINE G08AGF(N,X,XME,TAIL,ZEROS,RS,RSNOR,P,NZ1,WRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-922 (APR 1991).
C
C     G08AGF performs the one-sample Wilcoxon signed rank test.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, RS, RSNOR, XME
      INTEGER           IFAIL, N, NZ1
      CHARACTER*1       TAIL, ZEROS
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(3*N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  NWN, RMED, W, WD, WN, XF
      INTEGER           HIGH, I, IERROR, IFAIL2, IV, J, NZ
      LOGICAL           LOWER, NOZERO, TOTAIL, UPPER
C     .. Local Arrays ..
      INTEGER           IR(80)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF
      INTEGER           P01ABF
      EXTERNAL          S15ABF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G08AGY, G08AGZ, M01CBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, NINT, DBLE, SIGN, SQRT
C     .. Executable Statements ..
C
C     Parameter check
C
      TOTAIL = TAIL .EQ. 'T' .OR. TAIL .EQ. 't'
      UPPER = TAIL .EQ. 'U' .OR. TAIL .EQ. 'u'
      LOWER = TAIL .EQ. 'L' .OR. TAIL .EQ. 'l'
      NOZERO = ZEROS .EQ. 'N' .OR. ZEROS .EQ. 'n'
      IF ( .NOT. TOTAIL .AND. .NOT. UPPER .AND. .NOT. LOWER) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) TAIL
      ELSE IF (ZEROS.NE.'Y' .AND. ZEROS.NE.'y' .AND. .NOT. NOZERO) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) ZEROS
      ELSE IF (N.LT.1) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99997) N
      ELSE
C
C        First subtract XME, then store the sign, rank the absolute
C        values and then compute the Wilcoxon Rank-Sum Statistic, RS.
C
         IERROR = 0
         NZ = 0
         DO 20 I = 1, N
            WRK(I) = X(I) - XME
            IF (WRK(I).NE.0.0D0) THEN
               WRK(2*N+I) = SIGN(1.0D0,WRK(I))
            ELSE
               WRK(2*N+I) = 0.0D0
               NZ = NZ + 1
            END IF
            WRK(I) = ABS(WRK(I))
   20    CONTINUE
         IF (NZ.EQ.N) THEN
            IERROR = 3
            WRITE (P01REC,FMT=99996)
         ELSE
            XF = 0.0D0
            CALL G08AGY(WRK(1),WRK(N+1),N,1,XF)
            XF = XF - DBLE(NZ*(NZ**2-1))/12.0D0
            NZ1 = N - NZ
            RS = 0.0D0
            DO 40 I = 1, N
               IF (NOZERO) WRK(N+I) = WRK(N+I) - DBLE(NZ)
               IF (WRK(2*N+I).GT.0.0D0) RS = RS + WRK(N+I)
   40       CONTINUE
C
C              Compute normalized test statistic
C
            IF (NOZERO) THEN
               RMED = DBLE(NZ1*(NZ1+1)/2)
               WN = 2.0D0*RS - RMED
               WD = DBLE(NZ1*(NZ1+1))*(DBLE(2*NZ1+1)/6.0D0) - XF
            ELSE
               RMED = DBLE((N*(N+1)-NZ*(NZ+1))/2)
               WN = 2.0D0*RS - RMED
               WD = DBLE(N*(N+1))*(DBLE(2*N+1)/6.0D0) - DBLE(NZ*(NZ+1))
     *              *(DBLE(2*NZ+1)/6.0D0) - XF
            END IF
C
C           Check whether statistic is above or below the expected
C           statistic under the null hypothesis and apply
C           continuity correction
C
            IF (WN.GT.0.0D0) THEN
               NWN = WN - 1.0D0
            ELSE IF (WN.LT.0.0D0) THEN
               NWN = WN + 1.0D0
            ELSE
               NWN = WN
            END IF
            WD = SQRT(WD)
            RSNOR = NWN/WD
C
C           Approx significance level for large samples.
C           This depends on the choice of tail and because
C           the normalized correction above is always corrected
C           for to continuity TOWARDS the centre we may have to
C           re-adjust it before computing the significance level.
C
            IFAIL2 = 0
            IF (NZ1.GT.80) THEN
               IF (LOWER .AND. WN.LT.0.0D0) THEN
                  P = S15ABF(RSNOR,IFAIL2)
               ELSE IF (UPPER .AND. WN.GT.0.0D0) THEN
                  P = 1.0D0 - S15ABF(RSNOR,IFAIL2)
               ELSE IF (LOWER .AND. WN.GE.0.0D0) THEN
                  IF (RS.EQ.RMED) THEN
                     P = 1.0D0
                  ELSE
                     NWN = WN + 1.0D0
                     W = NWN/WD
                     P = S15ABF(W,IFAIL2)
                  END IF
               ELSE IF (UPPER .AND. WN.LE.0.0D0) THEN
                  IF (RS.EQ.0.0D0) THEN
                     P = 1.0D0
                  ELSE
                     NWN = WN - 1.0D0
                     W = NWN/WD
                     P = 1.0D0 - S15ABF(W,IFAIL2)
                  END IF
               ELSE IF (TOTAIL) THEN
                  P = S15ABF(RSNOR,IFAIL2)
                  P = 2.0D0*(MIN(P,1.0D0-P))
               END IF
            ELSE
C
C              Exact significance level for small samples.
C
               J = 0
               IF (XF.EQ.0.0D0) THEN
                  DO 60 I = 1, N
                     IF (WRK(2*N+I).NE.0.D0) THEN
                        J = J + 1
                        IR(J) = NINT(WRK(N+I))
                     END IF
   60             CONTINUE
                  IV = NINT(RS)
                  HIGH = NINT(RMED)
               ELSE
                  DO 80 I = 1, N
                     IF (WRK(2*N+I).NE.0.D0) THEN
                        J = J + 1
                        IR(J) = NINT(2.0D0*WRK(N+I))
                     END IF
   80             CONTINUE
                  IV = NINT(2.0D0*RS)
                  HIGH = NINT(2.0D0*RMED)
               END IF
               IFAIL2 = 0
               CALL M01CBF(IR,1,NZ1,'A',IFAIL2)
C
               IF (LOWER) THEN
                  IF (IV.EQ.HIGH) THEN
                     P = 1.0D0
                  ELSE IF ((2*IV).LE.HIGH) THEN
                     CALL G08AGZ(NZ1,IR,IV,P)
                  ELSE
                     IV = HIGH - IV - 1
                     CALL G08AGZ(NZ1,IR,IV,P)
                     P = 1.0D0 - P
                  END IF
               ELSE IF (UPPER) THEN
                  IF (IV.EQ.0) THEN
                     P = 1.0D0
                  ELSE IF ((2*IV).GT.HIGH) THEN
                     IV = HIGH - IV
                     CALL G08AGZ(NZ1,IR,IV,P)
                  ELSE
                     IV = IV - 1
                     CALL G08AGZ(NZ1,IR,IV,P)
                     P = 1.0D0 - P
                  END IF
               ELSE IF (TOTAIL) THEN
                  IF ((2*IV).LE.HIGH) THEN
                     CALL G08AGZ(NZ1,IR,IV,P)
                  ELSE
                     IV = HIGH - IV
                     CALL G08AGZ(NZ1,IR,IV,P)
                  END IF
                  P = 2.0D0*P
                  IF (P.GT.1.0D0) P = 1.0D0
               END IF
            END IF
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
99999 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99998 FORMAT (1X,'** On entry, ZEROS is not valid: ZEROS = ',A1)
99997 FORMAT (1X,'** On entry, N.lt.1 : N = ',I16)
99996 FORMAT (1X,'** All elements of the sample equal XME, i.e. (Varia',
     *       'nce = 0)')
      END
