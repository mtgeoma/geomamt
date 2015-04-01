      SUBROUTINE G12AAF(N,T,IC,FREQ,IFREQ,ND,TP,P,PSIG,IWK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G12AAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N, ND
      CHARACTER         FREQ
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N), PSIG(N), T(N), TP(N)
      INTEGER           IC(N), IFREQ(*), IWK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           D, I, IERROR, IFAULT, INC, NJ, NJ1, NREC
      LOGICAL           FR
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01DAF, M01ZAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
      IFAULT = 0
C
      IF (N.LT.2) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE
         DO 20 I = 1, N
            IF (IC(I).NE.0 .AND. IC(I).NE.1) THEN
               IERROR = 3
               WRITE (P01REC,FMT=99998) I
               GO TO 120
            END IF
   20    CONTINUE
C
         CALL M01DAF(T,1,N,'A',IWK,IFAULT)
         CALL M01ZAF(IWK,1,N,IFAULT)
C
         ND = 0
         NJ1 = 0
         SUM = 0.0D0
C
         IF (FREQ.EQ.'F' .OR. FREQ.EQ.'f') THEN
            FR = .TRUE.
         ELSE IF (FREQ.EQ.'S' .OR. FREQ.EQ.'s') THEN
            FR = .FALSE.
         ELSE
            IERROR = 2
            WRITE (P01REC,FMT=99996) FREQ
            GO TO 120
         END IF
C
C        If frequencies check values and find total number of
C        observations
C
         IF (FR) THEN
            DO 40 I = 1, N
               IF (IFREQ(I).LT.0) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99997) I
                  GO TO 120
               END IF
               NJ1 = NJ1 + IFREQ(I)
   40       CONTINUE
         ELSE
            NJ1 = N
         END IF
         INC = 1
         I = 0
   60    CONTINUE
C
C        Loop through data computing estimate
C
         NJ = NJ1
         D = 0
C
C        NJ = number still in study
C        D  = number failed at ith time
C
   80    CONTINUE
C
C        Check for ties
C
         I = I + 1
         IF (FR) INC = IFREQ(IWK(I))
         IF (I.GE.N) THEN
            GO TO 100
         ELSE IF (T(IWK(I)).EQ.T(IWK(I+1))) THEN
            NJ1 = NJ1 - INC
            IF (IC(IWK(I)).EQ.0) D = D + INC
            GO TO 80
         END IF
         NJ1 = NJ1 - INC
         IF (IC(IWK(I)).EQ.0) D = D + INC
         IF (D.NE.0) THEN
            ND = ND + 1
            IF (ND.EQ.1) THEN
               P(ND) = DBLE(NJ-D)/DBLE(NJ)
            ELSE
               P(ND) = P(ND-1)*DBLE(NJ-D)/DBLE(NJ)
            END IF
            TP(ND) = T(IWK(I))
            SUM = SUM + DBLE(D)/(DBLE(NJ)*DBLE(NJ-D))
            PSIG(ND) = SQRT((P(ND)**2)*SUM)
         END IF
         GO TO 60
C
C        Deal with final point
C
  100    IF (IC(IWK(I)).EQ.0) D = D + INC
         IF (D.GT.0) THEN
            IF ((NJ-D).GT.0) THEN
               ND = ND + 1
               IF (ND.EQ.1) THEN
                  P(ND) = DBLE(NJ-D)/DBLE(NJ)
               ELSE
                  P(ND) = P(ND-1)*DBLE(NJ-D)/DBLE(NJ)
               END IF
               TP(ND) = T(IWK(I))
               SUM = SUM + DBLE(D)/(DBLE(NJ)*DBLE(NJ-D))
               PSIG(ND) = SQRT((P(ND)**2)*SUM)
            ELSE IF ((NJ-D).EQ.0) THEN
               ND = ND + 1
               P(ND) = 0.0D0
               TP(ND) = T(IWK(I))
               PSIG(ND) = 0.0D0
            END IF
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.2 : N =',I16)
99998 FORMAT (' ** On entry, the ',I16,'th element of IC .ne. 0 or 1.')
99997 FORMAT (' ** On entry, the ',I16,'th element of IFREQ .lt. 0.')
99996 FORMAT (' ** On entry, FREQ is not valid: FREQ = ',A1)
      END
