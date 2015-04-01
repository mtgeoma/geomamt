      SUBROUTINE G02BYF(M,NY,NX,ISZ,R,LDR,P,LDP,WK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes partial correlation matrix from correlation matrix.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BYF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDP, LDR, M, NX, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  P(LDP,NY), R(LDR,M), WK(NY*NX+NX*(NX+1)/2)
      INTEGER           ISZ(M)
C     .. Local Scalars ..
      INTEGER           I, IERROR, INFO, IX, IY, J, JX, JY, K, KX, KXY,
     *                  L, NREC
      LOGICAL           ONE, ZERO
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DPPTRF, DSYRK, DTPTRS
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (M.LT.3) THEN
         WRITE (P01REC,FMT=99999) M
      ELSE IF (NY.LT.2) THEN
         WRITE (P01REC,FMT=99998) NY
      ELSE IF (NX.LT.1) THEN
         WRITE (P01REC,FMT=99997) NX
      ELSE IF (NY+NX.GT.M) THEN
         NREC = 2
         WRITE (P01REC,FMT=99996) NY, NX, M
      ELSE IF (LDR.LT.M) THEN
         WRITE (P01REC,FMT=99995) LDR, M
      ELSE IF (LDP.LT.NY) THEN
         WRITE (P01REC,FMT=99994) LDP, NY
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
C
C        Check ISZ
C
         K = 0
         L = 0
         DO 20 I = 1, M
            IF (ISZ(I).LT.0) THEN
               K = K + 1
            ELSE IF (ISZ(I).GT.0) THEN
               L = L + 1
            END IF
   20    CONTINUE
         IF (K.NE.NY) THEN
            NREC = 2
            IERROR = 2
            WRITE (P01REC,FMT=99993) K, NY
            GO TO 140
         END IF
         IF (L.NE.NX) THEN
            NREC = 2
            IERROR = 2
            WRITE (P01REC,FMT=99992) L, NX
            GO TO 140
         END IF
C
C        Copy re-ordered R into workspace
C
         KX = 0
         KXY = NX*(NX+1)/2
         JY = 0
         JX = 0
         DO 60 J = 1, M
            IF (ISZ(J).LT.0) JY = JY + 1
            IF (ISZ(J).GT.0) JX = JX + 1
            IY = 0
            IX = 0
            DO 40 I = 1, J
               IF (ISZ(J).LT.0 .AND. ISZ(I).LT.0) THEN
                  IY = IY + 1
                  P(IY,JY) = R(I,J)
               ELSE IF (ISZ(J).LT.0 .AND. ISZ(I).GT.0) THEN
                  IX = IX + 1
                  WK(KXY+(JY-1)*NX+IX) = R(I,J)
               ELSE IF (ISZ(J).GT.0 .AND. ISZ(I).GT.0) THEN
                  IX = IX + 1
                  KX = KX + 1
                  WK(KX) = R(I,J)
               ELSE IF (ISZ(J).GT.0 .AND. ISZ(I).LT.0) THEN
                  IY = IY + 1
                  WK(KXY+(IY-1)*NX+JX) = R(I,J)
               END IF
   40       CONTINUE
   60    CONTINUE
         KXY = KXY + 1
         CALL DPPTRF('U',NX,WK,INFO)
         IF (INFO.GT.0) THEN
            IERROR = 3
            WRITE (P01REC,FMT=99991)
            GO TO 140
         END IF
         CALL DTPTRS('U','T','N',NX,NY,WK,WK(KXY),NX,INFO)
         IF (INFO.GT.0) THEN
            IERROR = 3
            WRITE (P01REC,FMT=99990)
            GO TO 140
         END IF
         CALL DSYRK('U','T',NY,NX,-1.0D0,WK(KXY),NX,1.0D0,P,LDP)
C
C        Calculate partial correlations in upper triangle
C        and copy partial covariances in lower triangle
C
         ZERO = .FALSE.
         ONE = .FALSE.
         DO 80 I = 1, NY
            IF (P(I,I).LE.0.0D0) THEN
               ZERO = .TRUE.
               WK(I) = 1.0D0
               P(I,I) = 0.0D0
            ELSE
               WK(I) = 1.0D0/SQRT(P(I,I))
            END IF
   80    CONTINUE
         DO 120 J = 1, NY
            DO 100 I = 1, J - 1
               P(J,I) = P(I,J)
               P(I,J) = P(I,J)*WK(I)*WK(J)
               IF (P(I,J).GT.1.0D0) THEN
                  ONE = .TRUE.
                  P(I,J) = 1.0D0
               END IF
  100       CONTINUE
  120    CONTINUE
         IF (ONE .AND. ZERO) THEN
            IERROR = 4
            NREC = 2
            WRITE (P01REC,FMT=99989)
         ELSE IF (ZERO) THEN
            IERROR = 4
            WRITE (P01REC,FMT=99988)
         ELSE IF (ONE) THEN
            IERROR = 4
            WRITE (P01REC,FMT=99987)
         END IF
      END IF
  140 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, M .le. 2 : M = ',I16)
99998 FORMAT (' ** On entry, NY .le. 1 : NY = ',I16)
99997 FORMAT (' ** On entry, NX .lt. 1 : NX = ',I16)
99996 FORMAT (' ** On entry, NY + NX .gt. M :',/'              NY = ',
     *       I16,'  NX = ',I16,' M = ',I16)
99995 FORMAT (' ** On entry, LDR .lt. M : LDR = ',I16,'  M = ',I16)
99994 FORMAT (' ** On entry, LDP .lt. NY : LDP = ',I16,'  NY = ',I16)
99993 FORMAT (' ** On entry, there are ',I16,' elements of ISZ.lt.0',
     *       /'              rather than NY: NY = ',I16)
99992 FORMAT (' ** On entry, there are ',I16,' elements of ISZ.gt.0',
     *       /'              rather than NX: NX = ',I16)
99991 FORMAT (' ** On entry, the correlation matrix of the independent',
     *       ' variables is singular')
99990 FORMAT (' ** The square root of the correlation matrix of the in',
     *       'dependent variables is singular')
99989 FORMAT (' ** A diagonal element of the partial covariance matrix',
     *       ' is zero',/' and an element of the partial correlation ',
     *       'matrix is greater than one')
99988 FORMAT (' ** A diagonal element of the partial covariance matrix',
     *       ' is zero')
99987 FORMAT (' ** An element of the partial correlation matrix is gre',
     *       'ater than one')
      END
