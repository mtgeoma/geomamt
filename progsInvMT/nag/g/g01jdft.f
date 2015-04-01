      SUBROUTINE G01JDF(METHOD,N,RLAM,D,C,PROB,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Calculates the lower tail probability for a linear combination
C     of (central) chi-square variables.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01JDF')
      DOUBLE PRECISION  LAMTOL
      PARAMETER         (LAMTOL=0.01D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, D, PROB
      INTEGER           IFAIL, N
      CHARACTER         METHOD
C     .. Array Arguments ..
      DOUBLE PRECISION  RLAM(N), WORK(N+1)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, WI, WII, WM, XX
      INTEGER           I, I1, IERROR, IFAULT, NN, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF, G01JDV, G01JDY, G01JDZ
      INTEGER           P01ABF
      EXTERNAL          G01ECF, G01JDV, G01JDY, G01JDZ, P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
C
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (D.LT.0.0D0) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) D
      ELSE IF (METHOD.NE.'P' .AND. METHOD.NE.'I' .AND. METHOD.NE.
     *         'D' .AND. METHOD.NE.'p' .AND. METHOD.NE.'i' .AND.
     *         METHOD.NE.'d') THEN
         IERROR = 1
         WRITE (P01REC,FMT=99997) METHOD
      ELSE
         NN = 0
         DO 20 I = 1, N
            TEMP = RLAM(I) - D
            IF (TEMP.NE.0.0D0) THEN
               NN = NN + 1
               WORK(NN) = TEMP
            END IF
   20    CONTINUE
         IF (NN.LT.1) THEN
            IERROR = 2
            WRITE (P01REC,FMT=99996)
         ELSE IF (NN.NE.1) THEN
            CALL M01CAF(WORK,1,NN,'ASCENDING',IERROR)
            DO 40 I = 1, NN - 1
               I1 = I + 1
               WI = WORK(I)
               WII = WORK(I1)
               WM = MAX(ABS(WI),ABS(WII))
               TEMP = ABS(WII-WI)/WM
               IF (TEMP.LE.LAMTOL) GO TO 60
   40       CONTINUE
            IF (METHOD.EQ.'P' .OR. METHOD.EQ.'p') THEN
               PROB = G01JDZ(NN,WORK,C)
            ELSE IF (METHOD.EQ.'I' .OR. METHOD.EQ.'i') THEN
               IF (NN.GE.6) THEN
                  PROB = G01JDY(NN,WORK,C)
               ELSE
                  PROB = G01JDV(NN,WORK,C)
               END IF
            ELSE IF (NN.LE.60) THEN
               PROB = G01JDZ(NN,WORK,C)
            ELSE
               PROB = G01JDY(NN,WORK,C)
            END IF
            GO TO 80
   60       IF (METHOD.EQ.'P' .OR. METHOD.EQ.'p') THEN
               NREC = 2
               IERROR = 3
               WRITE (P01REC,FMT=99995)
            ELSE IF (NN.GE.6) THEN
               PROB = G01JDY(NN,WORK,C)
            ELSE
               PROB = G01JDV(NN,WORK,C)
            END IF
         ELSE IF (WORK(1).GT.0.0D0) THEN
            XX = C/WORK(1)
            IF (C.LE.0.0D0) THEN
               PROB = 0.0D0
            ELSE
               IFAULT = 1
               PROB = G01ECF('L',XX,1.0D0,IFAULT)
            END IF
         ELSE IF (WORK(1).LT.0.0D0) THEN
            XX = C/WORK(1)
            IF (C.LT.0.0D0) THEN
               IFAULT = 1
               PROB = G01ECF('U',XX,1.0D0,IFAULT)
            ELSE
               PROB = 1.0D0
            END IF
         END IF
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (' ** On entry, D.lt.0.0 : D = ',D13.5)
99997 FORMAT (' ** On entry, METHOD is an invalid character : METHOD = '
     *       ,A1)
99996 FORMAT (' ** On entry, all values of RLAM .eq. D')
99995 FORMAT (' ** On entry, METHOD.eq.''P'' but two successive values '
     *       ,/'    of LAMBDA-STAR were not 1% distinct. ')
      END
