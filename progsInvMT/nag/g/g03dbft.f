      SUBROUTINE G03DBF(EQUAL,MODE,NVAR,NG,GMEAN,LDG,GC,NOBS,M,ISX,X,
     *                  LDX,D,LDD,WK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C      computes generalized distance from point X to group centroids
C
C      for use after G03DAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03DBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDD, LDG, LDX, M, NG, NOBS, NVAR
      CHARACTER         EQUAL, MODE
C     .. Array Arguments ..
      DOUBLE PRECISION  D(LDD,NG), GC((NG+1)*NVAR*(NVAR+1)/2),
     *                  GMEAN(LDG,NVAR), WK(2*NVAR), X(LDX,*)
      INTEGER           ISX(*)
C     .. Local Scalars ..
      INTEGER           I, IERROR, IND, J, K, L, NC, NR, NREC
      LOGICAL           LOWER, MEAN
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DTPSV
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 1
      IF (NVAR.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NVAR
      ELSE IF (NG.LT.2) THEN
         WRITE (P01REC(1),FMT=99998) NG
      ELSE IF (LDG.LT.NG) THEN
         WRITE (P01REC(1),FMT=99995) LDG, NG
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 180
      IF (MODE.EQ.'S' .OR. MODE.EQ.'s') THEN
         IERROR = 1
         IF (NOBS.LT.1) THEN
            WRITE (P01REC(1),FMT=99997) NOBS
         ELSE IF (M.LT.NVAR) THEN
            WRITE (P01REC(1),FMT=99996) M, NVAR
         ELSE IF (LDX.LT.NOBS) THEN
            WRITE (P01REC(1),FMT=99994) LDX, NOBS
         ELSE IF (LDD.LT.NOBS) THEN
            WRITE (P01REC(1),FMT=99993) LDD, NOBS
         ELSE
            IERROR = 0
         END IF
         IF (IERROR.EQ.1) GO TO 180
         NR = NOBS
         LOWER = .FALSE.
         MEAN = .FALSE.
         K = 0
         DO 20 I = 1, M
            IF (ISX(I).GT.0) K = K + 1
   20    CONTINUE
         IF (K.NE.NVAR) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99989) K, NVAR
            GO TO 180
         END IF
      ELSE IF (MODE.EQ.'M' .OR. MODE.EQ.'m') THEN
         IF (LDD.LT.NG) THEN
            WRITE (P01REC(1),FMT=99987) LDD, NG
            IERROR = 1
            GO TO 180
         END IF
         NR = NG
         MEAN = .TRUE.
         IF (EQUAL.EQ.'E' .OR. EQUAL.EQ.'e') THEN
            LOWER = .TRUE.
         ELSE
            LOWER = .FALSE.
         END IF
      ELSE
         IERROR = 1
         WRITE (P01REC(1),FMT=99988) MODE
         GO TO 180
      END IF
C
C     equal covariance case
C
      IF (EQUAL.EQ.'E' .OR. EQUAL.EQ.'e') THEN
         IND = 0
C
C        check diagonals of R matrix
C
         L = 0
         DO 40 I = 1, NVAR
            L = L + I
            IF (GC(L).EQ.0.0D0) THEN
               IERROR = 2
               WRITE (P01REC(1),FMT=99991) I
               GO TO 180
            END IF
   40    CONTINUE
C
C        unequal covariance case
C
      ELSE IF (EQUAL.EQ.'U' .OR. EQUAL.EQ.'u') THEN
C
C        check diagonals of R matrices
C
         IND = NVAR*(NVAR+1)/2
         L = IND
         DO 80 J = 1, NG
            DO 60 I = 1, NVAR
               L = L + I
               IF (GC(L).EQ.0.0D0) THEN
                  NREC = 2
                  IERROR = 2
                  WRITE (P01REC,FMT=99990) I, J
                  GO TO 180
               END IF
   60       CONTINUE
   80    CONTINUE
      ELSE
         IERROR = 1
         WRITE (P01REC(1),FMT=99992) EQUAL
         GO TO 180
      END IF
C
C     compute distances
C
      NC = NG
      DO 160 I = 1, NR
         IF (MEAN) THEN
            CALL DCOPY(NG,GMEAN(I,1),LDG,WK(NVAR+1),1)
         ELSE
            L = NVAR
            DO 100 J = 1, M
               IF (ISX(J).GT.0) THEN
                  L = L + 1
                  WK(L) = X(I,J)
               END IF
  100       CONTINUE
         END IF
         IF (LOWER) NC = I - 1
         DO 140 J = 1, NC
            IF (MEAN .AND. I.EQ.J) GO TO 140
            DO 120 K = 1, NVAR
               WK(K) = WK(NVAR+K) - GMEAN(J,K)
  120       CONTINUE
            CALL DTPSV('U','T','N',NVAR,GC(J*IND+1),WK,1)
            D(I,J) = DDOT(NVAR,WK,1,WK,1)
  140    CONTINUE
  160 CONTINUE
  180 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, NVAR.lt.1 : NVAR = ',I16)
99998 FORMAT (' ** On entry, NG.lt.2 : NG = ',I16)
99997 FORMAT (' ** On entry, NOBS.lt.1 : NOBS = ',I16)
99996 FORMAT (' ** On entry, M.lt.NVAR : M = ',I16,' NVAR = ',I16)
99995 FORMAT (' ** On entry, LDG.lt.NG : LDG = ',I16,' NG = ',I16)
99994 FORMAT (' ** On entry, LDX.lt.NOBS : LDX = ',I16,' NOBS = ',I16)
99993 FORMAT (' ** On entry, LDD.lt.NOBS : LDD = ',I16,' NOBS = ',I16)
99992 FORMAT (' ** On entry, EQUAL is not a valid character: EQUAL = ',
     *       A1)
99991 FORMAT (' ** On entry, the ',I16,' th diagonal element of R.eq.0')
99990 FORMAT (' ** On entry, the ',I16,' th diagonal element of the ',
     *       /'                 ',I16,' th R matrix.eq.0')
99989 FORMAT (' ** On entry, ',I16,' values of ISX.gt.0, not NVAR = ',
     *       I16)
99988 FORMAT (' ** On entry MODE is not a valid character: MODE = ',A1)
99987 FORMAT (' ** On entry LDD .lt. NG : LDD = ',I16,' NG = ',I16)
      END
