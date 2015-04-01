      SUBROUTINE G04EAF(TYPE,N,LEVELS,IFACT,X,LDX,V,REP,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes dummy variables for general linear model.
C

C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G04EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDX, LEVELS, N
      CHARACTER         TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  REP(LEVELS), V(*), X(LDX,*)
      INTEGER           IFACT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, IERROR, IFAULT, NREC
      LOGICAL           FIRST, MEAN
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, G04EAX, G04EAY, G04EAZ, M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     .. Executable Statements ..
      IERROR = 1
      NREC = 1
      IF (LEVELS.LT.2) THEN
         WRITE (P01REC,FMT=99999) LEVELS
      ELSE IF (N.LT.LEVELS) THEN
         WRITE (P01REC,FMT=99998) N, LEVELS
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC,FMT=99997) LDX, N
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         IF (TYPE.EQ.'P' .OR. TYPE.EQ.'p') THEN
C
C        Check levels of IFACT
C
            CALL DCOPY(LEVELS,V,1,REP,1)
            IFAULT = 0
            CALL M01CAF(REP,1,LEVELS,'A',IFAULT)
            TEMP = REP(1)
            DO 20 I = 2, LEVELS
               IF (REP(I).LE.TEMP) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99992)
                  GO TO 160
               END IF
               TEMP = REP(I)
   20       CONTINUE
         END IF
         DO 40 I = 1, LEVELS
            REP(I) = 0.0D0
   40    CONTINUE
         DO 60 I = 1, N
            IF (IFACT(I).LE.0 .OR. IFACT(I).GT.LEVELS) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99995) I
               GO TO 160
            ELSE
               REP(IFACT(I)) = REP(IFACT(I)) + 1.0D0
            END IF
   60    CONTINUE
         DO 80 I = 1, LEVELS
            IF (REP(I).EQ.0) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99994)
               GO TO 160
            END IF
   80    CONTINUE
         IF (TYPE.EQ.'P' .OR. TYPE.EQ.'p') THEN
            DO 100 I = 1, N
               X(I,1) = V(IFACT(I))
  100       CONTINUE
            CALL G04EAZ(N,X,1,X,X,X,X,IFAULT)
            IF (IFAULT.NE.0) THEN
               IERROR = 3
               I = 1
               WRITE (P01REC,FMT=99993) I
               GO TO 160
            END IF
            CALL G04EAZ(N,X,2,X(1,2),X,X,X,IFAULT)
            IF (IFAULT.NE.0) THEN
               IERROR = 3
               I = 2
               WRITE (P01REC,FMT=99993) I
               GO TO 160
            END IF
            DO 120 I = 1, LEVELS - 1
               CALL G04EAZ(N,X,I,X(1,I),X(1,I-1),X(1,I-2),X,IFAULT)
               IF (IFAULT.NE.0) THEN
                  IERROR = 3
                  WRITE (P01REC,FMT=99993) I
                  GO TO 160
               END IF
  120       CONTINUE
         ELSE IF (TYPE.EQ.'H' .OR. TYPE.EQ.'h') THEN
            CALL G04EAY(N,IFACT,LEVELS,REP,X,LDX)
C
C        Restore REP
C
            TEMP = REP(1)
            DO 140 I = 2, LEVELS
               REP(I) = NINT(TEMP/REP(I))
               TEMP = TEMP + REP(I)
  140       CONTINUE
         ELSE IF (TYPE.EQ.'F' .OR. TYPE.EQ.'f') THEN
            MEAN = .FALSE.
            FIRST = .TRUE.
            CALL G04EAX(N,IFACT,LEVELS,MEAN,FIRST,X,LDX)
         ELSE IF (TYPE.EQ.'L' .OR. TYPE.EQ.'l') THEN
            MEAN = .FALSE.
            FIRST = .FALSE.
            CALL G04EAX(N,IFACT,LEVELS,MEAN,FIRST,X,LDX)
         ELSE IF (TYPE.EQ.'C' .OR. TYPE.EQ.'c') THEN
            MEAN = .TRUE.
            FIRST = .FALSE.
            CALL G04EAX(N,IFACT,LEVELS,MEAN,FIRST,X,LDX)
         ELSE
            IERROR = 1
            WRITE (P01REC,FMT=99996) TYPE
         END IF
      END IF
  160 CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, LEVELS .lt. 2: LEVELS = ',I16)
99998 FORMAT (' ** On entry, N .lt. LEVELS: N = ',I16,' LEVELS = ',I16)
99997 FORMAT (' ** On entry, LDX .lt. N: LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, TYPE is not valid: TYPE = ',A1)
99995 FORMAT (' ** On entry, the ',I16,'th value of IFACT is not valid')
99994 FORMAT (' ** On entry, not all levels are present in IFACT')
99993 FORMAT (' ** The ',I16,'th polynomial has all elements zero')
99992 FORMAT (' ** On entry, not all values of V are distinct')
      END
