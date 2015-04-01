      SUBROUTINE G03EHF(ORIENT,N,DORD,DMIN,DSTEP,NSYM,C,LENC,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Prints dendogram to character array
C
C     Uses results from G03ECF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03EHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DMIN, DSTEP
      INTEGER           IFAIL, LENC, N, NSYM
      CHARACTER         ORIENT
C     .. Array Arguments ..
      DOUBLE PRECISION  DORD(N)
      CHARACTER*(*)     C(LENC)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DL, DU
      INTEGER           CLEN, I, IDINC, IDLL, IDUL, IEND, IERROR, IST, J
      LOGICAL           LINE, MAX
      CHARACTER         LINC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      DL = DMIN
      DU = DL + DSTEP
      IERROR = 0
      IF (N.LT.2) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) N
      ELSE IF (NSYM.LT.1) THEN
         IERROR = 1
         WRITE (REC,FMT=99998) NSYM
      ELSE IF (DMIN.LT.0.0D0) THEN
         IERROR = 1
         WRITE (REC,FMT=99997) DMIN
      ELSE IF (DSTEP.LE.0.0D0) THEN
         IERROR = 1
         WRITE (REC,FMT=99996) DSTEP
      ELSE IF (ORIENT.NE.'W' .AND. ORIENT.NE.'w' .AND. ORIENT.NE.
     *         'N' .AND. ORIENT.NE.'n' .AND. ORIENT.NE.'S' .AND.
     *         ORIENT.NE.'s' .AND. ORIENT.NE.'E' .AND. ORIENT.NE.'e')
     *         THEN
         IERROR = 1
         WRITE (REC,FMT=99995) ORIENT
      ELSE IF ((ORIENT.EQ.'N' .OR. ORIENT.EQ.'n' .OR. ORIENT.EQ.'S' .OR.
     *         ORIENT.EQ.'s') .AND. (LENC.LT.NSYM)) THEN
         IERROR = 1
         WRITE (REC,FMT=99994) LENC
      ELSE IF ((ORIENT.EQ.'E' .OR. ORIENT.EQ.'e' .OR. ORIENT.EQ.'W' .OR.
     *         ORIENT.EQ.'w') .AND. (LENC.LT.N)) THEN
         IERROR = 1
         WRITE (REC,FMT=99994) LENC
      ELSE
         CLEN = LEN(C(1))
         IF ((ORIENT.EQ.'N' .OR. ORIENT.EQ.'n' .OR. ORIENT.EQ.'S' .OR.
     *       ORIENT.EQ.'s') .AND. (CLEN.LT.3*N)) THEN
            IERROR = 1
            WRITE (REC,FMT=99992)
         ELSE IF ((ORIENT.EQ.'E' .OR. ORIENT.EQ.'e' .OR. ORIENT.EQ.
     *            'W' .OR. ORIENT.EQ.'w') .AND. (CLEN.LT.NSYM)) THEN
            IERROR = 1
            WRITE (REC,FMT=99992)
         END IF
         MAX = .TRUE.
         DO 20 I = 1, N - 1
            MAX = MAX .AND. (DORD(N).GE.DORD(I))
   20    CONTINUE
         IF ( .NOT. MAX) THEN
            IERROR = 2
            WRITE (REC,FMT=99993)
         END IF
      END IF
      IF (IERROR.EQ.0) THEN
         DO 40 I = 1, LENC
            C(I) = ' '
   40    CONTINUE
         IF (ORIENT.EQ.'W' .OR. ORIENT.EQ.'w' .OR. ORIENT.EQ.'E' .OR.
     *       ORIENT.EQ.'e') THEN
            IF (ORIENT.EQ.'W' .OR. ORIENT.EQ.'w') THEN
               LINC = ')'
               IDUL = NSYM
               IDLL = 1
               IDINC = 1
            ELSE
               LINC = '('
               IDUL = 1
               IDLL = NSYM
               IDINC = -1
            END IF
            DO 80 J = IDLL, IDUL, IDINC
               LINE = .FALSE.
               IF (DORD(1).GE.DU) THEN
                  C(1) (J:J) = '.'
               ELSE IF (DORD(1).LT.DL) THEN
                  C(1) (J:J) = ' '
               ELSE
                  C(1) (J:J) = ' '
                  LINE = .TRUE.
               END IF
               DO 60 I = 2, N
                  D = DORD(I)
                  IF (LINE) THEN
                     IF (D.LT.DL) THEN
                        C(I) (J:J) = LINC
                     ELSE IF (D.GE.DU) THEN
                        C(I) (J:J) = LINC
                        LINE = .FALSE.
                     ELSE
                        C(I) (J:J) = LINC
                     END IF
                  ELSE
                     IF (D.GE.DU) THEN
                        C(I) (J:J) = '.'
                     ELSE IF (D.LT.DL) THEN
                        C(I) (J:J) = ' '
                     ELSE
                        C(I) (J:J) = ' '
                        LINE = .TRUE.
                     END IF
                  END IF
   60          CONTINUE
               DL = DL + DSTEP
               DU = DU + DSTEP
   80       CONTINUE
         ELSE
            IF (ORIENT.EQ.'S' .OR. ORIENT.EQ.'s') THEN
               IDLL = NSYM
               IDUL = 1
               IDINC = -1
            ELSE
               IDLL = 1
               IDUL = NSYM
               IDINC = 1
            END IF
            DO 120 J = IDLL, IDUL, IDINC
               LINE = .FALSE.
               IF (DORD(1).GE.DU) THEN
                  C(J) (1:3) = 'I  '
               ELSE IF (DORD(1).LT.DL) THEN
                  C(J) (1:3) = '   '
               ELSE
                  C(J) (1:3) = '---'
                  LINE = .TRUE.
               END IF
               DO 100 I = 2, N - 1
                  D = DORD(I)
                  IEND = I*3
                  IST = IEND - 2
                  IF (LINE) THEN
                     IF (D.LT.DL) THEN
                        C(J) (IST:IEND) = '---'
                     ELSE IF (D.GE.DU) THEN
                        C(J) (IST:IEND) = '*  '
                        LINE = .FALSE.
                     ELSE
                        C(J) (IST:IEND) = '---'
                     END IF
                  ELSE
                     IF (D.GE.DU) THEN
                        C(J) (IST:IEND) = 'I  '
                     ELSE IF (D.LT.DL) THEN
                        C(J) (IST:IEND) = '   '
                     ELSE
                        C(J) (IST:IEND) = '---'
                        LINE = .TRUE.
                     END IF
                  END IF
  100          CONTINUE
               D = DORD(N)
               IST = IEND + 1
               IEND = IST + 2
               IF (LINE) THEN
                  IF (D.GE.DU) THEN
                     C(J) (IST:IEND) = '*  '
                  ELSE
                     C(J) (IST:IEND) = '-  '
                  END IF
               ELSE
                  IF (D.GE.DU) THEN
                     C(J) (IST:IEND) = 'I  '
                  ELSE IF (D.LT.DL) THEN
                     C(J) (IST:IEND) = '   '
                  ELSE
                     C(J) (IST:IEND) = '-  '
                  END IF
               END IF
               DL = DL + DSTEP
               DU = DU + DSTEP
  120       CONTINUE
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.2: N = ',I16)
99998 FORMAT (1X,'** On entry, NSYM.lt.1: NSYM = ',I16)
99997 FORMAT (1X,'** On entry, DMIN.lt.0.0: DMIN = ',D13.5)
99996 FORMAT (1X,'** On entry, DSTEP.le.0.0: DSTEP = ',D13.5)
99995 FORMAT (1X,'** On entry, ORIENT is not valid: ORIENT = ',A1)
99994 FORMAT (1X,'** On entry, LENC is too small: LENC = ',I16)
99993 FORMAT (1X,'** On entry, DORD(I).gt.DORD(N), I.lt.N')
99992 FORMAT (1X,'** On entry, insufficient characters can be stored i',
     *       'n each element of vector C')
      END
