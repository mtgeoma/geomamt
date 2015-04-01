      SUBROUTINE G03CCF(METHOD,ROTATE,NVAR,NFAC,FL,LDFL,PSI,E,R,LDR,FS,
     *                  LDFS,WK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes factor score coefficients in FS from factor
C     unrotated loadings in FL.
C     Either Bartlett's method (METHOD='B') or the Regression
C     method (METHOD='R') is used.
C     Optionally rotations can be supplied in R
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03CCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDFL, LDFS, LDR, NFAC, NVAR
      CHARACTER         METHOD, ROTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  E(NVAR), FL(LDFL,NFAC), FS(LDFS,NFAC),
     *                  PSI(NVAR), R(LDR,*), WK(NVAR)
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST, SCALE
      INTEGER           I, IERROR, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           P01ABF
      EXTERNAL          DDOT, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FCF, F06FDF, DCOPY
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
C
C     Check inputs
C
      IF (NFAC.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NFAC
      ELSE IF (NVAR.LT.NFAC) THEN
         WRITE (P01REC(1),FMT=99998) NVAR, NFAC
      ELSE IF (LDFL.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99997) LDFL, NVAR
      ELSE IF (LDFS.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99996) LDFS, NVAR
      ELSE IF ((ROTATE.EQ.'R' .OR. ROTATE.EQ.'r') .AND. (LDR.LT.NFAC))
     *         THEN
         WRITE (P01REC(1),FMT=99995) LDR, NFAC
      ELSE IF (METHOD.NE.'B' .AND. METHOD.NE.'b' .AND. METHOD.NE.
     *         'R' .AND. METHOD.NE.'r') THEN
         WRITE (P01REC(1),FMT=99994) METHOD
      ELSE IF (ROTATE.NE.'R' .AND. ROTATE.NE.'r' .AND. ROTATE.NE.
     *         'U' .AND. ROTATE.NE.'u') THEN
         WRITE (P01REC(1),FMT=99993) ROTATE
      ELSE
         IERROR = 0
C
C        compute (lamda)*inv(gamma)
C
         IF (METHOD.EQ.'R' .OR. METHOD.EQ.'r') THEN
            CONST = 0.0D0
         ELSE
            CONST = 1.0D0
         END IF
         DO 20 I = 1, NFAC
            IF (E(I).LE.1.0D0) THEN
               IERROR = 2
               WRITE (P01REC(1),FMT=99992) I
               GO TO 120
            END IF
            SCALE = 1.0D0/(E(I)-CONST)
            CALL F06FDF(NVAR,SCALE,FL(1,I),1,FS(1,I),1)
   20    CONTINUE
C
C        compute inv(psi)*(lamda)*inv(gamma)
C
         DO 40 I = 1, NVAR
            IF (PSI(I).LE.0.0D0) THEN
               IERROR = 2
               WRITE (P01REC(1),FMT=99991) I
               GO TO 120
            END IF
            WK(I) = 1.0D0/PSI(I)
   40    CONTINUE
         DO 60 I = 1, NFAC
            CALL F06FCF(NVAR,WK,1,FS(1,I),1)
   60    CONTINUE
C
C        rotate if required
C
         IF (ROTATE.EQ.'R' .OR. ROTATE.EQ.'r') THEN
            DO 100 I = 1, NVAR
               CALL DCOPY(NFAC,FS(I,1),LDFS,WK,1)
               DO 80 J = 1, NFAC
                  FS(I,J) = DDOT(NFAC,WK,1,R(1,J),1)
   80          CONTINUE
  100       CONTINUE
         END IF
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, NFAC.lt.1 : NFAC = ',I16)
99998 FORMAT (' ** On entry, NVAR.lt.NFAC: NVAR = ',I16,' NFAC = ',I16)
99997 FORMAT (' ** On entry, LDFL.lt.NVAR: LDFL = ',I16,' NVAR = ',I16)
99996 FORMAT (' ** On entry, LDFS.lt.NVAR: LDFS = ',I16,' NVAR = ',I16)
99995 FORMAT (' ** On entry, LDR.lt.NFAC : LDR = ',I16,' NFAC = ',I16)
99994 FORMAT (' ** On entry, METHOD is not a valid character: METHOD = '
     *       ,A1)
99993 FORMAT (' ** On entry, ROTATE is not a valid character: ROTATE = '
     *       ,A1)
99992 FORMAT (' ** On entry, the ',I16,' th value of E.le.1.0')
99991 FORMAT (' ** On entry, the ',I16,' th value of PSI.le.0.0')
      END
