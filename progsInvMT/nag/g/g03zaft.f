      SUBROUTINE G03ZAF(N,M,X,LDX,NVAR,ISX,S,E,Z,LDZ,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes standardized values (z-scores) for data matrix X.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03ZAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDX, LDZ, M, N, NVAR
C     .. Array Arguments ..
      DOUBLE PRECISION  E(M), S(M), X(LDX,M), Z(LDZ,NVAR)
      INTEGER           ISX(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, SCALE
      INTEGER           I, IERROR, J, K, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (N.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (NVAR.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) NVAR
      ELSE IF (M.LT.NVAR) THEN
         WRITE (P01REC(1),FMT=99997) M, NVAR
      ELSE IF (LDX.LT.N) THEN
         WRITE (P01REC(1),FMT=99996) LDX, N
      ELSE IF (LDZ.LT.N) THEN
         WRITE (P01REC(1),FMT=99995) LDZ, N
      ELSE
         IERROR = 0
C
C        check values of ISX
C
         K = 0
         DO 20 I = 1, M
            IF (ISX(I).NE.0) K = K + 1
   20    CONTINUE
         IF (K.NE.NVAR) THEN
            NREC = 2
            IERROR = 2
            WRITE (P01REC,FMT=99994) K, NVAR
            GO TO 80
         END IF
C
C        Compute z-scores
C
         K = 0
         DO 60 J = 1, M
            IF (ISX(J).NE.0) THEN
               IF (S(J).LE.0.0D0) THEN
                  IERROR = 3
                  WRITE (P01REC(1),FMT=99993) J
                  GO TO 80
               END IF
               SCALE = 1.0D0/S(J)
               A = E(J)
               K = K + 1
               DO 40 I = 1, N
                  Z(I,K) = (X(I,J)-A)*SCALE
   40          CONTINUE
            END IF
   60    CONTINUE
      END IF
   80 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, N.lt.1: N = ',I16)
99998 FORMAT (' ** On entry, NVAR.lt.1: NVAR = ',I16)
99997 FORMAT (' ** On entry, M.lt.NVAR: M = ',I16,' NVAR = ',I16)
99996 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,' N = ',I16)
99995 FORMAT (' ** On entry, LDZ.lt.N : LDZ = ',I16,' N = ',I16)
99994 FORMAT (' ** On entry, ',I16,' values of ISX.gt.0,',/'          ',
     *       '   rather than NVAR values: NVAR =',I16)
99993 FORMAT (' ** On entry, the ',I16,' th value of S.le.0.0')
      END
