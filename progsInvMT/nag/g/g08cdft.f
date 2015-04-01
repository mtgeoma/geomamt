      SUBROUTINE G08CDF(N1,X,N2,Y,NTYPE,D,Z,P,SX,SY,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08CDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D, P, Z
      INTEGER           IFAIL, N1, N2, NTYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  SX(N1), SY(N2), X(N1), Y(N2)
C     .. Local Scalars ..
      DOUBLE PRECISION  DN, DP, DT, EN1, EN2, FN1, FN2
      INTEGER           I, IERROR, IF2, J1, J2, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  G08CDZ
      INTEGER           P01ABF
      EXTERNAL          G08CDZ, P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IF (N1.LT.1 .OR. N2.LT.1) THEN
         NREC = 2
         IERROR = 1
         WRITE (P01REC,FMT=99999) N1, N2
      ELSE IF (NTYPE.NE.1 .AND. NTYPE.NE.2 .AND. NTYPE.NE.3) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) NTYPE
      ELSE
         IERROR = 0
         DO 20 I = 1, N1
            SX(I) = X(I)
   20    CONTINUE
         DO 40 I = 1, N2
            SY(I) = Y(I)
   40    CONTINUE
         IF2 = 0
         CALL M01CAF(SX,1,N1,'A',IF2)
         CALL M01CAF(SY,1,N2,'A',IF2)
         EN1 = DBLE(N1)
         EN2 = DBLE(N2)
         J1 = 1
         J2 = 1
         FN1 = 0.0D0
         FN2 = 0.0D0
         D = 0.0D0
         DP = 0.0D0
         DN = 0.0D0
   60    CONTINUE
         IF (J1.LE.N1 .AND. J2.LE.N2) THEN
            IF (SX(J1).LT.SY(J2)) THEN
               FN1 = DBLE(J1)/EN1
               DT = FN1 - FN2
               IF (DT.GT.DP) DP = DT
               J1 = J1 + 1
            ELSE IF (SX(J1).GT.SY(J2)) THEN
               FN2 = DBLE(J2)/EN2
               DT = FN2 - FN1
               IF (DT.GT.DN) DN = DT
               J2 = J2 + 1
            ELSE
C
C              Take care of case of tied sample points
C
               FN1 = DBLE(J1)/EN1
               J1 = J1 + 1
               FN2 = DBLE(J2)/EN2
               J2 = J2 + 1
   80          CONTINUE
               IF (J1.LE.N1) THEN
                  IF (SX(J1).EQ.SX(J1-1)) THEN
                     FN1 = DBLE(J1)/EN1
                     J1 = J1 + 1
                     GO TO 80
                  END IF
               END IF
  100          CONTINUE
               IF (J2.LE.N2) THEN
                  IF (SY(J2).EQ.SY(J2-1)) THEN
                     FN2 = DBLE(J2)/EN2
                     J2 = J2 + 1
                     GO TO 100
                  END IF
               END IF
               DT = FN1 - FN2
               IF (DT.GT.DP) DP = DT
               DT = FN2 - FN1
               IF (DT.GT.DN) DN = DT
            END IF
            GO TO 60
         END IF
         DP = MAX(0.0D0,DP)
         DN = MAX(0.0D0,DN)
         D = MAX(DP,DN)
         IF2 = 1
         IF (NTYPE.EQ.1) THEN
            P = G08CDZ(N1,N2,D,IF2)
         ELSE IF (NTYPE.EQ.2) THEN
            D = DP
            P = 0.5D0*G08CDZ(N1,N2,D,IF2)
         ELSE IF (NTYPE.EQ.3) THEN
            D = DN
            P = 0.5D0*G08CDZ(N1,N2,D,IF2)
         END IF
         IF (IF2.NE.0) THEN
            IERROR = 3
            NREC = 1
            WRITE (P01REC,FMT=99997)
         END IF
         Z = SQRT(DBLE(N1+N2)/DBLE(N1*N2))*D
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry either N1.lt.1 or N2.lt.1',/'    N1 = ',
     *       I16,' and N2 = ',I16)
99998 FORMAT (1X,'** On entry NTYPE is not equal to 1, 2 or 3 : NTYPE ',
     *       '= ',I16)
99997 FORMAT (1X,'** The Kolmogorov approximation failed to converge.')
      END
