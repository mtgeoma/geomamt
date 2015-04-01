      SUBROUTINE G08CCF(N,X,CDF,NTYPE,D,Z,P,SX,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G08CCF performs the Kolmogorov-Smirnov one sample test. The
C     test statistic D is computed and the tail probability is
C     returned via P. This routine requires the user to provide a
C     distribuiton function through CDF.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08CCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D, P, Z
      INTEGER           IFAIL, N, NTYPE
C     .. Array Arguments ..
      DOUBLE PRECISION  SX(N), X(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  CDF
      EXTERNAL          CDF
C     .. Local Scalars ..
      DOUBLE PRECISION  DN, DP, DT, EN, FF, FFP, FN, FO
      INTEGER           I, IERROR, IF2, J, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  G08CBZ
      INTEGER           P01ABF
      EXTERNAL          G08CBZ, P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
C     .. Executable Statements ..
C
      NREC = 1
      IERROR = 0
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (NTYPE.NE.1 .AND. NTYPE.NE.2 .AND. NTYPE.NE.3) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) NTYPE
      ELSE
         DO 20 I = 1, N
            SX(I) = X(I)
   20    CONTINUE
         IF2 = 0
         CALL M01CAF(SX,1,N,'A',IF2)
         EN = DBLE(N)
         FO = 0.0D0
         FFP = 0.0D0
         D = 0.0D0
         DP = 0.0D0
         DN = 0.0D0
         DO 40 J = 1, N
            FN = DBLE(J)/EN
            FF = CDF(SX(J))
            IF (FF.LT.0.0D0 .OR. FF.GT.1.0D0) THEN
               IERROR = 3
               NREC = 2
               WRITE (P01REC,FMT=99997) SX(J), FF
               GO TO 60
            ELSE IF (FF.LT.FFP) THEN
               IERROR = 4
               NREC = 3
               WRITE (P01REC,FMT=99996) SX(J-1), FFP, SX(J), FF
               GO TO 60
            END IF
            FFP = FF
            DT = FN - FF
            IF (DT.GT.DP) DP = DT
            DT = FF - FO
            IF (DT.GT.DN) DN = DT
            FO = FN
   40    CONTINUE
         DP = MAX(0.0D0,DP)
         DN = MAX(0.0D0,DN)
         IF2 = 1
         IF (NTYPE.EQ.1) THEN
            D = MAX(DP,DN)
            P = 2.0D0*G08CBZ(N,D)
            IF (P.GT.1.0D0) P = 1.0D0
         ELSE IF (NTYPE.EQ.2) THEN
            D = DP
            P = G08CBZ(N,D)
         ELSE IF (NTYPE.EQ.3) THEN
            D = DN
            P = G08CBZ(N,D)
         END IF
         Z = SQRT(EN)*D
      END IF
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry N.lt.1: N = ',I16)
99998 FORMAT (1X,'** On entry NTYPE is not equal to 1, 2 or 3: NTYPE = '
     *       ,I16)
99997 FORMAT (1X,'** The supplied theoretical distribution returns a v',
     *       'alue less than zero or ',/4X,'greater than one : at x = ',
     *       D13.5,' ,  F(x) = ',D13.5)
99996 FORMAT (1X,'** The supplied theoretical distribution is not non-',
     *       'decreasing.',/4X,'At x = ',D13.5,' ,  F(x) = ',D13.5,' a',
     *       'nd',/4X,'at x = ',D13.5,' ,  F(x) = ',D13.5)
      END
