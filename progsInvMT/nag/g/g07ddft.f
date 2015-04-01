      SUBROUTINE G07DDF(N,X,ALPHA,TMEAN,WMEAN,TVAR,WVAR,K,SX,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Given data X(1), X(2), ..., X(N) in non-decreasing order,
C     this routine finds the alpha-trimmed mean TMEAN, and the
C     estimated variance of TMEAN, TVAR.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, HALF
      CHARACTER*6       SRNAME
      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,SRNAME='G07DDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, TMEAN, TVAR, WMEAN, WVAR
      INTEGER           IFAIL, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  SX(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, SUM2, Z, ZN, ZNT
      INTEGER           I, IERROR, N1, N2, NREC
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01CAF
C     .. Executable Statements ..
      NREC = 1
      IF (N.LE.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (ALPHA.LT.ZERO .OR. ALPHA.GE.HALF) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) ALPHA
      ELSE
         IERROR = 0
         DO 20 I = 1, N
            SX(I) = X(I)
   20    CONTINUE
         CALL M01CAF(SX,1,N,'A',IFAIL)
         ZN = N
C                       K  is number trimmed from each end
C                       N1 is number of lowest observation not trimmed
C                       N2 is number of highest observation not trimmed
         K = ALPHA*ZN
         N1 = K + 1
         N2 = N - K
         ZNT = K
C                       Trimmed Mean (TMEAN) &  Winsorized Mean (WMEAN)
         SUM = ZERO
         DO 40 I = N1, N2
            SUM = SUM + SX(I)
   40    CONTINUE
         Z = N - 2*K
         TMEAN = SUM/Z
         SUM = SUM + ZNT*(SX(N1)+SX(N2))
         WMEAN = SUM/ZN
C                       Winsorized sum of squares about TMEAN (TVAR)
C                                            &    about WMEAN (WVAR)
         SUM = ZERO
         SUM2 = ZERO
         DO 60 I = N1, N2
            SUM = SUM + (TMEAN-SX(I))*(TMEAN-SX(I))
            SUM2 = SUM2 + (WMEAN-SX(I))*(WMEAN-SX(I))
   60    CONTINUE
         IF (N1.EQ.N2) THEN
            TVAR = ZERO
            WVAR = ZERO
         ELSE
            IF (K.NE.0) THEN
               SUM = SUM + ZNT*((SX(N1)-TMEAN)*(SX(N1)-TMEAN)) +
     *               ZNT*((SX(N2)-TMEAN)*(SX(N2)-TMEAN))
               SUM2 = SUM2 + ZNT*((SX(N1)-WMEAN)*(SX(N1)-WMEAN)) +
     *                ZNT*((SX(N2)-WMEAN)*(SX(N2)-WMEAN))
            END IF
            TVAR = SUM/(ZN*ZN)
            WVAR = SUM2/(ZN*ZN)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry N.le.1 : N = ',I16)
99998 FORMAT (1X,'** On entry ALPHA.lt.0 or ALPHA.ge.0.5 : ALPHA = ',
     *       D13.5)
      END
