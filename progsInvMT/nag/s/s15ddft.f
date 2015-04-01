      COMPLEX*16  FUNCTION S15DDF(Z,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Evaluates the function w(z) = exp(-z*z)*erfc(-iz), for
C     complex z. Derived from ACM Algorithm 363. See W.Gautschi,
C     Efficient Computation of the Complex Error Function, SIAM
C     J. Numer. Anal. 7, 1970, pp 187-198.
C     .. Parameters ..
      DOUBLE PRECISION            SQTPI2, HALF, ONE, TWO, ZERO
      PARAMETER                   (SQTPI2=1.128379167095512573896D0,
     *                            HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,
     *                            ZERO=0.0D0)
      CHARACTER*6                 SRNAME
      PARAMETER                   (SRNAME='S15DDF')
C     .. Scalar Arguments ..
      COMPLEX*16                  Z
      INTEGER                     IFAIL
C     .. Local Scalars ..
      COMPLEX*16                  EXPMZ2, IZ, R, S, W, ZZ
      DOUBLE PRECISION            EPS, EXMZ2I, EXMZ2R, H, H0, HSAFE,
     *                            LAMBDA, RS, SAFE, WFUNI, WFUNR, X, X0,
     *                            Y, Y0
      INTEGER                     CAPN, IER, N, N0, N1, NREC, NU, NU0,
     *                            NU1, QUAD, TIFAIL
      LOGICAL                     B, FIRST
C     .. Local Arrays ..
      CHARACTER*80                REC(2)
C     .. External Functions ..
      COMPLEX*16                  S01EAF
      DOUBLE PRECISION            X02AJF, X02AMF
      INTEGER                     P01ABF
      EXTERNAL                    S01EAF, X02AJF, X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                   ABS, DIMAG, DCMPLX, DCONJG, EXP, LOG,
     *                            MAX, DBLE, SIGN, SQRT
C     .. Save statement ..
      SAVE                        X0, Y0, H0, N0, N1, NU0, NU1, SAFE,
     *                            HSAFE, FIRST
C     .. Data statements ..
      DATA                        FIRST/.TRUE./
C     .. Executable Statements ..
C
      IER = 0
      IF (FIRST) THEN
C        Initialise parameters depending on machine precision.
         FIRST = .FALSE.
         SAFE = ONE/X02AMF()
         HSAFE = HALF*SAFE
         EPS = X02AJF()
         IF (EPS.LT.0.5D-16) THEN
C           Use parameters for 18 digit accuracy.
            X0 = 10.15D0
            Y0 = 13.41D0
            H0 = 2.0D0
            N0 = 6
            N1 = 46
            NU0 = 9
            NU1 = 48
         ELSE IF (EPS.LT.0.5D-14) THEN
C           Use parameters for 16 digit accuracy.
            X0 = 8.72D0
            Y0 = 10.06D0
            H0 = 2.0D0
            N0 = 6
            N1 = 40
            NU0 = 9
            NU1 = 39
         ELSE IF (EPS.LT.0.5D-12) THEN
C           Use parameters for 14 digit accuracy.
            X0 = 7.44D0
            Y0 = 7.58D0
            H0 = 1.9D0
            N0 = 6
            N1 = 34
            NU0 = 9
            NU1 = 31
         ELSE IF (EPS.LT.0.5D-10) THEN
C           Use parameters for 12 digit accuracy.
            X0 = 6.31D0
            Y0 = 5.73D0
            H0 = 1.7D0
            N0 = 6
            N1 = 28
            NU0 = 9
            NU1 = 28
         ELSE IF (EPS.LT.0.5D-8) THEN
C           Use parameters for 10 digit accuracy.
            X0 = 5.33D0
            Y0 = 4.29D0
            H0 = 1.6D0
            N0 = 6
            N1 = 23
            NU0 = 9
            NU1 = 21
         ELSE
C           Use parameters for 8 digit accuracy.
            X0 = 4.48D0
            Y0 = 3.17D0
            H0 = 1.6D0
            N0 = 6
            N1 = 18
            NU0 = 9
            NU1 = 15
         END IF
      END IF
C
      X = DBLE(Z)
      Y = DIMAG(Z)
C     Determine the quadrant containing Z.
      IF (X.GE.ZERO) THEN
         IF (Y.GE.ZERO) THEN
            QUAD = 1
         ELSE
            QUAD = 4
         END IF
      ELSE
         IF (Y.GE.ZERO) THEN
            QUAD = 2
         ELSE
            QUAD = 3
         END IF
      END IF
C
      X = ABS(X)
      Y = ABS(Y)
C     Use Gautschi's algorithm to compute w(z) for z in first quadrant.
      IF (X.LT.X0 .AND. Y.LT.Y0) THEN
         RS = (ONE-Y/Y0)*SQRT(ONE-(X/X0)**2)
         H = H0*RS
         CAPN = N0 + N1*RS
         NU = NU0 + NU1*RS
         LAMBDA = (TWO*H)**CAPN
         B = LAMBDA .EQ. ZERO
      ELSE
         H = ZERO
         CAPN = 0
         NU = 8
         B = .TRUE.
      END IF
      R = ZERO
      S = ZERO
      IZ = DCMPLX(-Y,X)
      DO 20 N = NU, 0, -1
         R = HALF/(H-IZ+(N+1)*R)
         IF (H.GT.ZERO .AND. N.LE.CAPN) THEN
            S = R*(LAMBDA+S)
            LAMBDA = LAMBDA/(TWO*H)
         END IF
   20 CONTINUE
      IF (B) THEN
         W = SQTPI2*R
      ELSE
         W = SQTPI2*S
      END IF
      IF (Y.EQ.ZERO) THEN
         IF (ABS(X).GT.SQRT(LOG(SAFE))) THEN
            W = DCMPLX(ZERO,DIMAG(W))
         ELSE
            W = DCMPLX(EXP(-X*X),DIMAG(W))
         END IF
      END IF
C
C     Modify the result depending on the quadrant of z.
      IF (QUAD.EQ.1) THEN
         S15DDF = W
      ELSE IF (QUAD.EQ.2) THEN
         S15DDF = DCONJG(W)
      ELSE
C        We have to compute 2*exp(-z*z) - w(z).
         IF (MAX(X,Y).GE.SQRT(HSAFE)) THEN
C           Overflow is likely in computation of Z*Z.
            IF (X.GT.Y) THEN
               IF ((Y/X)*Y-X.LT.-LOG(SAFE)/X) THEN
                  EXPMZ2 = ZERO
               ELSE
                  IER = 5
               END IF
            ELSE
               IER = 5
            END IF
            TIFAIL = 0
         ELSE
            ZZ = DCMPLX(X,Y)
            TIFAIL = 1
            EXPMZ2 = S01EAF(-ZZ*ZZ,TIFAIL)
            IF (TIFAIL.GE.4) IER = TIFAIL
C           If TIFAIL is 4 or 5, we have lost precision in exp(-ZZ**2).
C           Otherwise we can ignore TIFAIL, because the surrogate result
C           returned by S01EAF will carry through below.
         END IF
         IF (IER.EQ.5) THEN
            S15DDF = ZERO
         ELSE
C           Handle real and imaginary parts separately to avoid
C           possible overflow.
            EXMZ2R = DBLE(EXPMZ2)
            EXMZ2I = DIMAG(EXPMZ2)
            IF (ABS(EXMZ2R).LE.HSAFE) THEN
               WFUNR = TWO*EXMZ2R - DBLE(W)
            ELSE
               IER = 1
               WFUNR = SIGN(SAFE,EXMZ2R)
            END IF
            IF (ABS(EXMZ2I).LE.HSAFE) THEN
               WFUNI = TWO*EXMZ2I - DIMAG(W)
            ELSE
               IER = IER + 2
               WFUNI = SIGN(SAFE,EXMZ2I)
            END IF
            IF (QUAD.EQ.3) THEN
               S15DDF = DCMPLX(WFUNR,WFUNI)
            ELSE
               S15DDF = DCMPLX(WFUNR,-WFUNI)
            END IF
         END IF
         IF (IER.EQ.1) THEN
            NREC = 2
            WRITE (REC,FMT=99999) Z
         ELSE IF (IER.EQ.2) THEN
            NREC = 2
            WRITE (REC,FMT=99998) Z
         ELSE IF (IER.EQ.3) THEN
            NREC = 2
            WRITE (REC,FMT=99997) Z
         ELSE IF (IER.EQ.4) THEN
            NREC = 2
            WRITE (REC,FMT=99996) Z
         ELSE IF (IER.EQ.5) THEN
            NREC = 2
            WRITE (REC,FMT=99995) Z
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** Real part of result overflows when entered with',
     *       /4X,'argument Z = (',1P,D13.5,',',D13.5,')')
99998 FORMAT (1X,'** Imaginary part of result overflows when entered w',
     *       'ith',/4X,'argument Z = (',1P,D13.5,',',D13.5,')')
99997 FORMAT (1X,'** Both real and imaginary parts of result overflow ',
     *       'when entered with',/4X,'argument Z = (',1P,D13.5,',',
     *       D13.5,')')
99996 FORMAT (1X,'** Result has less than half precision when entered ',
     *       'with',/4X,'argument Z = (',1P,D13.5,',',D13.5,')')
99995 FORMAT (1X,'** Result has no precision when entered with',/4X,
     *       'argument Z = (',1P,D13.5,',',D13.5,')')
      END
