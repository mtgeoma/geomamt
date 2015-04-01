      DOUBLE PRECISION FUNCTION G01ERF(T,VK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     --------------------------------------------------------
C     G01ERF returns the left tail area of the Von Mises
C     distribution, equal to the incomplete modified Bessel
C     function, of the first kind and zero-th order.
C
C     Parameters :
C
C         T = Angle in radians, treated as deviation from zero
C             (mean) reduced modulo 2PI to the range (-PI,+PI).
C
C        VK = concentration parameter, Kappa, VK.GE.0.0.
C
C     IFAIL = Error parameter
C
C             IFAIL = 1 : On input value of VK was .lt. 0.0.
C
C     --------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01ERF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 T, VK
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 A1, A2, A3, A4, C, C1, CN, P, PI,
     *                                 R, S, SN, TPI, U, V, Y, Z
      INTEGER                          IERROR, IFAULT, IP, N, NREC
      LOGICAL                          ASYMP
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S15ABF, X01AAF, X02AJF
      INTEGER                          P01ABF
      EXTERNAL                         S15ABF, X01AAF, X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        COS, DBLE, MAX, MIN, MOD, SIN,
     *                                 SQRT
C     .. Data statements ..
      DATA                             A1, A2, A3, A4/28.0D0, 0.5D0,
     *                                 100.0D0, 5.0D0/
C     .. Executable Statements ..
C
C
      G01ERF = 0.0D0
      PI = X01AAF(0.0D0)
      TPI = 2*PI
      NREC = 1
      IF (VK.LT.0.0D0) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) VK
      ELSE
         IERROR = 0
C
C        Check most quick exits : allows THETA = -PI.
C
         Y = MOD(T,TPI)
         IF (Y.EQ.-PI) THEN
            G01ERF = 0.0D0
         ELSE IF (Y.EQ.0.0D0) THEN
            G01ERF = 0.5D0
         ELSE IF (Y.EQ.PI) THEN
            G01ERF = 1.0D0
         ELSE
            IF (VK.GT.6.5D0) THEN
               IF (X02AJF().LT.1.0D-10) THEN
                  IF (VK.GT.50.0D0) THEN
                     ASYMP = .TRUE.
                     C1 = 50.1D0
                  ELSE
                     ASYMP = .FALSE.
                  END IF
               ELSE
                  ASYMP = .TRUE.
                  C1 = 62.0D0
               END IF
            ELSE
               ASYMP = .FALSE.
            END IF
            Z = VK
C
C           Convert angle Y, with shift allowed for U,
C           modulo 2PI to range (-PI,+PI).
C
            U = MOD(Y+PI,TPI)
            IF (U.LT.0.0D0) U = U + TPI
            Y = U - PI
            IF (ASYMP) THEN
C
C              For large VK compute the normal approximation
C              and left tail.
C
               C = 24.0D0*Z
               V = C - C1
               R = SQRT((54.0D0/(347.0D0/V+26.0D0-C)-6.0D0+C)/6.0D0)
               Z = SIN(Y*0.5D0)*R
               S = Z*Z
               V = V - S + 3.0D0
               Y = (C-S-S-16.0D0)/3.0D0
               Y = ((S+1.75D0)*S+83.5D0)/V - Y
               IFAULT = 1
               G01ERF = S15ABF((Z-S/(Y*Y)*Z),IFAULT)
            ELSE
               V = 0.0D0
               IF (Z.GT.0.0D0) THEN
C
C                 For small VK sum IP terms by backwards recursion.
C
                  IP = Z*A2 - A3/(Z+A4) + A1
                  P = DBLE(IP)
                  S = SIN(Y)
                  C = COS(Y)
                  Y = P*Y
                  SN = SIN(Y)
                  CN = COS(Y)
                  R = 0.0D0
                  Z = 2.0D0/Z
                  DO 20 N = 2, IP
                     P = P - 1.0D0
                     Y = SN
                     SN = SN*C - CN*S
                     CN = CN*C + Y*S
                     R = 1.0D0/(P*Z+R)
                     V = (SN/P+V)*R
   20             CONTINUE
               END IF
C
C              Calculate probability with shift on angle.
C
               G01ERF = (U*0.5D0+V)/PI
            END IF
C
C           Check for rounding error due to machine accuracy.
C
            G01ERF = MIN(MAX(0.0D0,G01ERF),1.0D0)
C
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, VK.lt.0.0: VK = ',D13.5)
      END
