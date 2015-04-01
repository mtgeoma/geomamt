      SUBROUTINE G07BEZ(N,NE,T,GAMMA,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     COMPUTES INITIAL VALUES FOR WEIBULL DISTRIBUTION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GAMMA
      INTEGER           IFAULT, N, NE
C     .. Array Arguments ..
      DOUBLE PRECISION  T(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  LMS, LMT, LPS, LPT, P, TT, UMS, UMT, UPS, UPT
      INTEGER           I, IEX, IML, IMU, IN, ND, NJ, NJ1, NP
      LOGICAL           EVEN
C     .. External Subroutines ..
      EXTERNAL          G07BEY, M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, MOD, DBLE
C     .. Executable Statements ..
C
C     SORT DATA
C
      EVEN = .FALSE.
      IFAULT = 0
      IF (NE.EQ.N) THEN
         CALL M01CAF(T,1,N,'A',IFAULT)
      ELSE
         CALL G07BEY(T,1,N,'A',IFAULT)
      END IF
      IEX = NE
      TT = T(1)
      DO 20 I = 2, N
         IF (T(I).GT.0.0D0) THEN
            IF (T(I).EQ.TT) THEN
               IEX = IEX - 1
            ELSE
               TT = T(I)
            END IF
         END IF
   20 CONTINUE
C
      IF (IEX.LT.2 .OR. (IEX.EQ.2 .AND. T(N).GT.0.0D0)) THEN
         IFAULT = 3
C
C        FIND QUANTILES
C
      ELSE
         IN = IEX/2
         IML = (IN+1)/2
         IMU = IN + IML
         IF (IEX.GT.3) THEN
            IF (MOD(IN,2).EQ.0) EVEN = .TRUE.
         END IF
C
         P = 1.0D0
         NP = 0
         NJ1 = N
         I = 0
   40    CONTINUE
C
C        LOOP THROUGHT CALCULATING SURVIVAL FUNTION
C
         NJ = NJ1
         ND = 0
   60    CONTINUE
         I = I + 1
         IF (ABS(T(I)).EQ.ABS(T(I+1))) THEN
            NJ1 = NJ1 - 1
            IF (T(I).GE.0) ND = ND + 1
            IF (I.LT.N) GO TO 60
         END IF
         NJ1 = NJ1 - 1
         IF (T(I).GE.0) ND = ND + 1
         IF (ND.NE.0) THEN
            NP = NP + 1
            P = P*DBLE(NJ-ND)/DBLE(NJ)
C
C           STORE VALUES FOR QUANTILE POINTS
C
            IF (NP.EQ.IML) THEN
               LMS = P
               LMT = T(I)
            ELSE IF (NP.EQ.IMU) THEN
               UMS = P
               UMT = T(I)
            END IF
            IF (EVEN) THEN
               IF (NP.EQ.IML+1) THEN
                  LMS = 0.5D0*(LMS+P)
                  LMT = 0.5D0*(LMT+T(I))
               ELSE IF (NP.EQ.IMU+1) THEN
                  UMS = 0.5D0*(UMS+P)
                  UMT = 0.5D0*(UMT+T(I))
               END IF
            END IF
         END IF
         IF (NP.LT.IMU .OR. (NP.LT.IMU+1 .AND. (EVEN))) GO TO 40
         IF (UMS.LE.0.0D0) THEN
            IFAULT = 3
C
C           COMPUTES ESTIMATES OF SLOPE
C
         ELSE
            LPS = LOG(-LOG(LMS))
            LPT = LOG(LMT)
            UPS = LOG(-LOG(UMS))
            UPT = LOG(UMT)
            GAMMA = (UPS-LPS)/(UPT-LPT)
         END IF
      END IF
      RETURN
      END
