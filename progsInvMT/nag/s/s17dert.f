      SUBROUTINE S17DER(Z,FNU,N,CY,TOL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-761 (DEC 1989).
C
C     Original name: CRATI
C
C     S17DER COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD
C     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
C     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
C     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
C     BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
C     BY D. J. SOOKNE.
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  FNU, TOL
      INTEGER           N
C     .. Array Arguments ..
      COMPLEX*16        CY(N)
C     .. Local Scalars ..
      COMPLEX*16        CDFNU, CONE, CZERO, P1, P2, PT, RZ, T1
      DOUBLE PRECISION  AK, AMAGZ, AP1, AP2, ARG, AZ, DFNU, FDNU, FLAM,
     *                  FNUP, RAP1, RHO, TEST, TEST1
      INTEGER           I, ID, IDNU, INU, ITIME, K, KK, MAGZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, INT, MAX, MIN, DBLE, SQRT
C     .. Data statements ..
      DATA              CZERO, CONE/(0.0D0,0.0D0), (1.0D0,0.0D0)/
C     .. Executable Statements ..
C
      AZ = ABS(Z)
      INU = INT(FNU)
      IDNU = INU + N - 1
      FDNU = IDNU
      MAGZ = INT(AZ)
      AMAGZ = MAGZ + 1
      FNUP = MAX(AMAGZ,FDNU)
      ID = IDNU - MAGZ - 1
      ITIME = 1
      K = 1
      RZ = (CONE+CONE)/Z
      T1 = DCMPLX(FNUP,0.0D0)*RZ
      P2 = -T1
      P1 = CONE
      T1 = T1 + RZ
      IF (ID.GT.0) ID = 0
      AP2 = ABS(P2)
      AP1 = ABS(P1)
C     ------------------------------------------------------------------
C     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
C     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
C     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
C     PREMATURELY.
C     ------------------------------------------------------------------
      ARG = (AP2+AP2)/(AP1*TOL)
      TEST1 = SQRT(ARG)
      TEST = TEST1
      RAP1 = 1.0D0/AP1
      P1 = P1*DCMPLX(RAP1,0.0D0)
      P2 = P2*DCMPLX(RAP1,0.0D0)
      AP2 = AP2*RAP1
   20 CONTINUE
      K = K + 1
      AP1 = AP2
      PT = P2
      P2 = P1 - T1*P2
      P1 = PT
      T1 = T1 + RZ
      AP2 = ABS(P2)
      IF (AP1.LE.TEST) THEN
         GO TO 20
      ELSE IF (ITIME.NE.2) THEN
         AK = ABS(T1)*0.5D0
         FLAM = AK + SQRT(AK*AK-1.0D0)
         RHO = MIN(AP2/AP1,FLAM)
         TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0D0))
         ITIME = 2
         GO TO 20
      END IF
      KK = K + 1 - ID
      AK = KK
      DFNU = FNU + N - 1
      CDFNU = DCMPLX(DFNU,0.0D0)
      T1 = DCMPLX(AK,0.0D0)
      P1 = DCMPLX(1.0D0/AP2,0.0D0)
      P2 = CZERO
      DO 40 I = 1, KK
         PT = P1
         P1 = RZ*(CDFNU+T1)*P1 + P2
         P2 = PT
         T1 = T1 - CONE
   40 CONTINUE
      IF (DBLE(P1).EQ.0.0D0 .AND. DIMAG(P1).EQ.0.0D0) P1 = DCMPLX(TOL,
     *    TOL)
      CY(N) = P2/P1
      IF (N.NE.1) THEN
         K = N - 1
         AK = K
         T1 = DCMPLX(AK,0.0D0)
         CDFNU = DCMPLX(FNU,0.0D0)*RZ
         DO 60 I = 2, N
            PT = CDFNU + T1*RZ + CY(K+1)
            IF (DBLE(PT).EQ.0.0D0 .AND. DIMAG(PT).EQ.0.0D0)
     *          PT = DCMPLX(TOL,TOL)
            CY(K) = CONE/PT
            T1 = T1 - CONE
            K = K - 1
   60    CONTINUE
      END IF
      RETURN
      END
