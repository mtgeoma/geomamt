      SUBROUTINE G08AGY(A,R,N,KI,XF)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PURPOSE
C     G08AGY RANKS A VECTOR OF VALUES AND CALCULATES
C     CORRECTION FACTOR DUE TO TIES
C
C     USAGE
C     CALL G08AGY(A,R,N,KI,XF)
C
C     DESCRIPTION OF PARAMETERS
C     A - INPUT VECTOR OF N VALUES
C     R - OUTPUT VECTOR OF LENGTH N. SMALLEST VALUE IS RANKED 1,
C         LARGEST IS RANKED N. TIES ARE ASSIGNED AVERAGE OF TIED
C         RANKS
C     N - NUMBER OF VALUES
C     KI - INPUT CODE FOR CALCULATION OF CORRECTION FACTOR
C          1   SOLVE XF=SUMME((NTIE**3-NTIE)/12.0)
C          2   SOLVE XF=SUMME(NTIE*(NTIE-1.0)/2.0)
C          3   SOLVE XF=SUMME(NTIE*(NTIE-1)*(2*NTIE+5)
C          4   SOLVE XF=SUMME(NTIE*(NTIE-1)*(NTIE-2))
C          WHERE NTIE IS THE NUMBER OF OBSERVATIONS TIED
C          FOR A GIVEN RANK
C     XF - CORRECTION FACTOR (OUTPUT)
C
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C     NONE
C
C     METHOD
C     VECTOR IS SEARCHED FOR SUCCESSIVELY LARGER ELEMENTS. IF TIES
C     OCCUR, THEY ARE LOCATED AND THEIR RANK VALUE COMPUTED.
C     FOR EXAMPLE, IF 2 VALUES ARE TIED FOR SIXTH RANK, THEY ARE
C     ASSIGNED A RANK OF 6.5 (=(6+7)/2)
C     AND CORRECTION FACTOR 1 OR 2 IS SUMMED.
C
C     -------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XF
      INTEGER           KI, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), R(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  P, X
      INTEGER           I, J, NTIE, NXLT
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE
C     .. Executable Statements ..
      DO 20 I = 1, N
         R(I) = 0.0D0
   20 CONTINUE
      XF = 0.0D0
C
C     Find rank of data
C
      DO 80 I = 1, N
C
C        Test whether data point is already ranked
C
         IF (R(I).LE.0.0D0) THEN
C
C           Data point to be ranked
C
            NXLT = 0
            NTIE = 0
            X = A(I)
            DO 40 J = 1, N
               IF (A(J).LT.X) THEN
C
C                 Count number of data points which are smaller
C
                  NXLT = NXLT + 1
C
               ELSE IF (A(J).EQ.X) THEN
C
C                 Count number of data points which are equal.
C                 Mark these by setting their ranks to -1.
C
                  NTIE = NTIE + 1
                  R(J) = -1.0D0
               END IF
   40       CONTINUE
C
C           Test for tie
C
            IF (NTIE.LE.1) THEN
C
C              Store rank of untied data points
C
               R(I) = DBLE(NXLT) + 1.0D0
C
C              Store rank of tied data points
C
            ELSE IF (NTIE.GT.1) THEN
               IF (MOD(NTIE,2).EQ.0) THEN
                  P = DBLE(NXLT) + DBLE(NTIE/2) + 0.5D0
               ELSE
                  P = DBLE(NXLT) + DBLE((NTIE+1)/2)
               END IF
               DO 60 J = I, N
                  IF (R(J).EQ.-1.0D0) R(J) = P
   60          CONTINUE
               IF (KI.EQ.1) THEN
                  XF = XF + DBLE(NTIE**3-NTIE)/12.0D0
               ELSE IF (KI.EQ.2) THEN
                  XF = XF + DBLE(NTIE*(NTIE-1)/2)
               ELSE IF (KI.EQ.3) THEN
                  XF = XF + DBLE(NTIE*(NTIE-1)*(2*NTIE+5))
               ELSE IF (KI.EQ.4) THEN
                  XF = XF + DBLE(NTIE*(NTIE-1)*(NTIE-2))
               END IF
            END IF
         END IF
   80 CONTINUE
      RETURN
      END
