      SUBROUTINE G13BEM(N4,MSN,MIS,MRN,NPE,MORD,MTYP,MSER,NXS)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEM DERIVES THE ORDER IN WHICH
C     THE BETA ARRAY MUST BE PROCESSED TO GIVE THE
C     PARAMETERS IN PARA ORDER. IT ALSO OUTPUTS A
C     DESCRIPTION OF EACH PARAMETER
C
C     .. Scalar Arguments ..
      INTEGER           N4, NPE, NXS
C     .. Array Arguments ..
      INTEGER           MIS(N4), MORD(NPE), MRN(N4), MSER(NPE), MSN(N4),
     *                  MTYP(NPE)
C     .. Local Scalars ..
      INTEGER           I, J, KSUB, KTYP
C     .. Executable Statements ..
      KSUB = 0
C
C     DERIVE ORDER AND DESCRIPTION FOR ARIMA PARAMETERS
C
      DO 20 I = 2, N4
         IF (MIS(I).NE.0) GO TO 20
         IF (MRN(I).GE.6) GO TO 20
         KSUB = KSUB + 1
         MTYP(KSUB) = MRN(I) - 1
         MORD(KSUB) = I - 1
         MSER(KSUB) = 0
   20 CONTINUE
      IF (NXS.LE.0) GO TO 140
C
C     PROCESS EACH INPUT SERIES IN TURN
C
      DO 120 J = 1, NXS
         DO 100 I = 2, N4
            IF (MIS(I).NE.J) GO TO 100
C
C           PROCESS SIMPLE OMEGA
C
            IF (MRN(I).NE.3) GO TO 40
            KTYP = 5
            GO TO 80
C
C           PROCESS T.F. OMEGA
C
   40       IF (MRN(I).NE.8) GO TO 60
            KTYP = 6
            GO TO 80
C
C           PROCESS T.F. DELTA
C
   60       IF (MRN(I).NE.1) GO TO 100
            KTYP = 7
C
C           DERIVE ORDER AND DESCRIPTION
C
   80       KSUB = KSUB + 1
            MTYP(KSUB) = KTYP
            MORD(KSUB) = I - 1
            MSER(KSUB) = J
  100    CONTINUE
  120 CONTINUE
C
C     DERIVE ORDER AND DESCRIPTION FOR CONSTANT
C
  140 DO 160 I = 2, N4
         IF (MRN(I).NE.11) GO TO 160
         KSUB = KSUB + 1
         MTYP(KSUB) = 8
         MORD(KSUB) = I - 1
         MSER(KSUB) = 0
         GO TO 180
  160 CONTINUE
  180 RETURN
      END
