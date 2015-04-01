      SUBROUTINE F01MAX(IN1,IN2,NZ,IP,N,A)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     IP( I )CONTAINS THE NUMBER OF NON-ZEROS IN ROW I OF THE INPUT
C     MATRIX. INITIALIZE IP( I )TO POINT JUST BEYOND WHERE THE LAST
C     ELEMENT OF THE ROW WILL BE STORED.
C
C     REORDER USING IN - PLACE SORT.
C
C     .. Scalar Arguments ..
      INTEGER           N, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NZ)
      INTEGER           IN1(NZ), IN2(NZ), IP(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2
      INTEGER           I, IC1, IC2, IDUMMY, IR1, IR2, KI
C     .. Executable Statements ..
      DO 80 I = 1, NZ
C
C        SAVE CURRENT ENTRY.
C
         IR1 = IN1(I)
C
C        IF IR1.LT.0 THE ELEMENT IS IN PLACE ALREADY.
C
         IF (IR1.LT.0) GO TO 80
         IC1 = IN2(I)
         A1 = A(I)
C
C        DETERMINE CORRECT POSITION.
C
         KI = IP(IR1) - 1
         DO 20 IDUMMY = 1, NZ
            IF (I.EQ.KI) GO TO 40
C
C           SAVE CONTENTS OF THAT POSITION.
C
            IR2 = IN1(KI)
            IC2 = IN2(KI)
C
C           STORE CURRENT ENTRY
C
            IN1(KI) = -IR1
            IN2(KI) = IC1
            IP(IR1) = KI
            IR1 = IR2
            IC1 = IC2
C
C           MAKE CORRESPONDING CHANGES FOR REALS IF REQUIRED.
C
            A2 = A(KI)
            A(KI) = A1
            A1 = A2
            KI = IP(IR1) - 1
   20    CONTINUE
   40    IF (IDUMMY.EQ.1) GO TO 60
C
C        IF CURRENT ENTRY IS IN PLACE IT IS STORED HERE.
C
         A(KI) = A1
         IN2(KI) = IC1
         IN1(KI) = -IR1
   60    IP(IR1) = I
   80 CONTINUE
      RETURN
C
C     END OF F01MAX. ( MA31E. )
C
      END
