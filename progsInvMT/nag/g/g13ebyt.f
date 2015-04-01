      SUBROUTINE G13EBY(N,P,A,LDA,C,LDC,U,LDU,WITHU,UPPER)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Based on SLICOT routine TB01BD
C
C     Reduces the pair (A,C) to lower or upper observer Hessenberg
C     form using (and optionally accumulating) unitary state-space
C     transformations.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDA, LDC, LDU, N, P
      LOGICAL           UPPER, WITHU
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(LDC,*), U(LDU,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DZ, TL
      INTEGER           II, J, JJ, N1, NJ, P1, PAR1, PAR2, PAR3, PAR4,
     *                  PAR5, PAR6
C     .. External Subroutines ..
      EXTERNAL          F06FSF, F06FUF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      TL = ZERO
      P1 = P + 1
      N1 = N - 1
C
C        Perform transformations involving both C and A.
C
      DO 120 J = 1, MIN(P,N1)
         NJ = N - J
         IF (UPPER) THEN
            PAR1 = P - J + 1
            PAR2 = NJ + 1
            PAR3 = 1
            PAR4 = P - J
            PAR5 = NJ
         ELSE
            PAR1 = J
            PAR2 = J
            PAR3 = J + 1
            PAR4 = P
            PAR5 = N
         END IF
         CALL F06FSF(NJ,C(PAR1,PAR2),C(PAR1,PAR3),LDC,TL,DZ)
C
C           Update A.
C
         DO 20 JJ = 1, N
            CALL F06FUF(NJ,C(PAR1,PAR3),LDC,DZ,A(PAR2,JJ),A(PAR3,JJ),1)
   20    CONTINUE
         DO 40 II = 1, N
            CALL F06FUF(NJ,C(PAR1,PAR3),LDC,DZ,A(II,PAR2),A(II,PAR3),
     *                  LDA)
   40    CONTINUE
C
         IF (WITHU) THEN
C
C              Update U.
C
            DO 60 JJ = 1, N
               CALL F06FUF(NJ,C(PAR1,PAR3),LDC,DZ,U(PAR2,JJ),U(PAR3,JJ),
     *                     1)
   60       CONTINUE
         END IF
         IF (J.NE.P) THEN
C
C              Update C.
C
            DO 80 JJ = PAR3, PAR4
               CALL F06FUF(NJ,C(PAR1,PAR3),LDC,DZ,C(JJ,PAR2),C(JJ,PAR3),
     *                     LDC)
   80       CONTINUE
         END IF
         DO 100 II = PAR3, PAR5
            C(PAR1,II) = ZERO
  100    CONTINUE
  120 CONTINUE
      IF (P1.LE.N1) THEN
         DO 220 J = P1, N1
C
C              Perform next transformations only involving A.
C
            NJ = N - J
            IF (UPPER) THEN
               PAR1 = N + P1 - J
               PAR2 = NJ + 1
               PAR3 = 1
               PAR4 = NJ
               PAR5 = 1
               PAR6 = N + P - J
            ELSE
               PAR1 = J - P
               PAR2 = J
               PAR3 = J + 1
               PAR4 = N
               PAR5 = J - P + 1
               PAR6 = N
            END IF
            IF (NJ.GT.0) THEN
               CALL F06FSF(NJ,A(PAR1,PAR2),A(PAR1,PAR3),LDA,TL,DZ)
C
C                 Update A.
C
               DO 140 JJ = 1, N
                  CALL F06FUF(NJ,A(PAR1,PAR3),LDA,DZ,A(PAR2,JJ),
     *                        A(PAR3,JJ),1)
  140          CONTINUE
               DO 160 II = PAR5, PAR6
                  CALL F06FUF(NJ,A(PAR1,PAR3),LDA,DZ,A(II,PAR2),
     *                        A(II,PAR3),LDA)
  160          CONTINUE
C
               IF (WITHU) THEN
C
C                    Update U.
C
                  DO 180 JJ = 1, N
                     CALL F06FUF(NJ,A(PAR1,PAR3),LDA,DZ,U(PAR2,JJ),
     *                           U(PAR3,JJ),1)
  180             CONTINUE
               END IF
               DO 200 II = PAR3, PAR4
                  A(PAR1,II) = ZERO
  200          CONTINUE
            END IF
  220    CONTINUE
      END IF
      RETURN
      END
