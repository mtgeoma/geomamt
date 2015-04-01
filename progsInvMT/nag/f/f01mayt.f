      SUBROUTINE F01MAY(A,IRN,IA,N,IK,IP,ROW,MA31J,MA31K)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Scalar Arguments ..
      INTEGER           IA, N
      LOGICAL           ROW
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA)
      INTEGER           IK(N), IP(N), IRN(IA), MA31J(5), MA31K(3)
C     .. Local Scalars ..
      INTEGER           IPI, J, K, KK, KL, KN, NN, NZ
C     .. Executable Statements ..
      MA31J(3) = MA31J(3) + 1
      DO 20 J = 1, N
C
C        STORE FIRST ELEMENT OF ENTRY IN IK( J ). THEN OVERWRITE IT
C        BY  -J
C
         NZ = IK(J)
         NN = MA31K(3)
         IF ( .NOT. ROW) NN = MA31K(2)
         IF (NZ.LE.0 .OR. IP(J).LT.NN) GO TO 20
         K = IP(J)
         IK(J) = IRN(K)
         IRN(K) = -J
   20 CONTINUE
C
C     KN IS POSITION OF NEXT ENTRY IN COMPRESSED FILE.
C
      KN = IA + 1
      IPI = IA + 1
      KL = IA - MA31K(2) + 1
      IF (ROW) KL = IA - MA31K(3) + 1
C
C     LOOP THROUGH OLD FILE SKIPPING ZERO ELEMENTS AND MOVING GENUINE
C     ELEMENTS FORWARD. THE ENTRY NUMBER BECOMES KNOWN ONLY WHEN
C     ITS END IS DETECTED BY PRESENCE OF A NEGATIVE INTEGER.
C
      DO 60 KK = 1, KL
         K = IA + 1 - KK
         IF (IRN(K).EQ.0) GO TO 60
         KN = KN - 1
         IF (ROW) A(KN) = A(K)
         IF (IRN(K).GE.0) GO TO 40
C
C        END OF ENTRY. RESTORE IRN( K ), SET POINTERS TO START OF
C        ENTRY AND STORE CURRENT KN IN IPI READY FOR USE WHEN NEXT
C        LAST ENTRY IS DETECTED.
C
         J = -IRN(K)
         IRN(K) = IK(J)
         IP(J) = KN
         IK(J) = IPI - KN
         IPI = KN
   40    IRN(KN) = IRN(K)
   60 CONTINUE
      IF (ROW) GO TO 80
      MA31K(2) = KN
      GO TO 100
   80 MA31K(3) = KN
  100 RETURN
C
C     END OF F01MAY. ( MA31D. )
C
      END
