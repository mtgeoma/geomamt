      SUBROUTINE G13BEN(MQAB,MQSG,JA,JB,NXSP,KRN,MSN,MSPA,MSPB,IMP,K,
     *                  KSPA,KSPB,MIS,MRN)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEN DERIVES THE ORDER IN WHICH THE A(T)
C     AND B(T) SETS RELATING TO
C     OMEGAS, DELTAS AND PRE-XS
C     WILL BE HELD IN THE BETA ARRAY FOR PROCESSING IN
C     S, G, AND H CALCULATIONS. IT GIVES THE RELATIVE
C     START POINTS.
C
C     .. Scalar Arguments ..
      INTEGER           IMP, JA, JB, K, KRN, KSPA, KSPB, NXSP
C     .. Array Arguments ..
      INTEGER           MIS(IMP), MQAB(8,NXSP), MQSG(8,NXSP), MRN(IMP),
     *                  MSN(IMP), MSPA(IMP), MSPB(IMP)
C     .. Local Scalars ..
      INTEGER           I, J, NSG
C     .. Executable Statements ..
      DO 60 J = JA, JB
         NSG = MQSG(KRN,J)
         IF (NSG.LE.0) GO TO 60
         DO 40 I = 1, NSG
            K = K + 1
C
C           MSN HOLDS THE SUBSCRIPTS OF THE RELEVANT A(T) AND B(T)
C           SETS. A NEGATIVE VALUE INDICATES A CHANGE OF SIGN WHEN
C           CALCULATIONS INVOLVE THESE SETS.
C
            MSN(K) = MQAB(KRN,J)
            MIS(K) = J
            MRN(K) = KRN
            IF (KRN.NE.8) GO TO 20
            IF (I.EQ.1) GO TO 20
            MSN(K) = -MSN(K)
C
C           MSPA AND MSPB HOLD RELATIVE START POINTS
C           FOR THE A(T) AND B(T) SETS.
C
   20       MSPA(K) = KSPA - KSPB*I
            MSPB(K) = MSPA(K)
   40    CONTINUE
   60 CONTINUE
      RETURN
      END
