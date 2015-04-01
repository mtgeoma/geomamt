      SUBROUTINE H02BBX(KM1,K,NXANOD,FREANO,ACTNOD,FATHER,RITSON,LFTSON,
     *                  FREDEL,NODTST)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     SUBROUTINE REMOVES A NODE FROM THE ACTIVE LIST AND THEN PARES
C     DELTA TREE AS APPROPRIATE.
C
C     .. Scalar Arguments ..
      INTEGER           FREANO, FREDEL, K, KM1, NODTST
C     .. Array Arguments ..
      INTEGER           ACTNOD(*), FATHER(*), LFTSON(*), NXANOD(0:*),
     *                  RITSON(*)
C     .. Local Scalars ..
      INTEGER           L, M
C     .. Executable Statements ..
C
      NXANOD(KM1) = NXANOD(K)
      NXANOD(K) = FREANO
      FREANO = K
      NODTST = NODTST - 1
      L = ACTNOD(K)
      K = NXANOD(KM1)
   20 CONTINUE
      IF (L.LE.1) GO TO 40
      M = FATHER(L)
      RITSON(L) = FREDEL
      FREDEL = L
      IF (LFTSON(M).EQ.L) THEN
         LFTSON(M) = 0
         IF (RITSON(M).NE.0) GO TO 40
      ELSE
         RITSON(M) = 0
         IF (LFTSON(M).NE.0) GO TO 40
      END IF
      L = M
      GO TO 20
   40 CONTINUE
      RETURN
      END
