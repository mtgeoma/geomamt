      SUBROUTINE D01JAX(IRAD2)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     D01JAX GENERATES, AT EACH CALL, A NEW PERMUTATION
C     IX(1)...IX(N) OF THE NUMBERS IX1,...,IXN, SATISFYING
C         IX1**2 + ... + IXN**2 = IRAD2    (1)
C     AND
C         IX1.GE.IX2.GE. ... .GE.IXN.GE.0. (2)
C     THE PERMUTATIONS ARE GENERATED IN INVERSE-LEXICOGRAPHIC ORDER.
C     LET K=(K(1)...K(N)) AND L=(L(1)...L(N)) BE TWO PERMUTATIONS OF
C     THE NUMBERS IX1...IXN, THEN K IS GENERATED BEFORE L IF THERE
C     IS SOME P, 1.LE.P.LT.N, SUCH THAT K(J)=L(J) FOR 1.LE.J.LT.P
C     BUT K(P).GT.L(P).
C     WHEN ALL THE PERMUTATIONS OF IX1,...,IXN ARE GENERATED, A NEW
C     COMBINATION IX1,...,IXN  SATISFYING (1) WILL BE GENERATED, AND
C     IX(1),...,IX(N) WILL BE SET TO IX1,...,IXN (BY THE ROUTINE
C     D01JAW).
C
C     MAJOR VARIABLES
C     ---------------
C
C     INDEX - SEE COMMENTS IN THE ROUTINE D01JAW.
C             WHEN INDEX=0, D01JAX IS CALLED FOR THE FIRST TIME
C             WITH THE GIVEN VALUE OF IRAD2.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IRAD2
C     .. Scalars in Common ..
      INTEGER           INDEX, INOG, IREST, N, NMIN, NPLUS
C     .. Arrays in Common ..
      INTEGER           IX(4), IX2(4)
C     .. Local Scalars ..
      INTEGER           I, IHULP, J, K, L
C     .. External Subroutines ..
      EXTERNAL          D01JAW
C     .. Common blocks ..
      COMMON            /BD01JA/N, NMIN, NPLUS, IX, IX2, INDEX, INOG,
     *                  IREST
C     .. Executable Statements ..
      IF (INDEX.EQ.0) GO TO 60
      I = NMIN
   20 IF (IX(I).GT.IX(I+1)) GO TO 80
      I = I - 1
      IF (I.NE.0) GO TO 20
C
C     ALL THE PERMUTATIONS  HAVE BEEN GENERATED. IX(1),...,IX(N) IS
C     RESET TO IX1,...,IXN AND THE ROUTINE D01JAW IS CALLED.
C
      K = 1
      L = N
   40 IHULP = IX(K)
      IX(K) = IX(L)
      IX(L) = IHULP
      K = K + 1
      L = L - 1
      IF (K.LT.L) GO TO 40
   60 CALL D01JAW(IRAD2)
      RETURN
   80 J = N
  100 IF (IX(J).LT.IX(I)) GO TO 120
      J = J - 1
      GO TO 100
  120 IHULP = IX(I)
      IX(I) = IX(J)
      IX(J) = IHULP
      K = I + 1
      L = N
  140 IF (K.GE.L) RETURN
      IHULP = IX(K)
      IX(K) = IX(L)
      IX(L) = IHULP
      K = K + 1
      L = L - 1
      GO TO 140
      END
