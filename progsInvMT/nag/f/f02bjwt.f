      SUBROUTINE F02BJW(N,A,IA,B,IB,MATZ,Z,IZ)
C     MARK 6 RELEASE  NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS SUBROUTINE IS THE FIRST STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM AND THE OTHER
C     TO UPPER TRIANGULAR FORM USING ORTHOGONAL TRANSFORMATIONS.
C     IT IS USUALLY FOLLOWED BY F02BJX, F02BJY AND, POSSIBLY,
C     F02BJZ.
C
C     ON INPUT
C
C     N IS THE ORDER OF THE MATRICES.
C
C     A CONTAINS A REAL GENERAL MATRIX.
C
C     IA MUST SPECIFY THE FIRST DIMENSION OF THE ARRAY A AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IA.GE.N
C
C     B CONTAINS A REAL GENERAL MATRIX.
C
C     IB MUST SPECIFY THE FIRST DIMENSION OF THE ARRAY B AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IB.GE.N
C
C     MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND
C     TRANSFORMATIONS
C     ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C     EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     IZ MUST SPECIFY THE FIRST DIMENSION OF THE ARRAY Z AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IZ.GE.N
C
C     ON OUTPUT
C
C     A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C     BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO.
C
C     B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C     BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO.
C
C     Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS IF
C     MATZ HAS BEEN SET TO .TRUE.  OTHERWISE, Z IS NOT REFERENCED.
C
C     IMPLEMENTED FROM EISPACK BY
C     W. PHILLIPS OXFORD UNIVERSITY COMPUTING SERVICE.
C
C     INITIALIZE Z
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IZ, N
      LOGICAL           MATZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, R, RHO, S, T, U1, U2, V1, V2, ZERO
      INTEGER           I, J, K, L, L1, LB, NK1, NM1, NM2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      IF ( .NOT. MATZ) GO TO 60
C
      DO 40 I = 1, N
C
         DO 20 J = 1, N
            Z(I,J) = ZERO
   20    CONTINUE
C
         Z(I,I) = ONE
   40 CONTINUE
C     REDUCE B TO UPPER TRIANGULAR FORM
   60 IF (N.LE.1) GO TO 420
      NM1 = N - 1
C
      DO 260 L = 1, NM1
         L1 = L + 1
         S = ZERO
C
         DO 80 I = L1, N
            S = S + ABS(B(I,L))
   80    CONTINUE
C
         IF (S.EQ.ZERO) GO TO 260
         S = S + ABS(B(L,L))
         R = ZERO
C
         DO 100 I = L, N
            B(I,L) = B(I,L)/S
            R = R + B(I,L)**2
  100    CONTINUE
C
         R = SQRT(R)
         IF (B(L,L).LT.ZERO) R = -R
         B(L,L) = B(L,L) + R
         RHO = R*B(L,L)
C
         DO 160 J = L1, N
            T = ZERO
C
            DO 120 I = L, N
               T = T + B(I,L)*B(I,J)
  120       CONTINUE
C
            T = -T/RHO
C
            DO 140 I = L, N
               B(I,J) = B(I,J) + T*B(I,L)
  140       CONTINUE
C
  160    CONTINUE
C
         DO 220 J = 1, N
            T = ZERO
C
            DO 180 I = L, N
               T = T + B(I,L)*A(I,J)
  180       CONTINUE
C
            T = -T/RHO
C
            DO 200 I = L, N
               A(I,J) = A(I,J) + T*B(I,L)
  200       CONTINUE
C
  220    CONTINUE
C
         B(L,L) = -S*R
C
         DO 240 I = L1, N
            B(I,L) = ZERO
  240    CONTINUE
C
  260 CONTINUE
C     REDUCE A TO UPPER HESSENBERG FORM, WHILE
C     KEEPING B TRIANGULAR
      IF (N.EQ.2) GO TO 420
      NM2 = N - 2
C
      DO 400 K = 1, NM2
         NK1 = NM1 - K
C        FOR L=N-1 STEP -1 UNTIL K+1 DO --
         DO 380 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C           ZERO A(L+1,K)
            S = ABS(A(L,K)) + ABS(A(L1,K))
            IF (S.EQ.ZERO) GO TO 380
            U1 = A(L,K)/S
            U2 = A(L1,K)/S
            R = SQRT(U1*U1+U2*U2)
            IF (U1.LT.ZERO) R = -R
            V1 = -(U1+R)/R
            V2 = -U2/R
            U2 = V2/V1
C
            DO 280 J = K, N
               T = A(L,J) + U2*A(L1,J)
               A(L,J) = A(L,J) + T*V1
               A(L1,J) = A(L1,J) + T*V2
  280       CONTINUE
C
            A(L1,K) = ZERO
C
            DO 300 J = L, N
               T = B(L,J) + U2*B(L1,J)
               B(L,J) = B(L,J) + T*V1
               B(L1,J) = B(L1,J) + T*V2
  300       CONTINUE
C           ZERO B(L+1,L)
            S = ABS(B(L1,L1)) + ABS(B(L1,L))
            IF (S.EQ.ZERO) GO TO 380
            U1 = B(L1,L1)/S
            U2 = B(L1,L)/S
            R = SQRT(U1*U1+U2*U2)
            IF (U1.LT.ZERO) R = -R
            V1 = -(U1+R)/R
            V2 = -U2/R
            U2 = V2/V1
C
            DO 320 I = 1, L1
               T = B(I,L1) + U2*B(I,L)
               B(I,L1) = B(I,L1) + T*V1
               B(I,L) = B(I,L) + T*V2
  320       CONTINUE
C
            B(L1,L) = ZERO
C
            DO 340 I = 1, N
               T = A(I,L1) + U2*A(I,L)
               A(I,L1) = A(I,L1) + T*V1
               A(I,L) = A(I,L) + T*V2
  340       CONTINUE
C
            IF ( .NOT. MATZ) GO TO 380
C
            DO 360 I = 1, N
               T = Z(I,L1) + U2*Z(I,L)
               Z(I,L1) = Z(I,L1) + T*V1
               Z(I,L) = Z(I,L) + T*V2
  360       CONTINUE
C
  380    CONTINUE
C
  400 CONTINUE
C
  420 RETURN
      END
