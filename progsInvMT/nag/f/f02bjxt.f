      SUBROUTINE F02BJX(N,A,IA,B,IB,EPS1,MATZ,Z,IZ,ITER,IERR)
C     MARK 6 RELEASE  NAG COPYRIGHT 1977
C     MARK 8 REVISED. IER-226 (MAR 1980).
C     MARK 9 REVISED. IER-326 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-508 (AUG 1986).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     THIS SUBROUTINE IS THE SECOND STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR
C     FORM.
C     IT REDUCES THE HESSENBERG MATRIX TO QUASI-TRIANGULAR FORM
C     USING
C     ORTHOGONAL TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR
C     FORM
C     OF THE OTHER MATRIX.  IT IS USUALLY PRECEDED BY F02BJW AND
C     FOLLOWED BY F02BJY AND, POSSIBLY, F02BJZ.
C
C     ON INPUT
C
C     N IS THE ORDER OF THE MATRICES.
C
C     A CONTAINS A REAL UPPER HESSENBERG MATRIX.
C
C     IA MUST SPECIFY THE FIRST DIMENSION OF ARRAY A AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IA.GE.N
C
C     B CONTAINS A REAL UPPER TRIANGULAR MATRIX.
C
C     IB MUST SPECIFY THE FIRST DIMENSION OF ARRAY B AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IB.GE.N
C
C     EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C     EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C     ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C     ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C     POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C     IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C     POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C     BUT LESS ACCURATE RESULTS.
C
C     MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND
C     TRANSFORMATIONS
C     ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C     EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     Z CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C     TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C     BY F02BJW, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C     IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     IZ MUST SPECIFY THE FIRST DIMENSION OF ARRAY Z AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IZ.GE.N
C
C     ON OUTPUT
C
C     A HAS BEEN REDUCED TO QUASI-TRIANGULAR FORM.  THE ELEMENTS
C     BELOW THE FIRST SUBDIAGONAL ARE STILL ZERO AND NO TWO
C     CONSECUTIVE SUBDIAGONAL ELEMENTS ARE NONZERO.
C
C     B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C     HAVE BEEN ALTERED.  THE LOCATION B(N,1) IS USED TO STORE
C     EPS1 TIMES THE NORM OF B FOR LATER USE BY F02BJY AND F02BJZ.
C
C     Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C     (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE..
C
C     IERR IS SET TO
C     ZERO       FOR NORMAL RETURN,
C     J          IF MORE THAN 30*N ITERATIONS HAVE BEEN
C                PERFORMED ALTOGETHER.
C
C     IMPLEMENTED FROM EISPACK BY
C     W. PHILLIPS OXFORD UNIVERSITY COMPUTING SERVICE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BJX')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS1
      INTEGER           IA, IB, IERR, IZ, N
      LOGICAL           MATZ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), Z(IZ,N)
      INTEGER           ITER(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A11, A12, A2, A21, A22, A3, A33, A34, A43,
     *                  A44, ANI, ANORM, B11, B12, B22, B33, B34, B44,
     *                  BIT, BNI, BNORM, EP, EPSA, EPSB, HALF, ONE, R,
     *                  S, SH, T, TWO, U1, U2, U3, V1, V2, V3, XXXX,
     *                  ZERO
      INTEGER           EN, ENM2, ENORN, I, ISAVE, ISH, ITSTOT, J, K,
     *                  K1, K2, KM1, L, L1, LD, LL, LM1, LOR1, NA
      LOGICAL           NOTLAS
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/, TWO/2.0D0/,
     *                  BIT/1.1605D0/, HALF/0.5D0/
C     .. Executable Statements ..
      ISAVE = IERR
      IERR = 0
C     COMPUTE EPSA,EPSB
      ANORM = ZERO
      BNORM = ZERO
      ITSTOT = 30*N
C
      DO 40 I = 1, N
         ITER(I) = 0
         ANI = ZERO
         IF (I.NE.1) ANI = ABS(A(I,I-1))
         BNI = ZERO
C
         DO 20 J = I, N
            ANI = ANI + ABS(A(I,J))
            BNI = BNI + ABS(B(I,J))
   20    CONTINUE
C
         IF (ANI.GT.ANORM) ANORM = ANI
         IF (BNI.GT.BNORM) BNORM = BNI
   40 CONTINUE
C
      IF (ANORM.EQ.ZERO) ANORM = ONE
      IF (BNORM.EQ.ZERO) BNORM = ONE
      EP = EPS1
      IF (EP.LE.ZERO) EP = X02AJF()
      EPSA = EP*ANORM
      EPSB = EP*BNORM
C     REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
C     KEEPING B TRIANGULAR
      LOR1 = 1
      ENORN = N
      EN = N
C     BEGIN QZ STEP
   60 IF (EN.LE.2) GO TO 580
      IF ( .NOT. MATZ) ENORN = EN
      NA = EN - 1
      ENM2 = NA - 1
   80 ISH = 2
C     CHECK FOR CONVERGENCE OR REDUCIBILITY.
C     FOR L=EN STEP -1 UNTIL 1 DO --
      DO 100 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L.EQ.1) GO TO 140
         IF (ABS(A(L,LM1)).LE.EPSA) GO TO 120
  100 CONTINUE
C
  120 A(L,LM1) = ZERO
      IF (L.LT.NA) GO TO 140
C     1-BY-1 OR 2-BY-2 BLOCK ISOLATED
      EN = LM1
      GO TO 60
C     CHECK FOR SMALL TOP OF B
  140 LD = L
  160 L1 = L + 1
      B11 = B(L,L)
      IF (ABS(B11).GT.EPSB) GO TO 200
      B(L,L) = ZERO
      S = ABS(A(L,L)) + ABS(A(L1,L))
      U1 = A(L,L)/S
      U2 = A(L1,L)/S
      R = SQRT(U1*U1+U2*U2)
      IF (U1.LT.ZERO) R = -R
      V1 = -(U1+R)/R
      V2 = -U2/R
      U2 = V2/V1
C
      DO 180 J = L, ENORN
         T = A(L,J) + U2*A(L1,J)
         A(L,J) = A(L,J) + T*V1
         A(L1,J) = A(L1,J) + T*V2
         T = B(L,J) + U2*B(L1,J)
         B(L,J) = B(L,J) + T*V1
         B(L1,J) = B(L1,J) + T*V2
  180 CONTINUE
C
      IF (L.NE.1) A(L,LM1) = -A(L,LM1)
      LM1 = L
      L = L1
      GO TO 120
  200 A11 = A(L,L)/B11
      A21 = A(L1,L)/B11
      IF (ISH.EQ.1) GO TO 240
C     ITERATION STRATEGY
      IF (ITSTOT.LE.0) GO TO 560
      IF (ITER(EN).EQ.10) GO TO 280
C     DETERMINE TYPE OF SHIFT
      B22 = B(L1,L1)
      IF (ABS(B22).LT.EPSB) B22 = EPSB
      B33 = B(NA,NA)
      IF (ABS(B33).LT.EPSB) B33 = EPSB
      B44 = B(EN,EN)
      IF (ABS(B44).LT.EPSB) B44 = EPSB
      A33 = A(NA,NA)
      A34 = A(NA,EN)
      A43 = A(EN,NA)
      A44 = A(EN,EN)
      B34 = B(NA,EN)
      T = HALF*(A43*B34-A33*B44-A44*B33)
      R = T*T + (A34*A43-A33*A44)*B33*B44
      IF (R.LT.ZERO) GO TO 260
C     DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A
      ISH = 1
      R = SQRT(R)
      SH = -T + R
      S = -T - R
      IF (ABS(S-A44*B33).LT.ABS(SH-A44*B44)) SH = S
      SH = (SH/B33)/B44
C     LOOK FOR TWO CONSECUTIVE SMALL
C     SUB-DIAGONAL ELEMENTS OF A.
C     FOR L=EN-2 STEP -1 UNTIL LD DO --
      DO 220 LL = LD, ENM2
         L = ENM2 + LD - LL
         IF (L.EQ.LD) GO TO 240
         LM1 = L - 1
         L1 = L + 1
         T = A(L,L)
         IF (ABS(B(L,L)).GT.EPSB) T = T - SH*B(L,L)
         IF (ABS(A(L,LM1)).LE.ABS(T/A(L1,L))*EPSA) GO TO 160
  220 CONTINUE
C
  240 A1 = A11 - SH
      A2 = A21
      IF (L.NE.LD) A(L,LM1) = -A(L,LM1)
      GO TO 300
C     DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A
  260 A12 = A(L,L1)/B22
      A22 = A(L1,L1)/B22
      B12 = B(L,L1)/B22
      A33 = A33/B33
      A34 = A34/B44
      A43 = A43/B33
      A44 = A44/B44
      B34 = B34/B44
      A1 = ((A33-A11)*(A44-A11)-A34*A43+A43*B34*A11)/A21 + A12 - A11*B12
      A2 = (A22-A11) - A21*B12 - (A33-A11) - (A44-A11) + A43*B34
      A3 = A(L1+1,L1)/B22
      GO TO 300
C     AD HOC SHIFT
  280 A1 = ZERO
      A2 = ONE
      A3 = BIT
  300 ITER(EN) = ITER(EN) + 1
      ITSTOT = ITSTOT - 1
      IF ( .NOT. MATZ) LOR1 = LD
C     MAIN LOOP
      DO 540 K = L, NA
         NOTLAS = K .NE. NA .AND. ISH .EQ. 2
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX(K-1,L)
         LL = MIN(EN,K1+ISH)
         IF (NOTLAS) GO TO 360
C        ZERO A(K+1,K-1)
         IF (K.EQ.L) GO TO 320
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
  320    S = ABS(A1) + ABS(A2)
         IF (S.EQ.ZERO) GO TO 80
         U1 = A1/S
         U2 = A2/S
         R = SQRT(U1*U1+U2*U2)
         IF (U1.LT.ZERO) R = -R
         V1 = -(U1+R)/R
         V2 = -U2/R
         U2 = V2/V1
C
         DO 340 J = KM1, ENORN
            T = A(K,J) + U2*A(K1,J)
            A(K,J) = A(K,J) + T*V1
            A(K1,J) = A(K1,J) + T*V2
            T = B(K,J) + U2*B(K1,J)
            B(K,J) = B(K,J) + T*V1
            B(K1,J) = B(K1,J) + T*V2
  340    CONTINUE
C
         IF (K.NE.L) A(K1,KM1) = ZERO
         GO TO 480
C        ZERO A(K+1,K-1) AND A(K+2,K-1)
  360    IF (K.EQ.L) GO TO 380
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
         A3 = A(K2,KM1)
  380    S = ABS(A1) + ABS(A2) + ABS(A3)
         IF (S.EQ.ZERO) GO TO 540
         U1 = A1/S
         U2 = A2/S
         U3 = A3/S
         R = SQRT(U1*U1+U2*U2+U3*U3)
         IF (U1.LT.ZERO) R = -R
         V1 = -(U1+R)/R
         V2 = -U2/R
         V3 = -U3/R
         U2 = V2/V1
         U3 = V3/V1
C
         DO 400 J = KM1, ENORN
            T = A(K,J) + U2*A(K1,J) + U3*A(K2,J)
            A(K,J) = A(K,J) + T*V1
            A(K1,J) = A(K1,J) + T*V2
            A(K2,J) = A(K2,J) + T*V3
            T = B(K,J) + U2*B(K1,J) + U3*B(K2,J)
            B(K,J) = B(K,J) + T*V1
            B(K1,J) = B(K1,J) + T*V2
            B(K2,J) = B(K2,J) + T*V3
  400    CONTINUE
C
         IF (K.EQ.L) GO TO 420
         A(K1,KM1) = ZERO
         A(K2,KM1) = ZERO
C        ZERO B(K+2,K+1) AND B(K+2,K)
  420    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))
         IF (S.EQ.ZERO) GO TO 480
         U1 = B(K2,K2)/S
         U2 = B(K2,K1)/S
         U3 = B(K2,K)/S
         R = SQRT(U1*U1+U2*U2+U3*U3)
         IF (U1.LT.ZERO) R = -R
         V1 = -(U1+R)/R
         V2 = -U2/R
         V3 = -U3/R
         U2 = V2/V1
         U3 = V3/V1
C
         DO 440 I = LOR1, LL
            T = A(I,K2) + U2*A(I,K1) + U3*A(I,K)
            A(I,K2) = A(I,K2) + T*V1
            A(I,K1) = A(I,K1) + T*V2
            A(I,K) = A(I,K) + T*V3
            T = B(I,K2) + U2*B(I,K1) + U3*B(I,K)
            B(I,K2) = B(I,K2) + T*V1
            B(I,K1) = B(I,K1) + T*V2
            B(I,K) = B(I,K) + T*V3
  440    CONTINUE
C
         B(K2,K) = ZERO
         B(K2,K1) = ZERO
         IF ( .NOT. MATZ) GO TO 480
C
         DO 460 I = 1, N
            T = Z(I,K2) + U2*Z(I,K1) + U3*Z(I,K)
            Z(I,K2) = Z(I,K2) + T*V1
            Z(I,K1) = Z(I,K1) + T*V2
            Z(I,K) = Z(I,K) + T*V3
  460    CONTINUE
C        ZERO B(K+1,K)
  480    S = ABS(B(K1,K1)) + ABS(B(K1,K))
         IF (S.EQ.ZERO) GO TO 540
         U1 = B(K1,K1)/S
         U2 = B(K1,K)/S
         R = SQRT(U1*U1+U2*U2)
         IF (U1.LT.ZERO) R = -R
         V1 = -(U1+R)/R
         V2 = -U2/R
         U2 = V2/V1
C
         DO 500 I = LOR1, LL
            T = A(I,K1) + U2*A(I,K)
            A(I,K1) = A(I,K1) + T*V1
            A(I,K) = A(I,K) + T*V2
            T = B(I,K1) + U2*B(I,K)
            B(I,K1) = B(I,K1) + T*V1
            B(I,K) = B(I,K) + T*V2
  500    CONTINUE
C
         B(K1,K) = ZERO
         IF ( .NOT. MATZ) GO TO 540
C
         DO 520 I = 1, N
            T = Z(I,K1) + U2*Z(I,K)
            Z(I,K1) = Z(I,K1) + T*V1
            Z(I,K) = Z(I,K) + T*V2
  520    CONTINUE
C
  540 CONTINUE
C     END QZ STEP
      GO TO 80
C     SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT
C     HAS BECOME NEGLIGIBLE AFTER 30*N ITERATIONS ALTOGETHER
  560 IERR = P01ABF(ISAVE,EN,SRNAME,0,P01REC)
C     SAVE EPSB FOR USE BY F02BJY AND F02BJZ
  580 IF (N.GT.1) B(N,1) = EPSB
      RETURN
      END
