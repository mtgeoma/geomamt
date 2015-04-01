      SUBROUTINE F02BJZ(N,A,IA,B,IB,ALFR,ALFI,BETA,Z,IZ)
C     MARK 6 RELEASE  NAG COPYRIGHT 1977
C     MARK 9A REVISED. IER-359 (NOV 1981)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS SUBROUTINE IS THE OPTIONAL FOURTH STEP OF THE QZ
C     ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL MATRICES, ONE OF THEM
C     IN
C     QUASI-TRIANGULAR FORM (IN WHICH EACH 2-BY-2 BLOCK CORRESPONDS
C     TO
C     A PAIR OF COMPLEX EIGENVALUES) AND THE OTHER IN UPPER
C     TRIANGULAR
C     FORM.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND
C     TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE
C     SYSTEM.
C     IT IS USUALLY PRECEDED BY F02BJW, F02BJX AND F02BJY.
C
C     ON INPUT
C
C     N IS THE ORDER OF THE MATRICES.
C
C     A CONTAINS A REAL UPPER QUASI-TRIANGULAR MATRIX.
C
C     IA MUST SPECIFY THE FIRST DIMENSION OF THE ARRAY A AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IA.GE.N
C
C     B CONTAINS A REAL UPPER TRIANGULAR MATRIX.  IN ADDITION,
C     LOCATION B(N,1) CONTAINS THE TOLERANCE QUANTITY (EPSB)
C     COMPUTED AND SAVED IN F02BJX.
C
C     IB MUST SPECIFY THE FIRST DIMENSION OF THE ARRAY B AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IB.GE.N
C
C     ALFR, ALFI, AND BETA  ARE VECTORS WITH COMPONENTS WHOSE
C     RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C     EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM F02BJY.
C
C     Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C     REDUCTIONS BY F02BJW, F02BJX AND F02BJY, IF PERFORMED.
C     IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C     DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     IZ MUST SPECIFY THE FIRST DIMENSION OF THE ARRAY Z AS
C     DECLARED IN THE CALLING (SUB)PROGRAM.  IZ.GE.N
C
C     ON OUTPUT
C
C     A IS UNALTERED.  ITS SUBDIAGONAL ELEMENTS PROVIDE INFORMATION
C     ABOUT THE STORAGE OF THE COMPLEX EIGENVECTORS.
C
C     B HAS BEEN DESTROYED.
C
C     ALFR, ALFI, AND BETA ARE UNALTERED.
C
C     Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C     IF ALFI(I) .EQ. 0.0, THE I-TH EIGENVALUE IS REAL AND
C     THE I-TH COLUMN OF Z CONTAINS ITS EIGENVECTOR.
C     IF ALFI(I) .NE. 0.0, THE I-TH EIGENVALUE IS COMPLEX.
C     IF ALFI(I) .GT. 0.0, THE EIGENVALUE IS THE FIRST OF
C     A COMPLEX PAIR AND THE I-TH AND (I+1)-TH COLUMNS
C     OF Z CONTAIN ITS EIGENVECTOR.
C     IF ALFI(I) .LT. 0.0, THE EIGENVALUE IS THE SECOND OF
C     A COMPLEX PAIR AND THE (I-1)-TH AND I-TH COLUMNS
C     OF Z CONTAIN THE CONJUGATE OF ITS EIGENVECTOR.
C     EACH EIGENVECTOR IS NORMALIZED SO THAT THE COMPONENT WITH
C     LARGEST MODULUS IS REAL AND THE SUM OF SQUARES OF THE
C     MODULI OF THE COMPONENTS IS EQUAL TO 1.0
C
C     IMPLEMENTED FROM EISPACK BY
C     W. PHILLIPS OXFORD UNIVERSITY COMPUTING SERVICE.
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), ALFI(N), ALFR(N), B(IB,N), BETA(N),
     *                  Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALFM, ALMI, ALMR, BETM, D, DI, DR, EPSB, ONE, Q,
     *                  R, RA, RR, S, SA, T, T1, T2, TI, TR, W, W1, X,
     *                  X1, Y, Z1, ZERO, ZZ
      INTEGER           EN, ENM2, I, II, ISW, J, JJ, K, M, NA, NN
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, F04JGW
      EXTERNAL          A02ABF, F04JGW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      EPSB = B(N,1)
      ISW = 1
C     FOR EN=N STEP -1 UNTIL 1 DO --
      DO 380 NN = 1, N
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW.EQ.2) GO TO 360
         IF (ALFI(EN).NE.ZERO) GO TO 140
C        REAL VECTOR
         M = EN
         B(EN,EN) = ONE
         IF (NA.EQ.0) GO TO 380
         ALFM = ALFR(M)
         BETM = BETA(M)
C        FOR I=EN-1 STEP -1 UNTIL 1 DO --
         DO 120 II = 1, NA
            I = EN - II
            W = BETM*A(I,I) - ALFM*B(I,I)
            R = ZERO
C
            DO 20 J = M, EN
               R = R + (BETM*A(I,J)-ALFM*B(I,J))*B(J,EN)
   20       CONTINUE
C
            IF (I.EQ.1 .OR. ISW.EQ.2) GO TO 40
            IF (BETM*A(I,I-1).EQ.ZERO) GO TO 40
            ZZ = W
            S = R
            GO TO 100
   40       M = I
            IF (ISW.EQ.2) GO TO 60
C           REAL 1-BY-1 BLOCK
            T = W
            IF (W.EQ.ZERO) T = EPSB
            B(I,EN) = -R/T
            GO TO 120
C           REAL 2-BY-2 BLOCK
   60       X = BETM*A(I,I+1) - ALFM*B(I,I+1)
            Y = BETM*A(I+1,I)
            Q = W*ZZ - X*Y
            T = (X*S-ZZ*R)/Q
            B(I,EN) = T
            IF (ABS(X).LE.ABS(ZZ)) GO TO 80
            B(I+1,EN) = (-R-W*T)/X
            GO TO 100
   80       B(I+1,EN) = (-S-Y*T)/ZZ
  100       ISW = 3 - ISW
  120    CONTINUE
C        END REAL VECTOR
         GO TO 380
C        COMPLEX VECTOR
  140    M = NA
         ALMR = ALFR(M)
         ALMI = ALFI(M)
         BETM = BETA(M)
C        LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C        EIGENVECTOR MATRIX IS TRIANGULAR
         Y = BETM*A(EN,NA)
         B(NA,NA) = -ALMI*B(EN,EN)/Y
         B(NA,EN) = (ALMR*B(EN,EN)-BETM*A(EN,EN))/Y
         B(EN,NA) = ZERO
         B(EN,EN) = ONE
         ENM2 = NA - 1
         IF (ENM2.EQ.0) GO TO 360
C        FOR I=EN-2 STEP -1 UNTIL 1 DO --
         DO 340 II = 1, ENM2
            I = NA - II
            W = BETM*A(I,I) - ALMR*B(I,I)
            W1 = -ALMI*B(I,I)
            RA = ZERO
            SA = ZERO
C
            DO 160 J = M, EN
               X = BETM*A(I,J) - ALMR*B(I,J)
               X1 = -ALMI*B(I,J)
               RA = RA + X*B(J,NA) - X1*B(J,EN)
               SA = SA + X*B(J,EN) + X1*B(J,NA)
  160       CONTINUE
C
            IF (I.EQ.1 .OR. ISW.EQ.2) GO TO 180
            IF (BETM*A(I,I-1).EQ.ZERO) GO TO 180
            ZZ = W
            Z1 = W1
            R = RA
            S = SA
            ISW = 2
            GO TO 340
  180       M = I
            IF (ISW.EQ.2) GO TO 260
C           COMPLEX 1-BY-1 BLOCK
            TR = -RA
            TI = -SA
  200       DR = W
            DI = W1
C           COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI)
  220       IF (ABS(DI).GT.ABS(DR)) GO TO 240
            RR = DI/DR
            D = DR + DI*RR
            T1 = (TR+TI*RR)/D
            T2 = (TI-TR*RR)/D
            GO TO (320,280) ISW
  240       RR = DR/DI
            D = DR*RR + DI
            T1 = (TR*RR+TI)/D
            T2 = (TI*RR-TR)/D
            GO TO (320,280) ISW
C           COMPLEX 2-BY-2 BLOCK
  260       X = BETM*A(I,I+1) - ALMR*B(I,I+1)
            X1 = -ALMI*B(I,I+1)
            Y = BETM*A(I+1,I)
            TR = Y*RA - W*R + W1*S
            TI = Y*SA - W*S - W1*R
            DR = W*ZZ - W1*Z1 - X*Y
            DI = W*Z1 + W1*ZZ - X1*Y
            IF (DR.EQ.ZERO .AND. DI.EQ.ZERO) DR = EPSB
            GO TO 220
  280       B(I+1,NA) = T1
            B(I+1,EN) = T2
            ISW = 1
            IF (ABS(Y).GT.ABS(W)+ABS(W1)) GO TO 300
            TR = -RA - X*B(I+1,NA) + X1*B(I+1,EN)
            TI = -SA - X*B(I+1,EN) - X1*B(I+1,NA)
            GO TO 200
  300       T1 = (-R-ZZ*B(I+1,NA)+Z1*B(I+1,EN))/Y
            T2 = (-S-ZZ*B(I+1,EN)-Z1*B(I+1,NA))/Y
  320       B(I,NA) = T1
            B(I,EN) = T2
  340    CONTINUE
C        END COMPLEX VECTOR
  360    ISW = 3 - ISW
  380 CONTINUE
C     END BACK SUBSTITUTION.
C     TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C     FOR J=N STEP -1 UNTIL 1 DO --
      DO 440 JJ = 1, N
         J = N + 1 - JJ
C
         DO 420 I = 1, N
            ZZ = ZERO
C
            DO 400 K = 1, J
               ZZ = ZZ + Z(I,K)*B(K,J)
  400       CONTINUE
C
            Z(I,J) = ZZ
  420    CONTINUE
  440 CONTINUE
C     NORMALIZE SO THAT, FOR EACH VECTOR, THE COMPONENT
C     WITH LARGEST MODULUS IS REAL AND THE SUM OF SQUARES OF
C     THE MODULI OF THE COMPONENTS IS 1.0
      J = 0
  460 J = J + 1
      IF (ALFI(J).NE.ZERO) GO TO 500
      D = F04JGW(N,Z(1,J))
      DO 480 I = 1, N
         Z(I,J) = Z(I,J)/D
  480 CONTINUE
      GO TO 560
  500 R = ZERO
      DO 520 I = 1, N
         W = A02ABF(Z(I,J),Z(I,J+1))
         IF (W.LE.R) GO TO 520
         R = W
         X = Z(I,J)
         Y = -Z(I,J+1)
  520 CONTINUE
      D = A02ABF(F04JGW(N,Z(1,J)),F04JGW(N,Z(1,J+1)))*A02ABF(X,Y)
      DO 540 I = 1, N
         W = Z(I,J)
         Z(I,J) = (Z(I,J)*X-Z(I,J+1)*Y)/D
         Z(I,J+1) = (Y*W+X*Z(I,J+1))/D
  540 CONTINUE
      J = J + 1
  560 IF (J.LT.N) GO TO 460
C
      RETURN
      END
