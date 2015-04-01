      SUBROUTINE F02GJX(N,AR,IAR,AI,IAI,BR,IBR,BI,IBI,MATZ,ZR,IZR,ZI,
     *                  IZI)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11A REVISED. IER-450 (JUN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-516 (AUG 1986).
C
C     THIS ROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE
C     PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS ROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND
C     NON-NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER
C     TRIANGULAR FORM USING UNITARY TRANSFORMATIONS. IT IS
C     USUALLY FOLLOWED BY F02GJY AND, POSSIBLY, F02GJZ.
C
C     ON INPUT
C
C     N IS THE ORDER OF THE MATRICES.
C
C     A = (AR,AI) CONTAINS A COMPLEX GENERAL MATRIX.
C
C     B = (BR,BI) CONTAINS A COMPLEX GENERAL MATRIX.
C
C     MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C     ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING EIGENVECTORS,
C     AND TO .FALSE. OTHERWISE.
C
C     IAR,IAI,IBR,IBI,IZR AND IZI MUST SPECIFY THE FIRST DIMENSIONS
C     OF THE ARRAYS AR,AI,BR,BI,ZR AND ZI AS DECLARED IN THE CALLING
C     (SUB)PROGRAM.  IAR,IAI,IBR,IBI,IZR,IZI.GE.N  .
C
C     ON OUTPUT
C
C     A HAS BEEN REDUCED TO UPPER HESSENBERG FORM. THE ELEMENTS
C     BELOW THE FIRST SUBDIAGONAL HAVE BEEN REDUCED TO ZERO,
C     AND THE SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL
C     (AND NON-NEGATIVE).
C
C     B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM. THE ELEMENTS
C     BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO.
C
C     Z = (ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C     TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C     OTHERWISE, Z IS NOT REFERENCED.
C
C     IMPLEMENTED FROM EISPACK BY
C     W. PHILLIPS OXFORD UNIVERSITY COMPUTING SERVICE.
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IBI, IBR, IZI, IZR, N
      LOGICAL           MATZ
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(IBI,N), BR(IBR,N),
     *                  ZI(IZI,N), ZR(IZR,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, R, RHO, S, T, TI, U1, U1I, U2, XI, XR, YI,
     *                  YR, ZERO
      INTEGER           I, J, K, K1, L, L1, LB, NK1, NM1
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF
      EXTERNAL          A02ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      IF ( .NOT. MATZ) GO TO 60
C
C     ---------- INITIALIZE Z ----------
      DO 40 I = 1, N
C
         DO 20 J = 1, N
            ZR(I,J) = ZERO
            ZI(I,J) = ZERO
   20    CONTINUE
C
         ZR(I,I) = ONE
   40 CONTINUE
C     ---------- REDUCE B TO UPPER TRIANGULAR FORM WITH
C     TEMPORARILY REAL DIAGONAL ELEMENTS ----------
   60 IF (N.LE.1) GO TO 500
      NM1 = N - 1
C
      DO 300 L = 1, NM1
         L1 = L + 1
         S = ZERO
C
         DO 80 I = L, N
            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   80    CONTINUE
C
         IF (S.EQ.ZERO) GO TO 300
         RHO = ZERO
C
         DO 100 I = L, N
            BR(I,L) = BR(I,L)/S
            BI(I,L) = BI(I,L)/S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
  100    CONTINUE
C
         R = SQRT(RHO)
         XR = A02ABF(BR(L,L),BI(L,L))
         IF (XR.EQ.ZERO) GO TO 120
         RHO = RHO + XR*R
         U1 = -BR(L,L)/XR
         U1I = -BI(L,L)/XR
         YR = R/XR + ONE
         BR(L,L) = YR*BR(L,L)
         BI(L,L) = YR*BI(L,L)
         GO TO 140
C
  120    BR(L,L) = R
         U1 = -ONE
         U1I = ZERO
C
  140    DO 200 J = L1, N
            T = ZERO
            TI = ZERO
C
            DO 160 I = L, N
               T = T + BR(I,L)*BR(I,J) + BI(I,L)*BI(I,J)
               TI = TI + BR(I,L)*BI(I,J) - BI(I,L)*BR(I,J)
  160       CONTINUE
C
            T = T/RHO
            TI = TI/RHO
C
            DO 180 I = L, N
               BR(I,J) = BR(I,J) - T*BR(I,L) + TI*BI(I,L)
               BI(I,J) = BI(I,J) - T*BI(I,L) - TI*BR(I,L)
  180       CONTINUE
C
            XI = U1*BI(L,J) - U1I*BR(L,J)
            BR(L,J) = U1*BR(L,J) + U1I*BI(L,J)
            BI(L,J) = XI
  200    CONTINUE
C
         DO 260 J = 1, N
            T = ZERO
            TI = ZERO
C
            DO 220 I = L, N
               T = T + BR(I,L)*AR(I,J) + BI(I,L)*AI(I,J)
               TI = TI + BR(I,L)*AI(I,J) - BI(I,L)*AR(I,J)
  220       CONTINUE
C
            T = T/RHO
            TI = TI/RHO
C
            DO 240 I = L, N
               AR(I,J) = AR(I,J) - T*BR(I,L) + TI*BI(I,L)
               AI(I,J) = AI(I,J) - T*BI(I,L) - TI*BR(I,L)
  240       CONTINUE
C
            XI = U1*AI(L,J) - U1I*AR(L,J)
            AR(L,J) = U1*AR(L,J) + U1I*AI(L,J)
            AI(L,J) = XI
  260    CONTINUE
C
         BR(L,L) = R*S
         BI(L,L) = ZERO
C
         DO 280 I = L1, N
            BR(I,L) = ZERO
            BI(I,L) = ZERO
  280    CONTINUE
C
  300 CONTINUE
C     REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C     ELEMENTS, WHILE KEEPING B TRIANGULAR
      DO 480 K = 1, NM1
         K1 = K + 1
C        SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL
         R = A02ABF(AR(N,K),AI(N,K))
         IF (R.EQ.ZERO) GO TO 340
         U1 = AR(N,K)/R
         U1I = AI(N,K)/R
         AR(N,K) = R
         AI(N,K) = ZERO
C
         DO 320 J = K1, N
            XI = U1*AI(N,J) - U1I*AR(N,J)
            AR(N,J) = U1*AR(N,J) + U1I*AI(N,J)
            AI(N,J) = XI
  320    CONTINUE
C
         XI = U1*BI(N,N) - U1I*BR(N,N)
         BR(N,N) = U1*BR(N,N) + U1I*BI(N,N)
         BI(N,N) = XI
  340    IF (K.EQ.NM1) GO TO 500
         NK1 = NM1 - K
C        ---------- FOR L=N-1 STEP -1 UNTIL K+1 DO -- ----------
         DO 460 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C           ---------- ZERO A(L+1,K) ----------
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
            IF (S.EQ.ZERO) GO TO 460
            U1 = AR(L,K)/S
            U1I = AI(L,K)/S
            U2 = AR(L1,K)/S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1/R
            U1I = U1I/R
            U2 = U2/R
            AR(L,K) = R*S
            AI(L,K) = ZERO
            AR(L1,K) = ZERO
C
            DO 360 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1*XR + U1I*XI + U2*YR
               AI(L,J) = U1*XI - U1I*XR + U2*YI
               AR(L1,J) = U1*YR - U1I*YI - U2*XR
               AI(L1,J) = U1*YI + U1I*YR - U2*XI
  360       CONTINUE
C
            XR = BR(L,L)
            BR(L,L) = U1*XR
            BI(L,L) = -U1I*XR
            BR(L1,L) = -U2*XR
C
            DO 380 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1*XR + U1I*XI + U2*YR
               BI(L,J) = U1*XI - U1I*XR + U2*YI
               BR(L1,J) = U1*YR - U1I*YI - U2*XR
               BI(L1,J) = U1*YI + U1I*YR - U2*XI
  380       CONTINUE
C           ---------- ZERO B(L+1,L) ----------
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
            IF (S.EQ.ZERO) GO TO 460
            U1 = BR(L1,L1)/S
            U1I = BI(L1,L1)/S
            U2 = BR(L1,L)/S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1/R
            U1I = U1I/R
            U2 = U2/R
            BR(L1,L1) = R*S
            BI(L1,L1) = ZERO
            BR(L1,L) = ZERO
C
            DO 400 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1*XR + U1I*XI + U2*YR
               BI(I,L1) = U1*XI - U1I*XR + U2*YI
               BR(I,L) = U1*YR - U1I*YI - U2*XR
               BI(I,L) = U1*YI + U1I*YR - U2*XI
  400       CONTINUE
C
            DO 420 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1*XR + U1I*XI + U2*YR
               AI(I,L1) = U1*XI - U1I*XR + U2*YI
               AR(I,L) = U1*YR - U1I*YI - U2*XR
               AI(I,L) = U1*YI + U1I*YR - U2*XI
  420       CONTINUE
C
            IF ( .NOT. MATZ) GO TO 460
C
            DO 440 I = L - K + 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1*XR + U1I*XI + U2*YR
               ZI(I,L1) = U1*XI - U1I*XR + U2*YI
               ZR(I,L) = U1*YR - U1I*YI - U2*XR
               ZI(I,L) = U1*YI + U1I*YR - U2*XI
  440       CONTINUE
C
  460    CONTINUE
C
  480 CONTINUE
C
  500 RETURN
      END
