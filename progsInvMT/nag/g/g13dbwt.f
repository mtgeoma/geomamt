      SUBROUTINE G13DBW(U,V,W,X,Y,Z,NZ,G,NSM,NS,WA,NWA)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        G13DBW CALCULATES U,V,Y, AND Z FROM G,W, AND X.
C
C             T          T  T           T T            T
C        1. WV =G   2. XU =G    3. Z=X-G V    4. Y=W-GU
C
C                   T
C        SOLVE FOR V
C     .. Scalar Arguments ..
      INTEGER           NS, NSM, NWA, NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  G(NS,NS), U(NSM,NSM), V(NSM,NSM), W(NSM,NS),
     *                  WA(NWA), X(NSM,NS), Y(NSM,NSM), Z(NZ,NZ)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, TEMP
      INTEGER           I, IERR, IFAIL1, ILIM, ITR, ITR1, J, J1, LTR,
     *                  NS1, NS2, NSQ, NSQ1
C     .. External Subroutines ..
      EXTERNAL          F04ABF, G13DBT, G13DBV
C     .. Executable Statements ..
      IFAIL1 = 1
      NS1 = NS + 1
      CALL F04ABF(W,NSM,G,NS,NS,NS,V,NSM,WA,WA(NS1),NS,IFAIL1)
C        RESTORE W TO SYMMETRY
      CALL G13DBV(W,NSM,W,NSM,NS)
C        TRANSPOSE G
      NSQ = NS*NS
      ILIM = NS - 1
      DO 40 ITR = 1, ILIM
         LTR = ITR + 1
         DO 20 ITR1 = LTR, NS
            TEMP = G(ITR1,ITR)
            G(ITR1,ITR) = G(ITR,ITR1)
            G(ITR,ITR1) = TEMP
   20    CONTINUE
   40 CONTINUE
C                   T
C        SOLVE FOR U
      IFAIL1 = 1
      CALL F04ABF(X,NSM,G,NS,NS,NS,U,NSM,WA,WA(NS1),NS,IFAIL1)
C        RESTORE X TO SYMMETRY
      CALL G13DBV(X,NSM,X,NSM,NS)
C                   T T
C        CALCULATE G V  AND STORE IN Z
      CALL G13DBT(Z,NZ,G,NS,V,NSM,NS,IERR)
C                     T T
C        CALCULATE X-G V  AND STORE IN Z
      DO 80 J = 1, NS
         DO 60 I = 1, NS
            Z(I,J) = X(I,J) - Z(I,J)
   60    CONTINUE
   80 CONTINUE
C        TRANSPOSE G BACK
      NSQ1 = NSQ + 1
      ILIM = NS - 1
      DO 120 ITR = 1, ILIM
         LTR = ITR + 1
         DO 100 ITR1 = LTR, NS
            TEMP = G(ITR1,ITR)
            G(ITR1,ITR) = G(ITR,ITR1)
            G(ITR,ITR1) = TEMP
  100    CONTINUE
  120 CONTINUE
C                    T
C        CALCULATE GU  AND STORE IN Y
      CALL G13DBT(Y,NSM,G,NS,U,NSM,NS,IERR)
C                     T
C        CACULATE W-GU  AND STORE IN Y
      DO 160 J = 1, NS
         DO 140 I = 1, NS
            Y(I,J) = W(I,J) - Y(I,J)
  140    CONTINUE
  160 CONTINUE
C        TRANSPOSE U AND V
      NS2 = NS - 1
      DO 200 J = 1, NS2
         J1 = J + 1
         DO 180 I = J1, NS
            A1 = U(I,J)
            A2 = U(J,I)
            U(I,J) = A2
            U(J,I) = A1
            A1 = V(I,J)
            A2 = V(J,I)
            V(I,J) = A2
            V(J,I) = A1
  180    CONTINUE
  200 CONTINUE
      RETURN
      END
