      SUBROUTINE E02DDY(A,F,N,M,NA,TOL,C,SQ,RANK,AA,FF,H)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     Subroutine E02DDY finds the minimum norm solution of a least-
C     squares problem in case of rank deficiency.
C
C     Input parameters:
C     A : array, which contains the non-zero elements of the observation
C         matrix after triangularization by Givens transformations.
C     F : array, which contains the transformed right hand side.
C     N : integer,wich contains the dimension of A.
C     M : integer, which denotes the bandwidth of A.
C     TOL : real value, giving a threshold to determine the rank of A.
C
C     Output parameters:
C     C : array, which contains the minimum norm solution.
C     SQ : real value, giving the contribution of reducing the rank
C        to the sum of squared residuals.
C     RANK : integer, which contains the rank of matrix A.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SQ, TOL
      INTEGER           M, N, NA, RANK
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NA,N+NA), AA(M,N), C(N), F(N), FF(N), H(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJ1, COS, FAC, PIV, SIN, STOR1, STOR2, STOR3,
     *                  STORE
      INTEGER           I, I1, I2, II, IJ, J, J1, J2, JJ, K, KK, M1, NL
C     .. Local Arrays ..
      DOUBLE PRECISION  YI(1)
C     .. External Subroutines ..
      EXTERNAL          F06AAF, DROT, DTBSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, SIGN
C     .. Executable Statements ..
      M1 = M - 1
C     The rank deficiency NL is considered to be the number of
C     sufficient small diagonal elements of A.
      NL = 0
      SQ = ZERO
      DO 100 I = 1, N
         IF (A(M,I).LE.TOL) THEN
C           If a sufficiently small diagonal element is found, we put
C           it to zero. The remainder of the row corresponding to that
C           zero diagonal element is then rotated into triangle by
C           Givens rotations. The rank deficiency is increased by one.
            NL = NL + 1
            IF (I.NE.N) THEN
               YI(1) = F(I)
               DO 20 J = 1, M1
                  H(J) = A(M-J,J+I)
   20          CONTINUE
               H(M) = ZERO
               I1 = I + 1
               DO 80 II = I1, N
                  I2 = MIN(N-II,M1)
                  PIV = ABS(H(1))
                  IF (PIV.NE.ZERO) THEN
                     AJ1 = ABS(A(M,II))
                     CALL F06AAF(AJ1,PIV,COS,SIN)
                     COS = SIGN(COS,A(M,II))
                     SIN = SIGN(SIN,H(1))
                     A(M,II) = AJ1
                     CALL DROT(1,YI(1),1,F(II),1,COS,-SIN)
                     CALL DROT(I2,H(2),1,A(M-1,II+1),NA-1,COS,-SIN)
                     DO 40 J = 1, I2
                        H(J) = H(J+1)
   40                CONTINUE
                  ELSE
                     DO 60 J = 1, I2
                        H(J) = H(J+1)
   60                CONTINUE
                  END IF
                  H(I2+1) = ZERO
   80          CONTINUE
C              Add to the sum of squared residuals the contribution
C              of deleting the row with small diagonal element.
               SQ = SQ + YI(1)*YI(1)
            END IF
         END IF
  100 CONTINUE
C     RANK denotes the rank of A.
      RANK = N - NL
C     Let B denote the (RANK*N) upper trapezoidal matrix which can be
C     obtained from the (N*N) upper triangular matrix A by deleting
C     the rows and interchanging the columns corresponding to a zero
C     diagonal element. If this matrix is factorized using Givens
C     transformations as  B = (R) (U)  where
C     R is a (RANK*RANK) upper triangular matrix,
C     U is a (RANK*N) orthonormal matrix
C     then the minimal least-squares solution C is given by C = B' V,
C     where V is the solution of the system  (R) (R)' V = G  and
C     G denotes the vector obtained from the old right hand side F, by
C     removing the elements corresponding to a zero diagonal
C     element of A.
C     Initialization.
      DO 140 I = 1, M
         DO 120 J = M - I + 1, MIN(M-I+RANK,N)
            AA(I,J) = ZERO
  120    CONTINUE
  140 CONTINUE
C     Form in AA the upper triangular matrix obtained from A by
C     removing rows and columns with zero diagonal elements. Form in FF
C     the new right hand side by removing the elements of the old right
C     hand side corresponding to a deleted row.
      II = 0
      DO 180 I = 1, N
         IF (A(M,I).GT.TOL) THEN
            II = II + 1
            FF(II) = F(I)
            AA(M,II) = A(M,I)
            JJ = II
            KK = 1
            J = I
            J1 = MIN(J-1,M1)
            DO 160 K = 1, J1
               J = J - 1
               IF (A(M,J).GT.TOL) THEN
                  KK = KK + 1
                  JJ = JJ - 1
                  IF (K+J.LE.N) THEN
                     AA(M+1-KK,KK-1+JJ) = A(M-K,K+J)
                  ELSE
                     AA(M+1-KK,KK-1+JJ) = ZERO
                  END IF
               END IF
  160       CONTINUE
         END IF
  180 CONTINUE
C     Form successively in H the columns of A with a zero diagonal
C     element.
      II = 0
      DO 300 I = 1, N
         II = II + 1
         IF (A(M,I).LE.TOL) THEN
            II = II - 1
            IF (II.NE.0) THEN
               JJ = 1
               J = I
               DO 200 K = 1, MIN(J-1,M1)
                  J = J - 1
                  IF (A(M,J).GT.TOL) THEN
                     IF (K+J.LE.N) THEN
                        H(JJ) = A(M-K,K+J)
                     ELSE
                        H(JJ) = ZERO
                     END IF
                     JJ = JJ + 1
                  END IF
  200          CONTINUE
               DO 220 KK = JJ, M
                  H(KK) = ZERO
  220          CONTINUE
C              Rotate this column into AA by Givens transformations.
               JJ = II
               DO 280 I1 = 1, II
                  J1 = MIN(JJ-1,M1)
                  PIV = ABS(H(1))
                  IF (PIV.NE.ZERO) THEN
                     AJ1 = ABS(AA(M,JJ))
                     CALL F06AAF(AJ1,PIV,COS,SIN)
                     COS = SIGN(COS,AA(M,JJ))
                     SIN = SIGN(SIN,H(1))
                     AA(M,JJ) = AJ1
                     CALL DROT(J1,H(2),1,AA(M-J1,JJ),-1,COS,-SIN)
                     DO 240 J2 = 1, J1
                        H(J2) = H(J2+1)
  240                CONTINUE
                  ELSE
                     DO 260 J2 = 1, J1
                        H(J2) = H(J2+1)
  260                CONTINUE
                  END IF
                  JJ = JJ - 1
                  H(J1+1) = ZERO
  280          CONTINUE
            END IF
         END IF
  300 CONTINUE
C     Solve the system (AA) (F1) = FF
      CALL DTBSV('U','N','N',RANK,M-1,AA,M,FF,1)
C     Solve the system  (AA)' (F2) = F1
      CALL DTBSV('U','T','N',RANK,M-1,AA,M,FF,1)
C     Premultiply F2 by the transpose of A.
      K = 0
      DO 340 I = 1, N
         STORE = ZERO
         IF (A(M,I).GT.TOL) K = K + 1
         J1 = MIN(I,M)
         KK = K
         IJ = I + 1
         DO 320 J = 1, J1
            IJ = IJ - 1
            IF (A(M,IJ).GT.TOL) THEN
               IF (J-1+IJ.LE.N) THEN
                  STOR1 = A(M+1-J,J-1+IJ)
               ELSE
                  STOR1 = ZERO
               END IF
               STOR2 = FF(KK)
               STORE = STORE + STOR1*STOR2
               KK = KK - 1
            END IF
  320    CONTINUE
         C(I) = STORE
  340 CONTINUE
C     Add to the sum of squared residuals the contribution of putting
C     to zero the small diagonal elements of matrix (A).
      STOR3 = ZERO
      DO 380 I = 1, N
         IF (A(M,I).LE.TOL) THEN
            STORE = F(I)
            DO 360 J = 1, MIN(N-I,M1)
               IJ = I + J
               STOR1 = C(IJ)
               STOR2 = A(M-J,IJ)
               STORE = STORE - STOR1*STOR2
  360       CONTINUE
            STOR1 = A(M,I)
            STOR2 = C(I)
            STOR1 = STOR1*STOR2
            STOR3 = STOR3 + STOR1*(STOR1-STORE-STORE)
         END IF
  380 CONTINUE
      FAC = STOR3
      SQ = SQ + FAC
      RETURN
      END
