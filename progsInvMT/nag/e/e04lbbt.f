      SUBROUTINE E04LBB(N,LEL,IS,EL,P)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     SOLVES LT*P = ES FOR P, WHERE LT IS THE TRANSPOSE OF L, AN N*N
C     UNIT LOWER-TRIANGULAR MATRIX STORED BY ROWS OMITTING THE UNIT
C     DIAGONAL IN THE ARRAY EL(I), I = 1(1)N*(N-1)/2. ES IS AN
C     N-VECTOR WITH ITS ISTH ELEMENT UNITY AND ALL OTHERS ZERO.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, SHEDMOND R.
C     GRAHAM AND MARGARET H. WRIGHT, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           IS, LEL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  EL(LEL), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, IP1, IQ, IS1, ISP1, IT, J, J1, K, L, NP1
C     .. Executable Statements ..
      IS1 = IS
      IT = IS1*(IS1-1)/2 + 1
      NP1 = N + 1
      ISP1 = IS1 + 1
      DO 80 J = 1, N
         I = NP1 - J
         IF (I.GT.IS1) GO TO 60
         IQ = IT
         IT = IT - 1
         SUM = 0.0D+0
         IF (I.EQ.IS1) SUM = 1.0D+0
         IP1 = I + 1
         IF (IP1.GT.IS1) GO TO 40
         L = ISP1 - IP1
         DO 20 J1 = 1, L
            K = ISP1 - J1
            SUM = SUM - P(K)*EL(IQ)
            IQ = IQ + 2 - K
   20    CONTINUE
   40    P(I) = SUM
         GO TO 80
   60    P(I) = 0.0D+0
   80 CONTINUE
      RETURN
C
C     END OF E04LBB (ELTSOL)
C
      END
