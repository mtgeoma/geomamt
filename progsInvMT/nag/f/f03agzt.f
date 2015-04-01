      SUBROUTINE F03AGZ(N,M,A,IA,RL,IL,MM,D1,ID,LFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C     Originally call F03AGF
C     CHOBANDDET
C     THE LOWER HALF OF A POSITIVE DEFINITE SYMMETRIC BAND MATRIX,
C     A, WITH M LINES ON EITHER SIDE OF THE DIAGONAL IS STORED AS
C     AN N*(M+1) ARRAY, A(I,K), I=1,N, K=1,M+1, A(I,M+1) BEING
C     THE DIAGONAL ELEMENTS. THE CHOLESKY DECOMPOSITION A=LU,
C     WHERE U IS THE TRANSPOSE OF L, IS PERFORMED AND L IS STORED
C     IN RL(I,K) IN THE SAME FORM AS A. THE RECIPROCALS OF THE
C     DIAGONAL ELEMENTS ARE STORED INSTEAD OF THE ELEMENTS
C     THEMSELVES. A IS RETAINED SO THAT THE SOLUTION OBTAINED CAN
C     SUBSEQUENTLY BE IMPROVED. HOWEVER, RL AND A CAN BE IDENTICAL
C     IN THE CALL OF THE SUBROUTINE. THE DETERMINANT, D1* 2.0**ID,
C     OF
C     A IS ALSO COMPUTED. THE SUBROUTINE WILL FAIL IF A, MODIFIED
C     BY ROUNDING ERRORS, IS NOT POSITIVE DEFINITE.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03AGZ')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D1
      INTEGER           IA, ID, IL, LFAIL, M, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,MM), RL(IL,MM)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y
      INTEGER           I, IP, IQ, IR, IS, ISAVE, J, K, M1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 0
      D1 = 1.0D0
      ID = 0
      M1 = M + 1
      DO 220 I = 1, N
         IP = MAX(1,M-I+2)
         IR = I - M + IP - 1
         DO 200 J = IP, M1
            IS = J - 1
            IQ = M + 1 - J + IP
            Y = A(I,J)
            IF (IS.LT.IP) GO TO 40
            DO 20 K = IP, IS
               Y = Y - RL(I,K)*RL(IR,IQ)
               IQ = IQ + 1
   20       CONTINUE
   40       IF (J.EQ.M1) GO TO 60
            RL(I,J) = Y*RL(IR,M1)
            GO TO 180
   60       D1 = D1*Y
            IF (Y.NE.0.0D0) GO TO 80
            ID = 0
            GO TO 140
   80       IF (ABS(D1).LT.1.0D0) GO TO 100
            D1 = D1*0.0625D0
            ID = ID + 4
            GO TO 80
  100       IF (ABS(D1).GE.0.0625D0) GO TO 120
            D1 = D1*16.0D0
            ID = ID - 4
            GO TO 100
  120       IF (Y.GE.0.0D0) GO TO 160
  140       LFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
            RETURN
  160       RL(I,J) = 1.0D0/SQRT(Y)
  180       IR = IR + 1
  200    CONTINUE
  220 CONTINUE
      RETURN
      END
