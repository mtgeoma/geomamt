      SUBROUTINE F01BUF(N,M1,K,A,IA,W,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     U*L*D*L**T*U**T DECOMPOSITION OF AN N*N POSITIVE DEFINITE
C     MATRIX A OF BAND WIDTH 2*M+1 WITH M1=M+1. K DETERMINES
C     THE SWITCHOVER POINT FROM U TO L AND MUST BE SUCH THAT
C     K.LE.N, K.GE.M. IF K=N THEN AN L*D*L**T DECOMPOSITION IS
C     OBTAINED. THE UPPER TRIANGULAR PART OF A IS
C     PRESENTED AS AN M1*N ARRAY A WITH THE DIAGONAL
C     ELEMENTS OF A IN THE M1-TH ROW. THE SPARE
C     LOCATIONS IN THE TOP LEFT ARE NOT USED. U, L AND D
C     ARE OVERWRITTEN ON A.
C     D(I) OVERWRITES A(I,I), L(I,J) OVERWRITES A(J,I) AND U(I,J)
C     OVERWRITES A(I,J), THE ELEMENT OF A(J,I) BEING STORED IN
C     POSITION (J+M+1-I,I) OF ARRAY A.
C     THE ALGORITHM WILL FAIL IF ANY COMPUTED ELEMENT OF D IS SUCH
C     THAT D(I).LE.0.0.
C     AN ADDITIONAL WORK SPACE OF M1 ELEMENTS IS USED.
C
C     PDK, DNAC, NPL, TEDDINGTON. DEC 1977.
C     NPL DNAC LIBRARY SUBROUTINE ULBAND.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BUF')
C     .. Scalar Arguments ..
      INTEGER           IA, IFAIL, K, M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), W(M1)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, ONE, U, ZERO
      INTEGER           I, I1, IM, IP, IR, IR1, ISAVE, J, J1, K1, M, MM,
     *                  MR
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      M = M1 - 1
      IF (K.LT.M .OR. K.GT.N) GO TO 200
      IF (K.EQ.N) GO TO 100
      IFAIL = 2
      K1 = K + 1
      DO 80 IR1 = K1, N
         IR = K1 + N - IR1
         D = A(M1,IR)
         IF (D.LE.ZERO) GO TO 200
         D = ONE/D
         IF (M.EQ.0) GO TO 80
         DO 20 J = 1, M
            W(J) = A(J,IR)
   20    CONTINUE
         DO 60 J = 1, M
            U = W(J)*D
            IP = IR
            IM = M1
            J1 = J + 1
            DO 40 I = J1, M1
               IP = IP - 1
               IM = IM - 1
               A(I,IP) = A(I,IP) - U*W(IM)
   40       CONTINUE
            A(J,IR) = U
   60    CONTINUE
   80 CONTINUE
  100 IF (K.EQ.0) GO TO 180
      IFAIL = 3
      MM = M
      DO 160 IR = 1, K
         D = A(M1,IR)
         IF (D.LE.ZERO) GO TO 200
         D = ONE/D
         IF (MM.LT.K) MM = MM + 1
         MR = M + IR - MM
         IR1 = IR + 1
         IF (MM.LT.IR1) GO TO 160
         DO 140 I1 = IR1, MM
            I = MM + IR1 - I1
            MR = MR + 1
            U = A(MR,I)*D
            IP = MR
            IM = M1
            DO 120 J = IR1, I
               IP = IP + 1
               IM = IM - 1
               A(IP,I) = A(IP,I) - U*A(IM,J)
  120       CONTINUE
            A(MR,I) = U
  140    CONTINUE
  160 CONTINUE
  180 IFAIL = 0
      RETURN
  200 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
