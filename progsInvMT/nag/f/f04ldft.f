      SUBROUTINE F04LDF(N,M1,M2,IR,A,IA,AL,IL,IN,B,IB,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-742 (DEC 1989).
C
C     PDK, DNAC, NPL, TEDDINGTON. JUNE 1976.
C     NPL DNAC LIBRARY SUBROUTINE BANSOL.
C     REVISED BY N.A.G. CENTRAL OFFICE, 1979.
C
C     SOLUTION OF THE BAND EQUATIONS AX=B WHERE A,AL AND IN
C     ARE AS SUPPLIED BY SUBROUTINE F01LBF AND B IS AN N*IR MATRIX.
C     ALTHOUGH B MUST BE A TWO-DIMENSIONAL MATRIX IT IS PERMISSIBLE
C     TO PUT IR=1. IA,IL AND IB ARE THE ROW DIMENSIONS OF A, AL AND
C     B IN THE CALLING PROGRAM. X IS OVERWRITTEN ON B.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04LDF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, IL, IR, M1, M2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), AL(IL,*), B(IB,IR)
      INTEGER           IN(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, II, IK, IW, J, JJ, K, KK, M
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IW = MIN(M1+M2+1,N)
      IF (N.LT.1 .OR. M1.LT.0 .OR. M1.GT.N-1 .OR. M2.LT.0 .OR. M2.GT.N-
     *    1 .OR. IA.LT.IW .OR. IL.LT.1 .OR. IL.LT.M1 .OR. IB.LT.N)
     *    GO TO 160
      M = M1
      DO 60 K = 1, N
         IK = K + 1
         M = MIN(N,M+1)
         J = IN(K)
         DO 40 JJ = 1, IR
            X = B(J,JJ)
            B(J,JJ) = B(K,JJ)
            B(K,JJ) = X
            IF (IK.GT.M) GO TO 40
            DO 20 I = IK, M
               II = I - K
               B(I,JJ) = B(I,JJ) - X*AL(II,K)
   20       CONTINUE
   40    CONTINUE
   60 CONTINUE
C
C     FORWARD SUBSTITUTION COMPLETED.
C
      DO 140 K = 1, N
         M = MIN(IW,K)
         I = N + 1 - K
         II = I - 1
         Y = A(1,I)
         DO 120 JJ = 1, IR
            X = B(I,JJ)
            IF (M.EQ.1) GO TO 100
            DO 80 J = 2, M
               KK = J + II
               X = X - A(J,I)*B(KK,JJ)
   80       CONTINUE
  100       B(I,JJ) = X*Y
  120    CONTINUE
  140 CONTINUE
C
C     BACKWARD SUBSTITUTION COMPLETED.
C
      IFAIL = 0
      RETURN
  160 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
