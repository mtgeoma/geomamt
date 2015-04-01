      SUBROUTINE G13DSW(K,P,Q,M,QQ,IK,COV,IM,DEL,X,D,TEMP,PP,N,PAR,NPAR,
     *                  PARHLD,WORK,INFO,GITAL,G,NOPARS,IERR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IERR, IK, IM, K, M, N, NPAR, P, Q
      LOGICAL           NOPARS
C     .. Array Arguments ..
      DOUBLE PRECISION  COV(IM,M*K*K), D(K), DEL(K,K), G(M*K*K,M*K*K),
     *                  GITAL(K,K), INFO(NPAR+1,NPAR), PAR(NPAR),
     *                  PP(K,K), QQ(IK,K), TEMP(K,M*K), WORK(K),
     *                  X(NPAR,M*K*K)
      LOGICAL           PARHLD(NPAR)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, SUM2, SUM3, TT
      INTEGER           A, B, I, IFAIL, J, K3, L, L2, MK2, NPARN
C     .. External Subroutines ..
      EXTERNAL          F01ADF, F02ABF, G13DSX, G13DSY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
C     This subroutine calculates the covariance matrix of the
C     residual cross-correlation matrices
C
C     IERR is returned as 1 if there is a breakdown in this subroutine
C     and 0 otherwise
C
C     reduce QQ to correlation form and store as DEL
C     (QQ has already been checked for positive-definiteness)
C
      IERR = 0
      DO 40 I = 1, K
         DO 20 J = 1, I - 1
            DEL(I,J) = QQ(I,J)/(SQRT(QQ(I,I))*SQRT(QQ(J,J)))
            DEL(J,I) = DEL(I,J)
   20    CONTINUE
         DEL(I,I) = 1.0D0
   40 CONTINUE
C
      K3 = K*K
      MK2 = M*K3
      TT = 1.0D0/DBLE(N)
      IF (NOPARS) THEN
         CALL G13DSY(K,M,DEL,COV,IM)
         DO 80 J = 1, MK2
            DO 60 I = J, MK2
               COV(I,J) = COV(I,J)*TT
               COV(J,I) = COV(I,J)
   60       CONTINUE
   80    CONTINUE
         RETURN
      END IF
C
C     construct X' matrix and store as X
C
      CALL G13DSX(K,M,P,Q,PAR,NPAR,DEL,X,G,COV,TEMP,PP,IM)
C
C     If any elements of the PAR array are fixed then delete the
C     appropriate rows of X
C
      L2 = 0
      DO 120 L = 1, NPAR
         IF ( .NOT. PARHLD(L)) THEN
C
C           Replace the L2(th) row of X by the L(th) row of X
C
            L2 = L2 + 1
            DO 100 J = 1, MK2
               X(L2,J) = X(L,J)
  100       CONTINUE
         END IF
  120 CONTINUE
      NPARN = L2
C
C     set Y = I ** (DEL ** DEL)  where ** denotes kronecker product
C
      CALL G13DSY(K,M,DEL,COV,IM)
C
C     construct GITAL matrix
C
      IFAIL = 1
      CALL F02ABF(DEL,K,K,D,PP,K,WORK,IFAIL)
      DO 140 I = 1, K
         IF (D(I).LE.0.0D0) IFAIL = 1
  140 CONTINUE
      IF (IFAIL.GT.0) THEN
         IERR = 1
         RETURN
      END IF
C                         -(1/2)
C     GITAL = PP'   *    D      *     PP
C
      DO 160 A = 1, K
         D(A) = 1.0D0/SQRT(D(A))
  160 CONTINUE
C
      DO 220 I = 1, K
         DO 200 J = 1, I
            SUM = 0.0D0
            DO 180 A = 1, K
               SUM = SUM + PP(A,I)*PP(A,J)*D(A)
  180       CONTINUE
            GITAL(I,J) = SUM
            GITAL(J,I) = SUM
  200    CONTINUE
  220 CONTINUE
C
C     construct G matrix
C
      CALL G13DSY(K,M,GITAL,G,MK2)
C
C     construct the information matrix INFO = X'(GG')X
C
      DO 320 I = 1, NPARN
         DO 300 J = 1, I
            SUM = 0.0D0
            DO 280 L = 1, MK2
               SUM2 = 0.0D0
               DO 260 B = 1, MK2
                  SUM3 = 0.0D0
                  DO 240 A = 1, MK2
                     SUM3 = SUM3 + X(I,A)*G(A,B)
  240             CONTINUE
                  SUM2 = SUM2 + SUM3*G(L,B)
  260          CONTINUE
               SUM = SUM + SUM2*X(J,L)
  280       CONTINUE
            INFO(I,J) = SUM
            INFO(J,I) = SUM
  300    CONTINUE
  320 CONTINUE
C
C     invert INFO
C
      IFAIL = 1
      CALL F01ADF(NPARN,INFO,NPAR+1,IFAIL)
      IF (IFAIL.GT.0) THEN
         IERR = 1
         RETURN
      END IF
C
      DO 360 J = 1, NPARN
         DO 340 I = J, NPARN
            INFO(I,J) = INFO(I+1,J)
  340    CONTINUE
  360 CONTINUE
C
      DO 400 I = 2, NPARN
         DO 380 J = 1, I - 1
            INFO(J,I) = INFO(I,J)
  380    CONTINUE
  400 CONTINUE
C
C     now calculate (Y - X * INFO * X') / n
C
      DO 480 J = 1, MK2
         DO 460 I = J, MK2
            SUM = 0.0D0
            DO 440 B = 1, NPARN
               SUM2 = 0.0D0
               DO 420 A = 1, NPARN
                  SUM2 = SUM2 + X(A,I)*INFO(A,B)
  420          CONTINUE
               SUM = SUM + SUM2*X(B,J)
  440       CONTINUE
            COV(I,J) = (COV(I,J)-SUM)*TT
            COV(J,I) = COV(I,J)
  460    CONTINUE
  480 CONTINUE
      RETURN
C
      END
