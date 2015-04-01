      DOUBLE PRECISION FUNCTION G08CDZ(N1,N2,D,IERROR)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 D
      INTEGER                          IERROR, N1, N2
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, ATERM, EPS1, FAC, P, SR, TERM,
     *                                 TP, W, X, XJ, Z
      INTEGER                          I, J, M, N
C     .. Local Arrays ..
      DOUBLE PRECISION                 U(2501)
C     .. External Functions ..
      DOUBLE PRECISION                 G08CBZ, X02AMF
      EXTERNAL                         G08CBZ, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, INT, LOG, MAX, MIN,
     *                                 DBLE, SQRT
C     .. Executable Statements ..
C
      IERROR = 0
      IF (D.LT.X02AMF()) THEN
         P = 1.0D0
      ELSE
         M = MIN(N1,N2)
         N = MAX(N1,N2)
         IF (M*N.LE.10000 .AND. N.LE.2500) THEN
            X = DBLE(M*N)*D - 0.5D0
            U(1) = 1.0D0
            DO 20 J = 1, N
               U(J+1) = 1.0D0
               IF (DBLE(M*J).GT.X) U(J+1) = 0.0D0
   20       CONTINUE
            DO 60 I = 1, M
               W = DBLE(I)/DBLE(I+N)
               U(1) = W*U(1)
               IF (DBLE(N*I).GT.X) U(1) = 0.0D0
               DO 40 J = 1, N
                  U(J+1) = U(J) + U(J+1)*W
                  IF (DBLE(ABS(N*I-M*J)).GT.X) U(J+1) = 0.0D0
   40          CONTINUE
   60       CONTINUE
            P = U(N+1)
            P = 1.0D0 - P
            P = MIN(1.0D0,P)
            P = MAX(0.0D0,P)
         ELSE IF (M.LT.INT(N/10) .AND. M.LT.80) THEN
            Z = D
            IF (M.NE.1) Z = Z - 0.5D0/DBLE(N)
            Z = MAX(0.0D0,Z)
            TP = 2.0D0*G08CBZ(M,Z)
            P = MIN(1.0D0,TP)
         ELSE
            Z = SQRT(DBLE(M*N)/DBLE(M+N))*D + 0.5D0/SQRT(DBLE(N))
            A = -2.0D0*Z*Z
            IF (-A.LT.X02AMF()) THEN
               P = 1.0D0
            ELSE
               SR = SQRT(LOG(X02AMF())/A)
               FAC = 2.0D0
               P = 0.0D0
               EPS1 = 0.000005D0
               DO 80 J = 1, 500
                  XJ = DBLE(J)
                  IF (XJ.LT.SR) THEN
                     TERM = FAC*EXP(A*XJ*XJ)
                     P = P + TERM
                     ATERM = ABS(TERM)
                     IF (ATERM.LT.EPS1*P) THEN
                        P = MIN(1.0D0,P)
                        GO TO 100
                     ELSE
                        FAC = -FAC
                     END IF
                  ELSE
                     P = MIN(1.0D0,P)
                     GO TO 100
                  END IF
   80          CONTINUE
C
C              If the above method fails to converge in the allowed 500
C              iterations we return an IFAIL warning and a probability
C              = 1.0
C
               P = 1.0D0
               IERROR = 3
            END IF
         END IF
      END IF
  100 G08CDZ = P
C
      END
