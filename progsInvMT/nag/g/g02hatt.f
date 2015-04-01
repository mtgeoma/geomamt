      SUBROUTINE G02HAT(ISIGMA,INDW,IPSI,N,D,WGT,BETA,MAXIT,TOL,WORK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Calculate value of BETA Huber's CHI function.
C
C     .. Parameters ..
      DOUBLE PRECISION  Q75
      PARAMETER         (Q75=0.6744897501962755D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BETA, D, TOL
      INTEGER           IFAIL, INDW, IPSI, ISIGMA, MAXIT, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WGT(N), WORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AMAXEX, ANORMC, B, D2, DC, DW, DW2, FUNC, FUNCP,
     *                  PC, PDFX, PDFX2, SMF, SMFP, UP, W2, XN
      INTEGER           I, IBIT, IFAIL2
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X01AAF, X02AKF
      EXTERNAL          S15ABF, X01AAF, X02AKF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, DBLE, SQRT
C     .. Executable Statements ..
      AMAXEX = -LOG(X02AKF())
      ANORMC = SQRT(2.0D0*X01AAF(0.0D0))
      XN = DBLE(N)
      BETA = 0.0D0
      D2 = D*D
      IF (ISIGMA.GT.0) THEN
C
C        CASE ISIGMA = POSITIVE
C
         IF (IPSI.EQ.0) THEN
            BETA = .5D0
         ELSE
            BETA = 0.0D0
            IF (INDW.GT.0) THEN
C
C              SCHWEPPE CASE
C
               IFAIL2 = 0
               DO 20 I = 1, N
                  W2 = WGT(I)*WGT(I)
                  DW = WGT(I)*D
                  PC = S15ABF(DW,IFAIL2)
                  DW2 = DW*DW
                  DC = 0.0D0
                  IF (DW2.LT.AMAXEX) DC = EXP(-DW2/2.0D0)/ANORMC
                  B = (-DW*DC+PC-0.5D0)/W2 + (1.0D0-PC)*D2
                  BETA = B*W2/XN + BETA
   20          CONTINUE
            ELSE
               IF (INDW.LT.0) THEN
C
C                 MALLOWS CASE
C
                  DO 40 I = 1, N
                     BETA = WGT(I) + BETA
   40             CONTINUE
                  BETA = BETA/XN
               ELSE
                  BETA = 1.0D0
               END IF
               IFAIL2 = 0
               PC = S15ABF(D,IFAIL2)
               DC = 0.0D0
               IF (D2.LT.AMAXEX) DC = EXP(-D2/2.0D0)/ANORMC
               BETA = (-D*DC+PC-0.5D0+(1.0D0-PC)*D2)*BETA
            END IF
         END IF
      END IF
      IF (ISIGMA.LT.0) THEN
C
C        CASE ISIGMA= NEGATIVE
C
         BETA = Q75
         IF (INDW.LT.0) THEN
C
C           MALLOWS CASE
C
            DO 60 I = 1, N
               WORK(I) = SQRT(WGT(I))
   60       CONTINUE
C
C           START ITERATIVE METHOD
C
            IBIT = 0
   80       IBIT = IBIT + 1
            SMF = 0.0D0
            IFAIL2 = 0
            SMFP = 0.0D0
            DO 100 I = 1, N
               IF (WORK(I).GT.0.0D0) THEN
                  PDFX = BETA/WORK(I)
                  PC = S15ABF(PDFX,IFAIL2)
                  PDFX2 = PDFX*PDFX
                  DC = 0.0D0
                  IF (PDFX2.LT.AMAXEX) DC = EXP(-PDFX2/2.0D0)/ANORMC
                  SMF = SMF + PC
                  SMFP = DC/WORK(I) + SMFP
               END IF
  100       CONTINUE
            FUNC = SMF/XN - 0.75D0
            FUNCP = SMFP/XN
            UP = FUNC/FUNCP
            BETA = BETA - UP
C
C           CHECK FOR CONVERGENCE
C
            IF (UP.GT.TOL) THEN
               IF (IBIT.LT.MAXIT) THEN
                  GO TO 80
               ELSE
                  IFAIL = 1
               END IF
            END IF
         END IF
      END IF
      RETURN
C
      END
