      SUBROUTINE G13EBZ(N,M,P,S,LDS,A,LDA,B,LDB,Q,LDQ,C,LDC,R,LDR,K,LDK,
     *                  H,LDH,IWORK,RWORK,TOLER,MULTBQ,IERR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Based on SLICOT routine FB01FD
C

C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOLER
      INTEGER           IERR, LDA, LDB, LDC, LDH, LDK, LDQ, LDR, LDS, M,
     *                  N, P
      LOGICAL           MULTBQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*), H(LDH,*),
     *                  K(LDK,*), Q(LDQ,*), R(LDR,*), RWORK(N+P,*),
     *                  S(LDS,*)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DZ1, RCOND
      INTEGER           I, I1, IN, IND, INFO, IPM, J, LDW, MINMP, MINPN,
     *                  NAXPY, P1, PI, PI1, PJ, PM1, PN, PNM
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DTRCON, DTRSV, F06FBF,
     *                  F06FSF, F06FUF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IERR = 0
C
C     Construction of the pre-array RWORK.
C
      MINPN = MIN(P,N)
      P1 = P + 1
      PM1 = P1 + M
      PN = P + N
      PNM = PN + M
      DO 20 J = 1, PNM
         CALL F06FBF(PN,ZERO,RWORK(1,J),1)
   20 CONTINUE
C
C     First part -  Storing lower triangular factor R in the (1,1) block
C                   of RWORK.
C
      DO 40 I = 1, P
         CALL DCOPY(P-I+1,R(I,I),1,RWORK(I,I),1)
   40 CONTINUE
C
C     Second part - Storing B x Q in the (2,2) block of RWORK.
C
      IF (MULTBQ) THEN
         DO 60 I = 1, M
            CALL DCOPY(N,B(1,I),1,RWORK(P1,P+I),1)
   60    CONTINUE
      ELSE
         DO 80 I = 1, M
            CALL DGEMV('N',N,M-I+1,ONE,B(1,I),LDB,Q(I,I),1,ZERO,
     *                 RWORK(P1,P+I),1)
   80    CONTINUE
      END IF
C
C     Third part -  Storing C x S in the (1,3) block of RWORK.
C
      DO 120 I = 1, MINPN
         IPM = I + P + M
         NAXPY = P - I + 1
         DO 100 J = I, MINPN
            CALL DAXPY(NAXPY,S(J,I),C(J,J),1,RWORK(J,IPM),1)
            NAXPY = NAXPY - 1
  100    CONTINUE
  120 CONTINUE
C
C     Fourth part  -  Storing A x S in the (2,3) block of RWORK.
C
      DO 160 I = 1, N
         IPM = I + P + M
         NAXPY = N
         IND = 1
         DO 140 J = I, N
            IF (J.GT.P1) THEN
               IND = IND + 1
               NAXPY = NAXPY - 1
            END IF
            CALL DAXPY(NAXPY,S(J,I),A(IND,J),1,RWORK(IND+P,IPM),1)
  140    CONTINUE
  160 CONTINUE
C
C     Triangularization (2 steps).
C
C     Step 1: eliminate the (1,3) block RWORK.
C
      LDW = N + P
      DO 200 I = 1, P
         IN = MIN(I,N)
         CALL F06FSF(IN,RWORK(I,I),RWORK(I,PM1),LDW,TOLER,DZ1)
         I1 = I + 1
         DO 180 J = I1, PN
            CALL F06FUF(IN,RWORK(I,PM1),LDW,DZ1,RWORK(J,I),RWORK(J,PM1),
     *                  LDW)
  180    CONTINUE
  200 CONTINUE
C
C     Step 2: triangularize the remaining (2,2) and (2,3) blocks of
C             RWORK.
C
      DO 240 I = 1, N
         MINMP = MIN(M+P,M+N-I)
         PI = P + I
         PI1 = PI + 1
         CALL F06FSF(MINMP,RWORK(PI,PI),RWORK(PI,PI1),LDW,TOLER,DZ1)
         IF (PI1.LE.PN) THEN
            DO 220 J = PI1, PN
               CALL F06FUF(MINMP,RWORK(PI,PI1),LDW,DZ1,RWORK(J,PI),
     *                     RWORK(J,PI1),LDW)
  220       CONTINUE
         END IF
  240 CONTINUE
C
C     Copy S; copy G into K
C
      DO 260 J = 1, N
         PJ = P + J
         CALL DCOPY(N-J+1,RWORK(PJ,PJ),1,S(J,J),1)
         CALL DCOPY(P,RWORK(P+J,1),LDW,K(J,1),LDK)
  260 CONTINUE
C
C        Copy H
C
      DO 280 J = 1, P
         CALL DCOPY(P-J+1,RWORK(J,J),1,H(J,J),1)
  280 CONTINUE
C
C     Calculate condition number for H
C
      CALL DTRCON('1-norm','Lower','N',P,H,LDH,RCOND,RWORK,IWORK,INFO)
      IF (RCOND.LT.TOLER) THEN
         IERR = 2
      ELSE
C
C     Calculate K
C
         DO 300 J = 1, N
            CALL DTRSV('L','T','N',P,H,LDH,K(J,1),LDK)
  300    CONTINUE
      END IF
C
      RETURN
      END
