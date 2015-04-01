      SUBROUTINE F04YCF(ICASE,N,X,ESTNRM,WORK,IWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15 REVISED. IER-913 (APR 1991).
C
C     F04YCF estimates the 1-norm of a square, real matrix A.
C     Reverse communication is used for evaluating matrix-vector
C     products.
C
C     On entry,
C
C        N       integer.  The order of the matrix.  N .ge. 1.
C
C        IWORK   integer(N). Used as workspace.
C
C        ICASE   integer. ICASE must be set to zero on initial entry.
C
C     On intermediate returns:
C
C        ICASE   = 1 or 2.
C
C        X       real(N), must be overwritten by
C
C                     A*X,             if ICASE=1,
C                     transpose(A)*X,  if ICASE=2,
C
C                and F04YCF must be re-called, with all the other
C                parameters unchanged.
C
C     On final return,
C
C        ICASE   = 0.
C
C        ESTNRM  real. Contains an estimate (a lower bound) for norm(A).
C
C        WORK    real(N). Contains vector V such that V = A*W, where
C                ESTNRM = norm(V)/norm(W) (W  is not returned).
C
C     Nick Higham, University of Manchester.
C
C     Reference:
C     N.J. Higham (1987) Fortran Codes for Estimating
C     the One-norm of a Real or Complex Matrix, with Applications
C     to Condition  Estimation, Numerical Analysis Report No. 135,
C     University of Manchester, Manchester M13 9PL, England.
C
C     .. Parameters ..
      INTEGER           ITMAX
      PARAMETER         (ITMAX=5)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04YCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ESTNRM
      INTEGER           ICASE, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(N), X(N)
      INTEGER           IWORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALTSGN, ESTOLD, TEMP
      INTEGER           I, IER, ITER, J, JLAST, JUMP, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DASUM
      INTEGER           IDAMAX, P01ABF
      EXTERNAL          DASUM, IDAMAX, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, NINT, DBLE, SIGN
C     .. Save statement ..
      SAVE              ITER, J, JUMP
C     .. Executable Statements ..
      IF (N.LT.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         ICASE = 0
         GO TO 340
      END IF
      IER = 0
      NREC = 0
      IF (ICASE.EQ.0) THEN
         DO 20 I = 1, N
            X(I) = 1.0D0/DBLE(N)
   20    CONTINUE
         ICASE = 1
         JUMP = 1
         GO TO 340
      END IF
C
      GO TO (40,80,140,220,300) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
C
   40 CONTINUE
      IF (N.EQ.1) THEN
         WORK(1) = X(1)
         ESTNRM = ABS(WORK(1))
C        ... QUIT
         GO TO 320
      END IF
      ESTNRM = DASUM(N,X,1)
C
      DO 60 I = 1, N
         X(I) = SIGN(1.0D0,X(I))
         IWORK(I) = NINT(X(I))
   60 CONTINUE
      ICASE = 2
      JUMP = 2
      GO TO 340
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
   80 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  100 CONTINUE
      DO 120 I = 1, N
         X(I) = 0.0D0
  120 CONTINUE
      X(J) = 1.0D0
      ICASE = 1
      JUMP = 3
      GO TO 340
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  140 CONTINUE
      CALL DCOPY(N,X,1,WORK,1)
      ESTOLD = ESTNRM
      ESTNRM = DASUM(N,WORK,1)
      DO 160 I = 1, N
         IF (NINT(SIGN(1.0D0,X(I))).NE.IWORK(I)) GO TO 180
  160 CONTINUE
C
C     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 240
C
  180 CONTINUE
C     TEST FOR CYCLING.
      IF (ESTNRM.LE.ESTOLD) GO TO 260
      DO 200 I = 1, N
         X(I) = SIGN(1.0D0,X(I))
         IWORK(I) = NINT(X(I))
  200 CONTINUE
      ICASE = 2
      JUMP = 4
      GO TO 340
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
  220 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
         ITER = ITER + 1
         GO TO 100
      END IF
C
C     ITERATION COMPLETE.  FINAL STAGE.
C
  240 CONTINUE
C
  260 ALTSGN = 1.0D0
      DO 280 I = 1, N
         X(I) = ALTSGN*(1+DBLE(I-1)/DBLE(N-1))
         ALTSGN = -ALTSGN
  280 CONTINUE
      ICASE = 1
      JUMP = 5
      GO TO 340
C
C     ................ ENTRY   (JUMP = 5)
C
  300 CONTINUE
      TEMP = 2.0D0*DASUM(N,X,1)/DBLE(3*N)
      IF (TEMP.GT.ESTNRM) THEN
         CALL DCOPY(N,X,1,WORK,1)
         ESTNRM = TEMP
      END IF
C
  320 ICASE = 0
  340 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N .lt. 1: N =',I16)
      END
