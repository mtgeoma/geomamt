      SUBROUTINE F04ZCF(ICASE,N,X,ESTNRM,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15 REVISED. IER-914 (APR 1991).
C
C     F04ZCF estimates the 1-norm of a square, complex matrix  A.
C     Reverse communication is used for evaluating matrix-vector
C     products.
C
C     On entry,
C
C        N       integer. The order of the matrix.  N .ge. 1.
C
C       ICASE    integer. ICASE must be set to zero on initial entry.
C
C     On intermediate returns
C
C       ICASE    = 1 or 2.
C
C        X       complex(N). X must be overwritten by
C
C                     A*X,          if ICASE=1,
C                     conjg(A')*X,  if ICASE=2,
C
C                and F04ZCF must be re-called, with all the other
C                parameters unchanged.
C
C     On final return,
C
C        ICASE    = 0.
C
C        ESTNRM  real. Contains an estimate (a lower bound) for norm(A).
C
C        WORK    complex(N). Contains vector V such that V = A*W, where
C                ESTNRM = norm(V)/norm(W) (W is not returned).
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
      PARAMETER         (SRNAME='F04ZCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ESTNRM
      INTEGER           ICASE, IFAIL, N
C     .. Array Arguments ..
      COMPLEX*16        WORK(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALTSGN, ESTOLD, SMAX, TEMP
      INTEGER           I, IER, ITER, J, JLAST, JUMP, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX, DBLE
C     .. Save statement ..
      SAVE              ITER, J, JUMP
C     .. Executable Statements ..
      IF (N.LT.1) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT=99999) N
         ICASE = 0
         GO TO 380
      END IF
      IER = 0
      NREC = 0
      IF (ICASE.EQ.0) THEN
         DO 20 I = 1, N
            X(I) = DCMPLX(1.0D0/DBLE(N),0.0D0)
   20    CONTINUE
         ICASE = 1
         JUMP = 1
         GO TO 380
      END IF
C
      GO TO (40,100,180,240,320) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
C
   40 CONTINUE
      IF (N.EQ.1) THEN
         WORK(1) = X(1)
         ESTNRM = ABS(WORK(1))
C        ... QUIT
         GO TO 360
      END IF
      ESTNRM = 0.0D0
      DO 60 I = 1, N
         ESTNRM = ESTNRM + ABS(X(I))
   60 CONTINUE
C
      DO 80 I = 1, N
         IF (ABS(X(I)).NE.0.0D0) THEN
            X(I) = X(I)/DCMPLX(ABS(X(I)),0.0D0)
         ELSE
            X(I) = (1.0D0,0.0D0)
         END IF
   80 CONTINUE
      ICASE = 2
      JUMP = 2
      GO TO 380
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
C
  100 CONTINUE
      J = 1
      SMAX = ABS(DBLE(X(1)))
      DO 120 I = 2, N
         IF (ABS(DBLE(X(I))).GT.SMAX) THEN
            SMAX = ABS(DBLE(X(I)))
            J = I
         END IF
  120 CONTINUE
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  140 CONTINUE
      DO 160 I = 1, N
         X(I) = (0.0D0,0.0D0)
  160 CONTINUE
      X(J) = (1.0D0,0.0D0)
      ICASE = 1
      JUMP = 3
      GO TO 380
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  180 CONTINUE
      CALL ZCOPY(N,X,1,WORK,1)
      ESTOLD = ESTNRM
      ESTNRM = 0.0D0
      DO 200 I = 1, N
         ESTNRM = ESTNRM + ABS(WORK(I))
  200 CONTINUE
C     TEST FOR CYCLING.
      IF (ESTNRM.LE.ESTOLD) GO TO 280
      DO 220 I = 1, N
         IF (ABS(X(I)).NE.0.0D0) THEN
            X(I) = X(I)/DCMPLX(ABS(X(I)),0.0D0)
         ELSE
            X(I) = (1.0D0,0.0D0)
         END IF
  220 CONTINUE
      ICASE = 2
      JUMP = 4
      GO TO 380
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
C
  240 CONTINUE
      JLAST = J
      J = 1
      SMAX = ABS(DBLE(X(1)))
      DO 260 I = 2, N
         IF (ABS(DBLE(X(I))).GT.SMAX) THEN
            SMAX = ABS(DBLE(X(I)))
            J = I
         END IF
  260 CONTINUE
      IF ((ABS(DBLE(X(JLAST))).NE.ABS(DBLE(X(J)))) .AND. (ITER.LT.ITMAX)
     *    ) THEN
         ITER = ITER + 1
         GO TO 140
      END IF
C
C     ITERATION COMPLETE.  FINAL STAGE.
C
  280 ALTSGN = 1.0D0
      DO 300 I = 1, N
         X(I) = ALTSGN*(1.0D0+(I-1.0D0)/(N-1.0D0))
         ALTSGN = -ALTSGN
  300 CONTINUE
      ICASE = 1
      JUMP = 5
      GO TO 380
C
C     ................ ENTRY   (JUMP = 5)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  320 CONTINUE
      TEMP = 0.0D0
      DO 340 I = 1, N
         TEMP = TEMP + ABS(X(I))
  340 CONTINUE
      TEMP = 2*TEMP/(3*N)
      IF (TEMP.GT.ESTNRM) THEN
         CALL ZCOPY(N,X,1,WORK,1)
         ESTNRM = TEMP
      END IF
C
  360 ICASE = 0
  380 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, N .lt. 1: N =',I16)
      END
