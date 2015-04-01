      SUBROUTINE E04NAW(UNITQ,QPHESS,N,NCOLR,NCOLZ,NCTOTL,NFREE,NHESS,
     *                  NQ,NROWH,NCOLH,NROWRT,NCOLRT,KFREE,HSIZE,HESS,
     *                  RT,SCALE,ZY,HZ,WRK)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C
C *********************************************************************
C     E04NAW  COMPUTES THE CHOLESKY FACTOR  R  OF THE PROJECTED HESSIAN
C     Z(T) H Z,  GIVEN  Z  AND ITS DIMENSIONS  NFREE BY NCOLZ.
C     IF THE PROJECTED HESSIAN IS INDEFINITE, A SMALLER CHOLESKY
C     FACTORIZATION  R1(T) R1 = Z1(T) H Z1  IS RETURNED, WHERE  Z1  IS
C     COMPOSED OF  NCOLR  COLUMNS OF  Z.  COLUMN INTERCHANGES ARE
C     USED TO MAXIMIZE  NCOLR.  THESE ARE APPLIED TO  Z.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     ORIGINAL VERSION OF JANUARY 1983.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HSIZE
      INTEGER           N, NCOLH, NCOLR, NCOLRT, NCOLZ, NCTOTL, NFREE,
     *                  NHESS, NQ, NROWH, NROWRT
      LOGICAL           UNITQ
C     .. Array Arguments ..
      DOUBLE PRECISION  HESS(NROWH,NCOLH), HZ(N), RT(NROWRT,NCOLRT),
     *                  SCALE(NCTOTL), WRK(N), ZY(NQ,NQ)
      INTEGER           KFREE(N)
C     .. Subroutine Arguments ..
      EXTERNAL          QPHESS
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
      LOGICAL           SCLDQP
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  D, DMAX, DMIN, EPSMCH, ONE, T, ZERO
      INTEGER           I, J, J1, JTHCOL, K, KMAX, KSAVE, LEN, NUM
C     .. Local Arrays ..
      CHARACTER*60      REC(3)
C     .. External Subroutines ..
      EXTERNAL          E04VDN, F06FBF, F06FCF, DAXPY, DCOPY, DSCAL,
     *                  X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
      COMMON            /AX02ZA/WMACH
      COMMON            /JE04VC/SCLDQP
C     .. Save statement ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
C     .. Executable Statements ..
      EPSMCH = WMACH(3)
C
      NCOLR = 0
      IF (NCOLZ.EQ.0) GO TO 300
C
C ---------------------------------------------------------------------
C     COMPUTE  Z(T) H Z  AND STORE THE UPPER-TRIANGULAR SYMMETRIC PART
C     IN THE FIRST  NCOLZ  COLUMNS OF  RT.
C ---------------------------------------------------------------------
      DO 80 K = 1, NCOLZ
         CALL F06FBF(N,ZERO,WRK,1)
         IF (UNITQ) GO TO 40
C
C        EXPAND THE COLUMN OF  Z  INTO AN  N-VECTOR.
C
         DO 20 I = 1, NFREE
            J = KFREE(I)
            WRK(J) = ZY(I,K)
   20    CONTINUE
         IF (SCLDQP) CALL F06FCF(N,SCALE,1,WRK,1)
         JTHCOL = 0
         GO TO 60
C
C        ONLY BOUNDS ARE IN THE WORKING SET.  THE  K-TH COLUMN OF  Z  IS
C        JUST A COLUMN OF THE IDENTITY MATRIX.
C
   40    JTHCOL = KFREE(K)
         WRK(JTHCOL) = ONE
C
C        SET  RT(*,K)  =  TOP OF   H * (COLUMN OF  Z).
C
   60    CALL QPHESS(N,NROWH,NCOLH,JTHCOL,HESS,WRK,HZ)
         NHESS = NHESS + 1
C
         IF (UNITQ .AND. SCLDQP) CALL DSCAL(N,SCALE(JTHCOL),HZ,1)
         IF (SCLDQP) CALL F06FCF(N,SCALE,1,HZ,1)
C
         CALL E04VDN(4,N,NFREE,NCOLZ,NFREE,NQ,UNITQ,KFREE,KFREE,HZ,ZY,
     *               WRK)
C
         CALL DCOPY(NCOLZ,HZ,1,RT(1,K),1)
C
C        UPDATE AN ESTIMATE OF THE SIZE OF THE PROJECTED HESSIAN.
C
         HSIZE = MAX(HSIZE,ABS(RT(K,K)))
   80 CONTINUE
C
C ---------------------------------------------------------------------
C     FORM THE CHOLESKY FACTORIZATION  R(T) R  =  Z(T) H Z  AS FAR AS
C     POSSIBLE, USING SYMMETRIC ROW AND COLUMN INTERCHANGES.
C ---------------------------------------------------------------------
      DMIN = EPSMCH*HSIZE
C
      DO 260 J = 1, NCOLZ
C
C        FIND THE MAXIMUM REMAINING DIAGONAL.
C
         KMAX = J
         DMAX = RT(J,J)
         DO 100 K = J, NCOLZ
            D = RT(K,K)
            IF (DMAX.GE.D) GO TO 100
            DMAX = D
            KMAX = K
  100    CONTINUE
C
C        SEE IF THE DIAGONAL IS BIG ENOUGH.
C
         IF (DMAX.LE.DMIN) GO TO 280
         NCOLR = J
C
C        PERMUTE THE COLUMNS OF  Z.
C
         IF (KMAX.EQ.J) GO TO 220
         IF (UNITQ) GO TO 120
         CALL DCOPY(NFREE,ZY(1,KMAX),1,WRK,1)
         CALL DCOPY(NFREE,ZY(1,J),1,ZY(1,KMAX),1)
         CALL DCOPY(NFREE,WRK,1,ZY(1,J),1)
         GO TO 140
C
C        Z  IS NOT STORED EXPLICITLY.
C
  120    KSAVE = KFREE(KMAX)
         KFREE(KMAX) = KFREE(J)
         KFREE(J) = KSAVE
C
C        INTERCHANGE ROWS AND COLUMNS OF THE PROJECTED HESSIAN.
C
  140    DO 160 I = 1, J
            T = RT(I,KMAX)
            RT(I,KMAX) = RT(I,J)
            RT(I,J) = T
  160    CONTINUE
C
         DO 180 K = J, KMAX
            T = RT(K,KMAX)
            RT(K,KMAX) = RT(J,K)
            RT(J,K) = T
  180    CONTINUE
C
         DO 200 K = KMAX, NCOLZ
            T = RT(KMAX,K)
            RT(KMAX,K) = RT(J,K)
            RT(J,K) = T
  200    CONTINUE
C
         RT(KMAX,KMAX) = RT(J,J)
C
C        SET THE DIAGONAL ELEMENT OF  R.
C
  220    D = SQRT(DMAX)
         RT(J,J) = D
         IF (J.EQ.NCOLZ) GO TO 260
C
C        SET THE ABOVE-DIAGONAL ELEMENTS OF THE K-TH ROW OF  R,
C        AND UPDATE THE ELEMENTS OF ALL REMAINING ROWS.
C
         J1 = J + 1
         DO 240 K = J1, NCOLZ
            T = RT(J,K)/D
            RT(J,K) = T
C
C           R(I,K)  =  R(I,K)  - T * R(J,I),   I = J1, K.
C
            NUM = K - J
            LEN = NROWRT*(NUM-1) + 1
            IF (T.NE.ZERO) CALL DAXPY(NUM,(-T),RT(J,J1),NROWRT,RT(J1,K),
     *                                1)
  240    CONTINUE
  260 CONTINUE
C
  280 IF (NCOLR.EQ.NCOLZ) GO TO 300
      IF (MSG.GE.80) THEN
         WRITE (REC,FMT=99999) NCOLR, NCOLZ
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         CALL X04BAF(NOUT,REC(3))
      END IF
C
  300 RETURN
C
C
C     END OF E04NAW (QPCRSH)
99999 FORMAT (/' //E04NAW//  INDEFINITE PROJECTED HESSIAN.',/' //E04NA',
     *  'W//  NCOLR =',I5,6X,'NCOLZ =',I5)
      END
