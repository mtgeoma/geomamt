      SUBROUTINE E04NAV(N,NROWH,NCOLH,CVEC,HESS,QPHESS,WRK,HX)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C *********************************************************************
C     E04NAV  PRINTS QUANTITIES DEFINING THE QUADRATIC FUNCTION GIVEN
C     TO  E04NAX.
C
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C     VERSION OF APRIL 1982.
C *********************************************************************
C
C     QPHESS
C     .. Scalar Arguments ..
      INTEGER           N, NCOLH, NROWH
C     .. Array Arguments ..
      DOUBLE PRECISION  CVEC(N), HESS(NROWH,NCOLH), HX(N), WRK(N)
C     .. Subroutine Arguments ..
      EXTERNAL          QPHESS
C     .. Scalars in Common ..
      INTEGER           ISTART, MSG, NOUT
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, ZERO
      INTEGER           I, J, K, NHESS
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Subroutines ..
      EXTERNAL          F06FBF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AE04VC/NOUT, MSG, ISTART
C     .. Data statements ..
      DATA              ZERO/0.0D+0/, ONE/1.0D+0/
C     .. Executable Statements ..
C
      CALL X04BAF(NOUT,' ')
      CALL X04BAF(NOUT,' ')
      CALL X04BAF(NOUT,' ')
      CALL X04BAF(NOUT,' ')
      WRITE (REC,FMT=99999)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      CALL X04BAF(NOUT,REC(3))
      WRITE (REC,FMT=99998)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      DO 20 I = 1, N, 5
         WRITE (REC,FMT=99997) (CVEC(J),J=I,MIN(I+4,N))
         CALL X04BAF(NOUT,REC(1))
   20 CONTINUE
C
C     PRINT  HESS  UNLESS IT APPEARS TO BE IMPLICIT.
C
      WRITE (REC,FMT=99996) NROWH, NCOLH
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      IF (NROWH.EQ.1 .AND. NCOLH.EQ.1) GO TO 100
      IF (NCOLH.EQ.1) THEN
         WRITE (REC,FMT=99995)
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 40 I = 1, NROWH, 5
            WRITE (REC,FMT=99997) (HESS(J,1),J=I,MIN(NROWH,I+4))
            CALL X04BAF(NOUT,REC(1))
   40    CONTINUE
      END IF
      IF (NCOLH.EQ.1) GO TO 100
      NHESS = MIN(NCOLH,N)
      DO 80 J = 1, NHESS
         WRITE (REC,FMT=99994) J
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 60 I = 1, NHESS, 5
            WRITE (REC,FMT=99997) (HESS(K,J),K=I,MIN(NHESS,I+4))
            CALL X04BAF(NOUT,REC(1))
   60    CONTINUE
   80 CONTINUE
C
C     CALL  QPHESS  TO COMPUTE EACH COLUMN OF THE HESSIAN.
C
  100 WRITE (REC,FMT=99993)
      CALL X04BAF(NOUT,REC(1))
      CALL X04BAF(NOUT,REC(2))
      CALL X04BAF(NOUT,REC(3))
      CALL F06FBF(N,ZERO,WRK,1)
      DO 140 J = 1, N
         WRK(J) = ONE
         CALL QPHESS(N,NROWH,NCOLH,J,HESS,WRK,HX)
         WRITE (REC,FMT=99992) J
         CALL X04BAF(NOUT,REC(1))
         CALL X04BAF(NOUT,REC(2))
         DO 120 I = 1, N
            WRITE (REC,FMT=99997) (HX(K),K=I,MIN(N,I+4))
            CALL X04BAF(NOUT,REC(1))
  120    CONTINUE
         WRK(J) = ZERO
  140 CONTINUE
      RETURN
C
C
C     END OF E04NAV  ( QPDUMP )
99999 FORMAT (' ',/' OUTPUT FROM E04NAV',/' ******************')
99998 FORMAT (/' CVEC ...')
99997 FORMAT (5G15.6)
99996 FORMAT (/' NROWH =',I6,'    NCOLH =',I6)
99995 FORMAT (/' HESS ...')
99994 FORMAT (/' COLUMN',I6,'  OF  HESS ...')
99993 FORMAT (//' THE FOLLOWING IS RETURNED BY  QPHESS.')
99992 FORMAT (/' COLUMN',I6,'  FROM  QPHESS ...')
      END
