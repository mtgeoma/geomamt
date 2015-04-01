      SUBROUTINE C06GSF(M,N,X,U,V,IFAIL)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06GSF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  U(0:M*N-1), V(0:M*N-1), X(0:M*N-1)
C     .. Local Scalars ..
      INTEGER           I, IB, IERROR, IF, IM, J, NJ, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      IF (M.LT.1) GO TO 100
      IF (N.LT.1) GO TO 120
C
      IERROR = 0
      DO 20 I = 0, M - 1
         U(I) = X(I)
         V(I) = 0.0D0
   20 CONTINUE
C
      DO 60 J = 1, (N-1)/2
         DO 40 I = 0, M - 1
            NJ = N - J
            IF = J*M + I
            IB = NJ*M + I
            U(IF) = X(IF)
            V(IF) = X(IB)
            U(IB) = X(IF)
            V(IB) = -V(IF)
   40    CONTINUE
   60 CONTINUE
C
      IF (MOD(N,2).EQ.0) THEN
         DO 80 I = 0, M - 1
            IM = (N/2)*M + I
            U(IM) = X(IM)
            V(IM) = 0.0D0
   80    CONTINUE
      END IF
      GO TO 140
C
  100 IERROR = 1
      WRITE (REC(1),FMT=99999) M
      GO TO 140
  120 IERROR = 2
      WRITE (REC(1),FMT=99998) N
C
  140 NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
      END
