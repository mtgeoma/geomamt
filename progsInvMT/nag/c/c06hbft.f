      SUBROUTINE C06HBF(M,N,X,INIT,TRIG,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14C REVISED. IER-872 (NOV 1990).
C     MARK 15 REVISED. IER-895 (APR 1991).
C
C     C06HBF computes multiple Fourier cosine transforms of sequences
C     of real data using the multiple real transform kernel C06FPX,
C     and pre- and post-processing steps described by Swarztrauber.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06HBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N), WORK(1:M,0:N-1), X(1:M,0:N)
C     .. Local Scalars ..
      DOUBLE PRECISION  PI, PIBYN, ROOT2, ROOTN, SUM
      INTEGER           I, IERROR, J, N2, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FPX
C     .. Intrinsic Functions ..
      INTRINSIC         COS, MOD, DBLE, SIN, SQRT
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0) THEN
         PIBYN = X01AAF(PI)/DBLE(N)
         ROOTN = SQRT(DBLE(N))
         ROOT2 = SQRT(2.0D0)
C
C        Calculate the auxiliary array D and store it in the WORK array
C
C        Set the first column of the array to be transformed to the
C        required value computed from the input array
C
         DO 20 I = 1, M
            WORK(I,0) = 0.5D0*(X(I,0)+X(I,N))
   20    CONTINUE
C
C        Set up the rest of the D array
C
         DO 60 J = 1, N - 1
            DO 40 I = 1, M
               WORK(I,J) = 0.5D0*(X(I,J)+X(I,N-J)) - SIN(DBLE(J)*PIBYN)
     *                     *(X(I,J)-X(I,N-J))
   40       CONTINUE
   60    CONTINUE
C
C        The first element of the transform is generated here and
C        stored in X(1:M,N), which is not used as workspace by C06FPX
C
         DO 100 I = 1, M
            SUM = 0.5D0*X(I,0)
            DO 80 J = 1, N - 1
               SUM = SUM + X(I,J)*COS(DBLE(J)*PIBYN)
   80       CONTINUE
            SUM = SUM - 0.5D0*X(I,N)
            X(I,N) = ROOT2*SUM/ROOTN
  100    CONTINUE
C
C        Transform the D array, stored in WORK, using X as the workspace
C        array
C
         CALL C06FPX(WORK,X,M,N,Q,NQ,TRIG)
C
C        Now extract the cosine transform from the Fourier transform in
C        WORK, putting the required transform in X
C
         DO 120 I = 1, M
            X(I,0) = ROOT2*WORK(I,0)
  120    CONTINUE
         DO 140 I = 1, M
            X(I,1) = X(I,N)
  140    CONTINUE
C
         DO 180 J = 1, (N-1)/2
            DO 160 I = 1, M
               X(I,2*J) = ROOT2*WORK(I,J)
               X(I,2*J+1) = X(I,2*J-1) - ROOT2*WORK(I,N-J)
  160       CONTINUE
  180    CONTINUE
         IF (MOD(N,2).EQ.0) THEN
            N2 = N/2
            DO 200 I = 1, M
               X(I,N) = ROOT2*WORK(I,N2)
  200       CONTINUE
         END IF
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,' but TRIG array never initialised')
99995 FORMAT (' ** INIT = ',A1,' but N and TRIG array incompatible')
      END
