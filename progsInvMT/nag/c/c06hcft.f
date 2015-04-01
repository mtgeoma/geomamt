      SUBROUTINE C06HCF(DIRECT,M,N,X,INIT,TRIG,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 13B REVISED. IER-650 (AUG 1988).
C     MARK 14C REVISED. IER-873 (NOV 1990).
C     MARK 15 REVISED. IER-896 (APR 1991).
C
C     C06HCF computes multiple quarter-wave Fourier sine transforms
C     of sequences of real data using the multiple real and
C     Hermitian transform kernels C06FPX and C06FQX, and pre- and
C     post-processing steps described by Swarztrauber.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06HCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       DIRECT, INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N), WORK(1:M,0:N-1), X(1:M,1:N)
C     .. Local Scalars ..
      DOUBLE PRECISION  PI, PIBY2N, ROOT2
      INTEGER           I, IERROR, J, JFAIL, K, N2, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FPX, C06FQX, C06GQF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, MOD, DBLE, SIN, SQRT
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      PIBY2N = X01AAF(PI)/DBLE(2*N)
      IF (DIRECT.NE.'F' .AND. DIRECT.NE.'f' .AND. DIRECT.NE.'B' .AND.
     *    DIRECT.NE.'b') IERROR = 6
      IF (IERROR.EQ.0) THEN
         IF (DIRECT.EQ.'F' .OR. DIRECT.EQ.'f') THEN
            ROOT2 = SQRT(2.0D0)
C
C           Calculate the auxiliary array D and store it in
C           the WORK array
C
            DO 20 I = 1, M
               WORK(I,0) = X(I,N)
   20       CONTINUE
            DO 60 J = 1, N - 1
               DO 40 I = 1, M
                  WORK(I,J) = COS(DBLE(J)*PIBY2N)*(X(I,J)+X(I,N-J)) +
     *                        SIN(DBLE(J)*PIBY2N)*(X(I,J)-X(I,N-J))
   40          CONTINUE
   60       CONTINUE
C
C           Transform the D array, stored in WORK, using X as the
C           workspace array
C
            CALL C06FPX(WORK,X,M,N,Q,NQ,TRIG)
C
C           Now extract the quarter-wave sine transform from the Fourier
C           transform in WORK, putting the required transform in X
C
            DO 80 I = 1, M
               X(I,1) = 0.5D0*WORK(I,0)
   80       CONTINUE
            DO 120 J = 1, (N-1)/2
               DO 100 I = 1, M
                  X(I,2*J) = -0.5D0*(WORK(I,J)+WORK(I,N-J))
                  X(I,2*J+1) = 0.5D0*(WORK(I,J)-WORK(I,N-J))
  100          CONTINUE
  120       CONTINUE
            IF (MOD(N,2).EQ.0) THEN
               N2 = N/2
               DO 140 I = 1, M
                  X(I,N) = -0.5D0*WORK(I,N2)
  140          CONTINUE
            END IF
         ELSE IF (DIRECT.EQ.'B' .OR. DIRECT.EQ.'b') THEN
C
C           Form the Fourier transform of the d array from the
C           quarter-wave sine transform
C
            DO 200 I = 1, M
               WORK(I,0) = 2.0D0*X(I,1)
  200       CONTINUE
            DO 240 K = 1, (N-1)/2
               DO 220 I = 1, M
                  WORK(I,K) = X(I,2*K+1) - X(I,2*K)
                  WORK(I,N-K) = -X(I,2*K+1) - X(I,2*K)
  220          CONTINUE
  240       CONTINUE
            IF (MOD(N,2).EQ.0) THEN
               N2 = N/2
               DO 300 I = 1, M
                  WORK(I,N2) = -2.0D0*X(I,N)
  300          CONTINUE
            END IF
C
C           Calculate the D array by taking the inverse transform of
C           the Re + i*Im arrays
C
            JFAIL = 0
            CALL C06GQF(M,N,WORK,JFAIL)
            CALL C06FQX(WORK,X,M,N,Q,NQ,TRIG)
C
C           Now calculate the transform from the D array
C
            DO 340 J = 1, N - 1
               DO 320 I = 1, M
                  X(I,J) = 0.5D0*(COS(DBLE(J)*PIBY2N)*(WORK(I,J)
     *                     -WORK(I,N-J))+SIN(DBLE(J)*PIBY2N)*(WORK(I,J)
     *                     +WORK(I,N-J)))
  320          CONTINUE
  340       CONTINUE
            DO 360 I = 1, M
               X(I,N) = WORK(I,0)
  360       CONTINUE
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
      ELSE IF (IERROR.EQ.6) THEN
         WRITE (REC(1),FMT=99994) DIRECT
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
99994 FORMAT (' ** ',A1,' is an invalid value of DIRECT')
      END
