      SUBROUTINE E01SAF(M,X,Y,F,TRIANG,GRADS,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Creates a Thiessen triangulation of (M,X,Y) and computes
C     derivatives required for interpolation.
C
C     Input arguments:
C     M is the number of data points (nodes).
C     X, Y, F are the scattered data to be interpolated, F = F(X,Y).
C
C     Output arguments:
C     TRIANG contains the data structure defining the triangulation;
C      used as parameters IADJ and IEND to subroutine E01SAZ.
C     GRADS contains the estimated partial derivatives at the nodes;
C      first row contains derivatives with respect to X, second row
C      with respect to Y. Used as parameter ZXZY to subroutine E01SAZ.
C
C     Parameters:
C     TOL is a convergence criterion for estimating gradients at nodes;
C      TOL .ge. 0.0; TOL = 0.01 is sufficient.
C     MAXIT is the maximum number of iterations allowed for computation
C      of estimated gradients; MAXIT .ge. 0.
C
C     M, X, Y, F, TRIANG and GRADS should be used as input to NAG
C     library routine E01SBF to compute interpolated values of F(X,Y).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01SAF')
      DOUBLE PRECISION  TOL, ZERO
      PARAMETER         (TOL=0.01D0,ZERO=0.0D0)
      INTEGER           MAXIT
      PARAMETER         (MAXIT=6)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), GRADS(2,M), X(M), Y(M)
      INTEGER           TRIANG(7*M)
C     .. Local Scalars ..
      INTEGER           COINC1, COINC2, I, IER, NIT, NREC, WK1, WK2
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E01SAY, E01SAZ
C     .. Executable Statements ..
      IER = 1
      NREC = 1
      IF (M.LT.3) THEN
         WRITE (REC,FMT=99999) M
      ELSE
C        Split the triangulation information array.
         WK1 = 1
         WK2 = 6*M + 1
C        Set up the Thiessen triangulation.
         IER = 0
         NREC = 0
         CALL E01SAY(M,X,Y,TRIANG(WK1),TRIANG(WK2),IER)
         IF (IER.EQ.2) THEN
            NREC = 1
            WRITE (REC,FMT=99998)
         ELSE
C           Initialise gradients to zero for E01SAZ.
            DO 20 I = 1, M
               GRADS(1,I) = ZERO
               GRADS(2,I) = ZERO
   20       CONTINUE
C           Estimate the gradients at the nodes.
            NIT = MAXIT
            CALL E01SAZ(M,X,Y,F,TRIANG(WK1),TRIANG(WK2),TOL,NIT,GRADS,
     *                  COINC1,COINC2,IER)
            IF (IER.EQ.3) THEN
               NREC = 2
               WRITE (REC,FMT=99997) COINC1, COINC2, X(COINC1),
     *           Y(COINC1)
            END IF
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, M .lt. 3: M =',I16,'.')
99998 FORMAT (1X,'** On entry, all (X,Y) pairs are collinear.')
99997 FORMAT (1X,'** On entry, (X(I),Y(I)) = (X(J),Y(J)) for I,J =',2I8,
     *       /4X,'X(I), Y(I) = ',1P,2D13.5,' .')
      END
