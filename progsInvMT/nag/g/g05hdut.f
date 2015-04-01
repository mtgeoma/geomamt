      SUBROUTINE G05HDU(Z,N,R,NR,IERR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     G05HDU returns a multivariate Normal vector from the
C     parameters converted by routine G05HDV.
C     If we want y - Normal(mu,Sigma) where mu is a mean vector
C     and Sigma the covariance matrix, then by a spectral decomposition
C       Sigma = V'DV  where the columns of V' store the eigenvectors
C                     and D is a diagonal matrix storing the eigenvalues
C       Then z = V'sqrt(D)z where z is the vector of standard normals.
C
C     Adapted from G05EZF by M. Eagle, NAG Ltd., December 1990
C     G05EZF written by N.M.MACLAREN
C                       UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IERR, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR), Z(N)
C     .. Local Scalars ..
      INTEGER           I, IRANK, N2, NIR, NN
C     .. External Functions ..
      DOUBLE PRECISION  G05DDF
      EXTERNAL          G05DDF
C     .. External Subroutines ..
      EXTERNAL          F06EFF, F06PAF
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Executable Statements ..
C
C     Check the arguments for sensible values and consistency.
C
      IF ((N.LT.1) .OR. (INT(R(1)).NE.N)) THEN
         IERR = 1
      ELSE IF (NR.LT.(N+1)*(N+2)) THEN
         IERR = 2
      ELSE
         IERR = 0
         N2 = 2 + N
         NN = N*N
C        Recover the rank
         IRANK = INT(R(2))
         NIR = N - IRANK
C
C        Generate N standard normal numbers.
C
         DO 20 I = 1, IRANK
            R(N2+NN+NIR+I) = G05DDF(ZERO,ONE)
   20    CONTINUE
C
C        Initialise the result vector to the mean
C
         CALL F06EFF(N,R(3),1,Z,1)
C
C        Multiply the vector of N standard Normals by the scaled matrix
C        of eigenvectors and add to the mean.
C
         CALL F06PAF('N',N,IRANK,ONE,R(N2+NIR*N+1),N,R(N2+NN+NIR+1),1,
     *               ONE,Z,1)
      END IF
      RETURN
      END
