      SUBROUTINE G05HDV(A,N,C,LDC,TOL,R,NR,IERR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1683 (JUN 1995).
C
C     G05HDV CONVERTS A VECTOR OF MEANS AND A COVARIANCE
C     MATRIX INTO A FORM THAT IS EFFICIENT FOR THE GENERATION OF
C     MULTIVARIATE NORMAL VECTORS WITH THE SPECIFIED PARAMETERS.
C     THE REFERENCE VECTOR HAS THE FOLLOWING VALUES AT THE
C     FOLLOWING LOCATIONS -
C     1) THE DIMENSION OF THE DISTRIBUTION (N)
C     2) THE RANK OF THE COVARIANCE MATRIX
C     3) THE MEAN OF THE DISTRIBUTION
C     3+N) THE MATRIX OF EIGENVECTORS, STORED BY COLUMNS, (FROM THE
C          SPECTRAL DECOMPOSITION) OF THE COVARIANCE MATRIX. ON EXIT
C          THEY ARE SCALED BY THE SQUARE ROOT OF THE EIGENVECTORS
C          AND THE SIGNS ARE CHANGED TO ENSURE THAT THE FIRST NON-ZERO
C          ELEMENT OF EACH EIGENVECTOR IS POSITIVE
C     3+N+N*N) THE EIGENVALUES IN ASCENDING ORDER OF THE COVARIANCE
C              MATRIX
C     3+N+N*N+N) WORKSPACE, STORES EIGENVALUES IN DESCENDING ODER
C
C     Adapted from G05EAF which was written by N.M.MACLAREN
C        UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY.
C
C     .. Parameters ..
      DOUBLE PRECISION  HALF, ZERO
      PARAMETER         (HALF=0.5D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IERR, LDC, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), C(LDC,N), R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, SQEPS, T1, T2
      INTEGER           I, IFAIL, IRANK, J1, J2, J3, L, N2, NIR, NN, NN2
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           IDAMAX, F06KLF
      EXTERNAL          X02AJF, IDAMAX, F06KLF
C     .. External Subroutines ..
      EXTERNAL          F02ABF, DSCAL, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Executable Statements ..
C
C     Check the arguments for sensible values and consistency.
C     These arguments do not need checking when called from G05HDF.
C
C      IF (N.LT.1) THEN
C         IERR = 1
C      ELSE IF (NR.LT.(N+1)*(N+2)) THEN
C         IERR = 2
C      ELSE IF (LDC.LT.N) THEN
C         IERR = 3
C      ELSE IF (TOL.LT.ZERO .OR. TOL.GT.1.0D0/DBLE(10*N)) THEN
C         IERR = 4
C      ELSE
      IERR = 0
C
C        copy the means
C
      R(1) = DBLE(N) + HALF
      CALL DCOPY(N,A,1,R(3),1)
C
C        Now perform the spectral decomposition of the covariance matrix
C
      N2 = 2 + N
      NN = N*N
      NN2 = N2 + NN + N
      IFAIL = 1
      CALL F02ABF(C,LDC,N,R(N2+NN+1),R(N2+1),N,R(NN2+1),IFAIL)
      IF (IFAIL.EQ.1) THEN
         IERR = 5
         RETURN
      END IF
C
C        Find the rank of the matrix by examining the eigenvalues.
C
      DO 20 I = 1, N
         R(NN2+I) = R(NN2-I+1)
   20 CONTINUE
      IRANK = F06KLF(N,R(NN2+1),1,TOL)
      R(2) = DBLE(IRANK) + HALF
      NIR = N - IRANK
C
C     Scale the eigenvectors by the square root of the corresponding
C     eigenvector, changing the sign so that the first non-zero element
C     is always positive.
C
      SQEPS = SQRT(X02AJF())
      DO 40 I = 1, IRANK
         L = N2 + (NIR+I-1)*N
         ALPHA = SQRT(R(N2+NN+NIR+I))
         J1 = IDAMAX(N,R(L+1),1)
         IF (J1.EQ.1) THEN
            IF (R(L+1).LT.ZERO) ALPHA = -ALPHA
         ELSE IF (J1.EQ.N) THEN
            J2 = IDAMAX(N-1,R(L+1),1)
            T1 = R(L+N)
            T2 = R(L+J2)
            IF ((ABS(T1)-ABS(T2)).LT.SQEPS) THEN
               IF (T2.LT.ZERO) ALPHA = -ALPHA
            ELSE
               IF (T1.LT.ZERO) ALPHA = -ALPHA
            END IF
         ELSE
            J2 = IDAMAX(J1-1,R(L+1),1)
            J3 = IDAMAX(N-J1,R(L+J1+1),1) + J1
            IF (ABS(R(L+J2)).LT.ABS(R(L+J3))) J2 = J3
            T1 = R(L+J1)
            T2 = R(L+J2)
            IF ((ABS(T1)-ABS(T2)).LT.SQEPS) THEN
               IF (J1.LT.J2) THEN
                  IF (T1.LT.ZERO) ALPHA = -ALPHA
               ELSE
                  IF (T2.LT.ZERO) ALPHA = -ALPHA
               END IF
            ELSE
               IF (T1.LT.ZERO) ALPHA = -ALPHA
            END IF
         END IF
         CALL DSCAL(N,ALPHA,R(N2+(NIR+I-1)*N+1),1)
   40 CONTINUE
C      END IF
      RETURN
      END
