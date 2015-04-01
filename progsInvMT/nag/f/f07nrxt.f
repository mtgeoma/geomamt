      SUBROUTINE F07NRX(A,B,C,RT1,RT2,EVSCAL,CS1,SN1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1704 (JUN 1995).
C     ENTRY             ZLAESY(A,B,C,RT1,RT2,EVSCAL,CS1,SN1)
C
C  Purpose
C  =======
C
C  ZLAESY computes the eigendecomposition of a 2x2 symmetric matrix
C     ( ( A, B );( B, C ) )
C  provided the norm of the matrix of eigenvectors is larger than
C  some threshold value.
C
C  RT1 is the eigenvalue of larger absolute value, and RT2 of
C  smaller absolute value.  If the eigenvectors are computed, then
C  on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence
C
C  [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]
C  [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]
C
C  Arguments
C  =========
C
C  A       (input) COMPLEX
C          The ( 1, 1 ) entry of input matrix.
C
C  B       (input) COMPLEX
C          The ( 1, 2 ) entry of input matrix.  The ( 2, 1 ) entry is
C          also given by B, since the 2 x 2 matrix is symmetric.
C
C  C       (input) COMPLEX
C          The ( 2, 2 ) entry of input matrix.
C
C  RT1     (output) COMPLEX
C          The eigenvalue of larger modulus.
C
C  RT2     (output) COMPLEX
C          The eigenvalue of smaller modulus.
C
C  EVSCAL  (output) COMPLEX
C          The complex value by which the eigenvector matrix was scaled
C          to make it orthonormal.  If EVSCAL is zero, the eigenvectors
C          were not computed.  This means one of two things:  the 2 x 2
C          matrix could not be diagonalized, or the norm of the matrix
C          of eigenvectors before scaling was larger than the threshold
C          value THRESH (set below).
C
C  CS1     (output) COMPLEX
C  SN1     (output) COMPLEX
C          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector
C          for RT1.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=(1.0D0,0.0D0))
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
      DOUBLE PRECISION  THRESH
      PARAMETER         (THRESH=0.1D0)
C     .. Scalar Arguments ..
      COMPLEX*16        A, B, C, CS1, EVSCAL, RT1, RT2, SN1
C     .. Local Scalars ..
      COMPLEX*16        S, T, TMP
      DOUBLE PRECISION  BABS, EVNORM, TABS, Z
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C
C     Special case:  The matrix is actually diagonal.
C     To avoid divide by zero later, we treat this case separately.
C
      IF (ABS(B).EQ.ZERO) THEN
         RT1 = A
         RT2 = C
         IF (ABS(RT1).LT.ABS(RT2)) THEN
            TMP = RT1
            RT1 = RT2
            RT2 = TMP
            CS1 = ZERO
            SN1 = ONE
         ELSE
            CS1 = ONE
            SN1 = ZERO
         END IF
      ELSE
C
C        Compute the eigenvalues and eigenvectors.
C        The characteristic equation is
C           lamba **2 - (A+C) lamba + (A*C - B*B)
C        and we solve it using the quadratic formula.
C
         S = (A+C)*HALF
         T = (A-C)*HALF
C
C        Take the square root carefully to avoid over/under flow.
C
         BABS = ABS(B)
         TABS = ABS(T)
         Z = MAX(BABS,TABS)
         IF (Z.GT.ZERO) T = Z*SQRT((T/Z)**2+(B/Z)**2)
C
C        Compute the two eigenvalues.  RT1 and RT2 are exchanged
C        if necessary so that RT1 will have the greater magnitude.
C
         RT1 = S + T
         RT2 = S - T
         IF (ABS(RT1).LT.ABS(RT2)) THEN
            TMP = RT1
            RT1 = RT2
            RT2 = TMP
         END IF
C
C        Choose CS1 = 1 and SN1 to satisfy the first equation, then
C        scale the components of this eigenvector so that the matrix
C        of eigenvectors X satisfies  X * X' = I .  (No scaling is
C        done if the norm of the eigenvalue matrix is less than THRESH.)
C
         SN1 = (RT1-A)/B
         TABS = ABS(SN1)
         IF (TABS.GT.ONE) THEN
            T = TABS*SQRT((ONE/TABS)**2+(SN1/TABS)**2)
         ELSE
            T = SQRT(CONE+SN1*SN1)
         END IF
         EVNORM = ABS(T)
         IF (EVNORM.GE.THRESH) THEN
            EVSCAL = CONE/T
            CS1 = EVSCAL
            SN1 = SN1*EVSCAL
         ELSE
            EVSCAL = ZERO
         END IF
      END IF
      RETURN
C
C     End of F07NRX (ZLAESY)
C
      END
