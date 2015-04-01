      SUBROUTINE F08HEW(F,G,CS,SN,R)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLARTG(F,G,CS,SN,R)
C
C  Purpose
C  =======
C
C  DLARTG generate a plane rotation so that
C
C     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
C     [ -SN  CS  ]     [ G ]     [ 0 ]
C
C  This is a faster version of the BLAS1 routine DROTG, except for
C  the following differences:
C     F and G are unchanged on return.
C     If G=0, then CS=1 and SN=0.
C     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
C        floating point operations (saves work in DBDSQR when
C        there are zeros on the diagonal).
C
C  Arguments
C  =========
C
C  F       (input) DOUBLE PRECISION
C          The first component of vector to be rotated.
C
C  G       (input) DOUBLE PRECISION
C          The second component of vector to be rotated.
C
C  CS      (output) DOUBLE PRECISION
C          The cosine of the rotation.
C
C  SN      (output) DOUBLE PRECISION
C          The sine of the rotation.
C
C  R       (output) DOUBLE PRECISION
C          The nonzero component of the rotated vector.
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
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CS, F, G, R, SN
C     .. Local Scalars ..
      DOUBLE PRECISION  T, TT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
C
      IF (G.EQ.ZERO) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF (F.EQ.ZERO) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         IF (ABS(F).GT.ABS(G)) THEN
            T = G/F
            TT = SQRT(ONE+T*T)
            CS = ONE/TT
            SN = T*CS
            R = F*TT
         ELSE
            T = F/G
            TT = SQRT(ONE+T*T)
            SN = ONE/TT
            CS = T*SN
            R = G*TT
         END IF
      END IF
      RETURN
C
C     End of F08HEW (DLARTG)
C
      END
