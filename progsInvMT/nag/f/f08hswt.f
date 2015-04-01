      SUBROUTINE F08HSW(F,G,CS,SN,R)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1650 (JUN 1995).
C     ENTRY             ZLARTG(F,G,CS,SN,R)
C
C  Purpose
C  =======
C
C  ZLARTG generates a plane rotation so that
C
C     [  CS  SN  ]     [ F ]     [ R ]
C     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
C     [ -SN  CS  ]     [ G ]     [ 0 ]
C
C  This is a faster version of the BLAS1 routine ZROTG, except for
C  the following differences:
C     F and G are unchanged on return.
C     If G=0, then CS=1 and SN=0.
C     If F=0, then CS=0 and SN is chosen so that R is real.
C
C  Arguments
C  =========
C
C  F       (input) COMPLEX*16
C          The first component of vector to be rotated.
C
C  G       (input) COMPLEX*16
C          The second component of vector to be rotated.
C
C  CS      (output) DOUBLE PRECISION
C          The cosine of the rotation.
C
C  SN      (output) COMPLEX*16
C          The sine of the rotation.
C
C  R       (output) COMPLEX*16
C          The nonzero component of the rotated vector.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=0.0D0)
C     .. Scalar Arguments ..
      COMPLEX*16        F, G, R, SN
      DOUBLE PRECISION  CS
C     .. Local Scalars ..
      COMPLEX*16        FS, GS, T
      DOUBLE PRECISION  D, F1, F2, G1, G2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCONJG, DIMAG, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  ABS1, ABSSQ
C     .. Statement Function definitions ..
      ABS1(T) = ABS(DBLE(T)) + ABS(DIMAG(T))
      ABSSQ(T) = DBLE(T)**2 + DIMAG(T)**2
C     .. Executable Statements ..
C
      IF (G.EQ.CZERO) THEN
         CS = ONE
         SN = CZERO
         R = F
      ELSE IF (F.EQ.CZERO) THEN
         CS = ZERO
         SN = DCONJG(G)/ABS(G)
         R = ABS(G)
      ELSE
         F1 = ABS1(F)
         G1 = ABS1(G)
         IF (F1.GE.G1) THEN
            GS = G/F1
            G2 = ABSSQ(GS)
            FS = F/F1
            F2 = ABSSQ(FS)
            D = SQRT(ONE+G2/F2)
            CS = ONE/D
            R = F*D
            SN = DCONJG(GS)*FS*(CS/F2)
         ELSE
            FS = F/G1
            F2 = ABSSQ(FS)
            GS = G/G1
            G2 = ABSSQ(GS)
            D = G1*SQRT(F2+G2)
            F1 = ABS(F)
            FS = F/F1
            CS = F1/D
            R = FS*D
            SN = (DCONJG(G)/D)*FS
         END IF
      END IF
      RETURN
C
C     End of F08HSW (ZLARTG)
C
      END
