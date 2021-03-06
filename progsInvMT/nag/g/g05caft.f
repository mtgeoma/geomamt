      DOUBLE PRECISION FUNCTION G05CAF(X)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     Returns a pseudo-random number uniformly distributed between
C     A and B.
C
C     Pseudo-random numbers are generated by the auxiliary routine
C     G05CAY, 63 at a time, and stored in the array RV in common block
C     CG05CA. G05CAF copies one number from the array RV into X,
C     calling G05CAY to replenish RV when necessary.
C
C     This revised version of G05CAF has been introduced for
C     compatibility with the new routines G05FAF, G05FBF and G05FDF,
C     introduced at Mark 14.
C
C     Jeremy Du Croz, NAG Ltd, June 1989.
C
C     .. Parameters ..
      INTEGER                          LV
      PARAMETER                        (LV=63)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Scalars in Common ..
      INTEGER                          KV
C     .. Arrays in Common ..
      DOUBLE PRECISION                 RV(LV)
C     .. Local Scalars ..
      LOGICAL                          INIT
C     .. External Subroutines ..
      EXTERNAL                         G05CAY, G05CAZ
C     .. Common blocks ..
      COMMON                           /CG05CA/RV, KV
C     .. Save statement ..
      SAVE                             INIT, /CG05CA/
C     .. Data statements ..
      DATA                             INIT/.TRUE./
C     .. Executable Statements ..
C
C     Ensure that KV in common block /CG05CA/ has been initialized
C
      IF (INIT) CALL G05CAZ(INIT)
C
C     Replenish the buffer if necessary
C
      IF (KV.GE.LV) CALL G05CAY(.FALSE.)
C
      KV = KV + 1
      G05CAF = RV(KV)
      RETURN
      END
