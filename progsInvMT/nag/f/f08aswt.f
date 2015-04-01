      SUBROUTINE F08ASW(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLARF(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
C
C  Purpose
C  =======
C
C  ZLARF applies a complex elementary reflector H to a complex m by n
C  matrix C, from either the left or the right. H is represented in the
C  form
C
C        H = I - tau * v * v'
C
C  where tau is a complex scalar and v is a complex vector.
C
C  If tau = 0, then H is taken to be the unit matrix.
C
C  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
C  tau.
C
C  Arguments
C  =========
C
C  SIDE    (input) CHARACTER*1
C          = 'L': form  H * C
C          = 'R': form  C * H
C
C  M       (input) INTEGER
C          The number of rows of the matrix C.
C
C  N       (input) INTEGER
C          The number of columns of the matrix C.
C
C  V       (input) COMPLEX*16 array, dimension
C                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
C                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
C          The vector v in the representation of H. V is not used if
C          TAU = 0.
C
C  INCV    (input) INTEGER
C          The increment between elements of v. INCV <> 0.
C
C  TAU     (input) COMPLEX*16
C          The value tau in the representation of H.
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
C          On entry, the m-by-n matrix C.
C          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
C          or C * H if SIDE = 'R'.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDA >= max(1,M).
C
C  WORK    (workspace) COMPLEX*16 array, dimension
C                         (N) if SIDE = 'L'
C                      or (M) if SIDE = 'R'
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      COMPLEX*16        TAU
      INTEGER           INCV, LDC, M, N
      CHARACTER         SIDE
C     .. Array Arguments ..
      COMPLEX*16        C(LDC,*), V(*), WORK(*)
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERC
C     .. Executable Statements ..
C
      IF ((SIDE.EQ.'L' .OR. SIDE.EQ.'l')) THEN
C
C        Form  H * C
C
         IF (TAU.NE.ZERO) THEN
C
C           w := C' * v
C
            CALL ZGEMV('Conjugate transpose',M,N,ONE,C,LDC,V,INCV,ZERO,
     *                 WORK,1)
C
C           C := C - v * w'
C
            CALL ZGERC(M,N,-TAU,V,INCV,WORK,1,C,LDC)
         END IF
      ELSE
C
C        Form  C * H
C
         IF (TAU.NE.ZERO) THEN
C
C           w := C * v
C
            CALL ZGEMV('No transpose',M,N,ONE,C,LDC,V,INCV,ZERO,WORK,1)
C
C           C := C - w * v'
C
            CALL ZGERC(M,N,-TAU,WORK,1,V,INCV,C,LDC)
         END IF
      END IF
      RETURN
C
C     End of F08ASW (ZLARF)
C
      END
