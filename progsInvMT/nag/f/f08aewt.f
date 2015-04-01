      SUBROUTINE F08AEW(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLARF(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
C
C  Purpose
C  =======
C
C  DLARF applies a real elementary reflector H to a real m by n matrix
C  C, from either the left or the right. H is represented in the form
C
C        H = I - tau * v * v'
C
C  where tau is a real scalar and v is a real vector.
C
C  If tau = 0, then H is taken to be the unit matrix.
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
C  V       (input) DOUBLE PRECISION array, dimension
C                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
C                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
C          The vector v in the representation of H. V is not used if
C          TAU = 0.
C
C  INCV    (input) INTEGER
C          The increment between elements of v. INCV <> 0.
C
C  TAU     (input) DOUBLE PRECISION
C          The value tau in the representation of H.
C
C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C          On entry, the m by n matrix C.
C          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
C          or C * H if SIDE = 'R'.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M).
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension
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
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TAU
      INTEGER           INCV, LDC, M, N
      CHARACTER         SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,*), V(*), WORK(*)
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER
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
            CALL DGEMV('Transpose',M,N,ONE,C,LDC,V,INCV,ZERO,WORK,1)
C
C           C := C - v * w'
C
            CALL DGER(M,N,-TAU,V,INCV,WORK,1,C,LDC)
         END IF
      ELSE
C
C        Form  C * H
C
         IF (TAU.NE.ZERO) THEN
C
C           w := C * v
C
            CALL DGEMV('No transpose',M,N,ONE,C,LDC,V,INCV,ZERO,WORK,1)
C
C           C := C - w * v'
C
            CALL DGER(M,N,-TAU,WORK,1,V,INCV,C,LDC)
         END IF
      END IF
      RETURN
C
C     End of F08AEW (DLARF)
C
      END
