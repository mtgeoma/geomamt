      SUBROUTINE F07NRW(N,CX,INCX,CY,INCY,C,S)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLACRT(N,CX,INCX,CY,INCY,C,S)
C
C  Purpose
C  =======
C
C  ZLACRT applies a plane rotation, where the cos and sin (C and S) are
C  complex and the vectors CX and CY are complex.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of elements in the vectors CX and CY.
C
C  CX      (input/output) COMPLEX array, dimension (N)
C          On input, the vector X.
C          On output, CX is overwritten with C*X + S*Y.
C
C  INCX    (input) INTEGER
C          The increment between successive values of CY.  INCX <> 0.
C
C  CY      (input/output) COMPLEX array, dimension (N)
C          On input, the vector Y.
C          On output, CY is overwritten with -S*X + C*Y.
C
C  INCY    (input) INTEGER
C          The increment between successive values of CY.  INCX <> 0.
C
C  C       (input) COMPLEX
C  S       (input) COMPLEX
C          C and S define a complex rotation
C             [  C   S  ]
C             [ -S   C  ]
C          where C*C + S*S = 1.0.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Scalar Arguments ..
      COMPLEX*16        C, S
      INTEGER           INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16        CX(*), CY(*)
C     .. Local Scalars ..
      COMPLEX*16        CTEMP
      INTEGER           I, IX, IY
C     .. Executable Statements ..
C
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 40
C
C     Code for unequal increments or equal increments not equal to 1
C
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 20 I = 1, N
         CTEMP = C*CX(IX) + S*CY(IY)
         CY(IY) = C*CY(IY) - S*CX(IX)
         CX(IX) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   20 CONTINUE
      RETURN
C
C     Code for both increments equal to 1
C
   40 CONTINUE
      DO 60 I = 1, N
         CTEMP = C*CX(I) + S*CY(I)
         CY(I) = C*CY(I) - S*CX(I)
         CX(I) = CTEMP
   60 CONTINUE
      RETURN
      END
