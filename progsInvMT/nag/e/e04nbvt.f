      SUBROUTINE E04NBV(N,NU,NRANK,LDR,LENV,LENW,R,U,V,W,C,S)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-589 (MAR 1988).
C     MARK 16 REVISED. IER-1057 (JUL 1993).
C
C     ==================================================================
C     E04NBV  modifies the  nrank*n  upper-triangular matrix  R  so that
C     Q*(R + v*w')  is upper triangular,  where  Q  is orthogonal,
C     v  and  w  are vectors, and the modified  R  overwrites the old.
C     Q  is the product of two sweeps of plane rotations (not stored).
C     If required,  the rotations are applied to the NU columns of
C     the matrix  U.
C
C     The matrix v*w' is an (LENV) by (LENW) matrix.
C     The vector v is overwritten.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version   October  1984.
C     Level-2 matrix routines added 22-Apr-1988.
C     This version of  E04NBV  dated 22-Apr-1988.
C     ==================================================================
C     .. Scalar Arguments ..
      INTEGER           LDR, LENV, LENW, N, NRANK, NU
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), R(LDR,*), S(N), U(N,*), V(N), W(N)
C     .. Local Scalars ..
      INTEGER           J
C     .. External Subroutines ..
      EXTERNAL          DAXPY, F06FQF, F06QSF, F06QWF, F06QXF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
      J = MIN(LENV,NRANK)
      IF (NRANK.GT.0) THEN
C        ---------------------------------------------------------------
C        Reduce  v to beta*e( j )  using a backward sweep of rotations
C        in planes (j-1, j), (j-2, j), ..., (1, j).
C        ---------------------------------------------------------------
         CALL F06FQF('Fixed','Backwards',J-1,V(J),V,1,C,S)
C
C        ---------------------------------------------------------------
C        Apply the sequence of rotations to U.
C        ---------------------------------------------------------------
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Backwards',J,NU,1,J,
     *                            C,S,U,N)
C
C        ---------------------------------------------------------------
C        Apply the sequence of rotations to R. This generates a spike in
C        the j-th row of R, which is stored in s.
C        ---------------------------------------------------------------
         CALL F06QWF('Left',N,1,J,C,S,R,LDR)
C
C        ---------------------------------------------------------------
C        Form  beta*e(j)*w' + R.  This a spiked matrix, with a row
C        spike in row j.
C        ---------------------------------------------------------------
         CALL DAXPY(MIN(J-1,LENW),V(J),W,1,S,1)
         CALL DAXPY(LENW-J+1,V(J),W(J),1,R(J,J),LDR)
C
C        ---------------------------------------------------------------
C        Eliminate the spike using a forward sweep of rotations in
C        planes (1, j), (2, j), ..., (j-1, j).
C        ---------------------------------------------------------------
         CALL F06QSF('Left',N,1,J,C,S,R,LDR)
C
C        ---------------------------------------------------------------
C        Apply the rotations to U.
C        ---------------------------------------------------------------
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Forwards',J,NU,1,J,C,
     *                            S,U,N)
      END IF
C
C     End of  E04NBV. (CMR1MD)
C
      END
