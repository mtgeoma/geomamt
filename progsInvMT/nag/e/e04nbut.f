      SUBROUTINE E04NBU(N,NU,NRANK,LDR,I,J,R,U,C,S)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     This version dated 8-June-1988. (F06 routines included.)
C
C
C***********************************************************************
C     E04NBU  interchanges the  I-th  and  J-th  (I .LT. J)  columns of
C     an  NRANK*N  upper-trapezoidal matrix  R   and restores the
C     resulting matrix to upper-trapezoidal form using two sweeps of
C     plane rotations applied on the left.  R is overwritten.
C
C     If NU .GT. 0,  the rotations are applied to the  nu  columns of
C     the matrix  U.
C
C     Systems Optimization Laboratory, Stanford University.
C     Original version written 31-October-1984.
C     Level-2 matrix routines added 13-May-1988.
C     This version of  E04NBU  dated  30-May-1988.
C***********************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           I, J, LDR, N, NRANK, NU
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), R(LDR,*), S(N), U(N,*)
C     .. Local Scalars ..
      INTEGER           LENJ
C     .. External Subroutines ..
      EXTERNAL          F06FBF, F06FQF, F06QSF, F06QWF, F06QXF, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     Swap the elements of the i-th and j-th columns of R on, or above,
C     the main diagonal.
C
      CALL DSWAP(MIN(I,NRANK),R(1,I),1,R(1,J),1)
      LENJ = MIN(J,NRANK)
C
      IF (LENJ.GT.I) THEN
C        ---------------------------------------------------------------
C        Reduce elements  r(i+1,j), ..., r(lenj,j)  to  beta*e(lenj)
C        using a backward sweep in planes
C        (lenj-1,lenj), (lenj-2,lenj), ..., (i+1,lenj).
C        If required, apply the sequence of rotations to U.
C        ---------------------------------------------------------------
         CALL F06FQF('Fixed','Backwards',LENJ-I-1,R(LENJ,J),R(I+1,J),1,
     *               C(I+1),S(I+1))
C
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Backwards',N,NU,I+1,
     *                            LENJ,C,S,U,N)
C
C        Put zeros into the j-th column of R in positions corresponding
C        to the sub-diagonals of the i-th column.
C
         S(I) = R(LENJ,J)
         CALL F06FBF(LENJ-I,ZERO,R(I+1,J),1)
C
C        Apply the sequence of rotations to R.  This generates a spike
C        in the lenj-th row of R, which is stored in S.
C
         CALL F06QWF('Left',N,I+1,LENJ,C,S,R,LDR)
C
C        Eliminate the spike using a forward sweep in planes
C        (i,lenj), (i+1,lenj), ..., (lenj-1,lenj).
C        If necessary, apply the sequence of rotations to U.
C
         CALL F06QSF('Left',N,I,LENJ,C,S,R,LDR)
C
         IF (NU.GT.0) CALL F06QXF('Left','Bottom','Forwards',LENJ,NU,I,
     *                            LENJ,C,S,U,N)
      END IF
C
      RETURN
C
C     End of  E04NBU
C
      END
