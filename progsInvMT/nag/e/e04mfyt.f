      SUBROUTINE E04MFY(NRZ,LDR,R,RZZ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ******************************************************************
C     E04MFY  loads the last column of the  NRZ x NRZ  triangular factor
C     Rz  with the multiple  RZZ  of the  NRZ-th unit vector.
C
C     Original version written by PEG,  23-Jul-87.
C     This version of  E04MFY  dated 17-Jul-90.
C     ******************************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RZZ
      INTEGER           LDR, NRZ
C     .. Array Arguments ..
      DOUBLE PRECISION  R(LDR,*)
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Executable Statements ..
C
      IF (NRZ.EQ.0) RETURN
C
      CALL F06FBF(NRZ-1,ZERO,R(1,NRZ),1)
      R(NRZ,NRZ) = RZZ
C
      RETURN
C
C     End of  E04MFY.  (LPCOLR)
C
      END
