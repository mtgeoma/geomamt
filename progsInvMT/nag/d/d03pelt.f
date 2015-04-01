      SUBROUTINE D03PEL(TIME,NIP,NPDE,X,U,FMON)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C     Dummy monitor function routine for D03PRF.
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TIME
      INTEGER           NIP, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  FMON(NIP), U(NPDE,*), X(NIP)
C     .. Executable Statements ..
      RETURN
      END
