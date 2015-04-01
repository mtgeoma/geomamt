      SUBROUTINE D03PCL(TIME,NIP,NPDE,X,U,R,FMON)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-979 (JUN 1993).
C--------------------------------------------------------------
C   Dummy  monitor  function routine for remeshing.
C   R(j,i) contains the flux for jth PDE. at point Y(i)
C   where Y(1) = X(1) , Y(NIP+1) = X(NIP)
C   Y(I) =(X(I+1) - X(I)) * 0.5   I = 2,..., NIP
C   The monitor function should be returned in the array FMON.
C   VP changed.
C--------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TIME
      INTEGER           NIP, NPDE
C     .. Array Arguments ..
      DOUBLE PRECISION  FMON(NIP), R(NPDE,*), U(NPDE,*), X(NIP)
C     .. Executable Statements ..
      RETURN
      END
