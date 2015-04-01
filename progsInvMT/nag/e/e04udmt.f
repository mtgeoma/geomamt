      SUBROUTINE E04UDM(MODE,NCNLN,N,NROWJ,NEEDC,X,C,CJAC,NSTATE,IUSER,
     *                  USER)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     E04UDM  is a dummy  CONFUN  routine for use with  E04UCF  when
C     there are no constraints present.
C
C     .. Scalar Arguments ..
      INTEGER           MODE, N, NCNLN, NROWJ, NSTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  C(*), CJAC(NROWJ,*), USER(*), X(N)
      INTEGER           IUSER(*), NEEDC(*)
C     .. Executable Statements ..
      RETURN
C
C     End of  E04UDM.
C
      END
