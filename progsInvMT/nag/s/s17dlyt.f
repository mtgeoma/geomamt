      SUBROUTINE S17DLY(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-782 (DEC 1989).
C
C     Original name: CBUNK
C
C     S17DLY COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL.
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
C     IN S18DCZ AND THE EXPANSION FOR H(2,FNU,Z) IN S18DCY
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, TOL
      INTEGER           KODE, MR, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  AX, AY, XX, YY
C     .. External Subroutines ..
      EXTERNAL          S18DCY, S18DCZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DBLE
C     .. Executable Statements ..
C
      NZ = 0
      XX = DBLE(Z)
      YY = DIMAG(Z)
      AX = ABS(XX)*1.7321D0
      AY = ABS(YY)
      IF (AY.GT.AX) THEN
C        ---------------------------------------------------------------
C        ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C        APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C        AND HPI=PI/2
C        ---------------------------------------------------------------
         CALL S18DCY(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
      ELSE
C        ---------------------------------------------------------------
C        ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
C        -PI/3.LE.ARG(Z).LE.PI/3
C        ---------------------------------------------------------------
         CALL S18DCZ(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
      END IF
      RETURN
      END
