      SUBROUTINE D03PDG(DINPDF,DINPJF,NPDE,NPTS,X,U,NV,V)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C      ROUTINE FOR P.D.E. INITIAL VALUES.
C     .. Scalar Arguments ..
      INTEGER           NPDE, NPTS, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPDE,NPTS), V(*), X(NPTS)
C     .. Subroutine Arguments ..
      EXTERNAL          DINPDF, DINPJF
C     .. Executable Statements ..
C
      CALL DINPDF(NPDE,NPTS,X,U)
      CALL DINPJF(NPDE,NPTS,X,U,NV,V)
C
      RETURN
C
      END
