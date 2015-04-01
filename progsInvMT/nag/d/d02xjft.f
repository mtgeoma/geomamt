      SUBROUTINE D02XJF(XSOL,SOL,M,W,NEQMAX,IW,NEQ,X,NQ,HU,H,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-545 (FEB 1987).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02XJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HU, X, XSOL
      INTEGER           IFAIL, IW, M, NEQ, NEQMAX, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  SOL(M), W(NEQMAX,IW)
C     .. Local Scalars ..
      INTEGER           ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02XJZ
C     .. Executable Statements ..
      ISAVE = 1
      IF (M.LT.1 .OR. NEQMAX.LT.1 .OR. NEQ.LT.1 .OR. M.GT.NEQ .OR.
     *    NEQ.GT.NEQMAX .OR. NQ.LT.1 .OR. IW.LT.NQ+1) GO TO 20
      IF (H.EQ.0.0D0 .OR. HU.EQ.0.0D0) THEN
         ISAVE = 2
         GO TO 20
      END IF
      CALL D02XJZ(XSOL,0,W,NEQMAX,SOL,ISAVE,M,H,X,HU,NQ)
      IF ((XSOL-X)*H.GT.0.0D0) ISAVE = 3
   20 IFAIL = P01ABF(IFAIL,ISAVE,SRNAME,0,P01REC)
      RETURN
      END
