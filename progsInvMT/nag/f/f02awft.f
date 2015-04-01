      SUBROUTINE F02AWF(AR,IAR,AI,IAI,N,WR,WK1,WK2,WK3,IFAIL)
C     MARK 3 RELEASE. NAG COPYRIGHT 1972.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-738 (DEC 1989).
C
C     1ST. APRIL 1972.
C     EIGENVALUES OF A COMPLEX HERMITIAN MATRX
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AWF')
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), WK1(N), WK2(N), WK3(N),
     *                  WR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XXXX
      INTEGER           ISAVE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01BCF, F02AVF
C     .. Executable Statements ..
      ISAVE = IFAIL
      CALL F01BCF(N,XXXX,AR,IAR,AI,IAI,WR,WK1,WK2,WK3)
      IFAIL = 1
      CALL F02AVF(N,X02AJF(),WR,WK1,IFAIL)
      IF (IFAIL.NE.0) IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
