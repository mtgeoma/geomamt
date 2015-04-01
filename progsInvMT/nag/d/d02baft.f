      SUBROUTINE D02BAF(X,XEND,N,Y,TOL,FCN,W,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     INTEGRATES THE DIFFERENTIAL EQUATIONS DEFINED BY
C     FCN AND WITH INITIAL VALUES IN Y FROM X TO XEND.
C     THE LOCAL ERROR IS CONTROLLED IN A MIXED ERROR TEST
C     USING TOL.
C     FCN
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02BAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL, X, XEND
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(N,7), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      INTEGER           I, IND
C     .. Local Arrays ..
      DOUBLE PRECISION  C(6), COMM(5), CON(3), COUT(14)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02PAF
C     .. Executable Statements ..
      IF (TOL.GT.0.D0 .AND. N.GT.0) GO TO 20
C     INPUT ERROR
      IND = 1
      GO TO 40
   20 C(1) = 0.D0
      I = 1
      CALL D02PAF(X,XEND,N,Y,C,TOL,FCN,COMM,CON,COUT,W,N,7,I)
      IND = 4
      IF (I.NE.0 .AND. I.NE.2 .AND. I.NE.3 .AND. I.NE.4) GO TO 40
      IND = I
      IF (IND.GT.2) IND = IND - 1
      IF (IND.NE.0 .AND. IND.NE.2) GO TO 40
      IF (3.D0*COUT(3).GE.COUT(8)) TOL = -TOL
      IF (IND.EQ.2) GO TO 40
      IF (X.NE.XEND .OR. C(1).NE.2.D0) IND = 4
   40 IFAIL = P01ABF(IFAIL,IND,SRNAME,0,P01REC)
      RETURN
      END
