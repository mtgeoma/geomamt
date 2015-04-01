      SUBROUTINE E01SFF(M,X,Y,F,RNW,FNODES,PX,PY,PF,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     This function serves as the evaluation module in a quadratic
C     Shepard interpolation process. This module should be called
C     after a call to subroutine E01SEF.
C
C     Routine created - December 1986
C     Author          - Richard Franke
C                       Naval Postgraduate School Monterey,
C                       California  93940
C                       Adapted for Nag by H.Scullion (Leic Univ.)
C                       and I. Gladwell (Nag Ltd.)
C
C     Input Parameters:
C
C            M  -  The number of data points.
C                  Unchanged on exit.
C
C        X,Y,F  -  The data points, (X(I),Y(I),F(I),I=1,M).
C                  Unchanged on exit.
C
C          RNW  -  The radius for the weight functions. As
C                  defined in E01SEF.
C                  Unchanged on exit.
C
C       FNODES  -  Real array of dimension at least (5*M).
C                  This array is used to store the coefficients for
C                  the nodal functions.
C                  Unchanged on exit.
C
C        PX,PY  -  The point (PX,PY) where the function is to be
C                  evaluated.
C                  Unchanged on exit.
C
C     Output Parameter:
C
C           PF  -  The value of the interpolant at point (PX,PY).
C
C     On exit, if IFAIL = 0, normal return.
C                    = 1, M is .lt. 3.
C                    = 2, evaluation point is outside support region.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E01SFF')
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PF, PX, PY, RNW
      INTEGER           IFAIL, M
C     .. Array Arguments ..
      DOUBLE PRECISION  F(M), FNODES(5,M), X(M), Y(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  BOT, D, TOP, W, XD, YD
      INTEGER           IER, K, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IER = 0
      NREC = 0
      PF = ZERO
      IF (M.LT.3) THEN
         IER = 1
         NREC = 1
         WRITE (REC,FMT=99999) M
      ELSE
         TOP = ZERO
         BOT = ZERO
         DO 20 K = 1, M
            XD = PX - X(K)
            YD = PY - Y(K)
            D = SQRT(XD**2+YD**2)
            IF (D.EQ.ZERO) THEN
C              Point (PX,PY) is equal to a data point.
               PF = F(K)
               GO TO 40
            ELSE IF (D.LT.RNW) THEN
C              We are inside the region of influence of data point K,
C              so add in its contribution to the evaluation.
               W = ((D-RNW)/(D*RNW))**2
               TOP = TOP + (F(K)+(FNODES(1,K)+FNODES(3,K)*XD+FNODES(4,K)
     *               *YD)*XD+(FNODES(2,K)+FNODES(5,K)*YD)*YD)*W
               BOT = BOT + W
            END IF
   20    CONTINUE
C
         IF (BOT.EQ.ZERO) THEN
C           The evaluation point is outside the region of influence of
C           every data point.
            IER = 2
            NREC = 2
            WRITE (REC,FMT=99998) PX, PY, RNW
         ELSE
            PF = TOP/BOT
         END IF
      END IF
C
   40 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, M .lt. 3: M =',I16,'.')
99998 FORMAT (1X,'** On entry, the evaluation point (',1P,D13.5,',',
     *       D13.5,') is outside the',/4X,'support region of the data ',
     *       'points defined by RNW =',D13.5,' .')
      END
