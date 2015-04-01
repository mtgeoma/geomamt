      SUBROUTINE E02BDF(NCAP7,K,C,DEFINT,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************
C     *                                                *
C     *         NAG LIBRARY SUBROUTINE  E02BDF         *
C     *                                                *
C     *  DEFINITE INTEGRAL OF CUBIC SPLINE FROM ITS    *
C     *  B-SPLINE REPRESENTATION                       *
C     *                                                *
C     *  ROUTINE CREATED ... 17 NOV 1977               *
C     *  LATEST UPDATE ....  24 APR 1978               *
C     *  RELEASE NUMBER ...  01                        *
C     *  AUTHORS ... MAURICE G. COX AND                *
C     *              J. GEOFFREY HAYES, N.P.L.         *
C     *                                                *
C     **************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02BDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DEFINT
      INTEGER           IFAIL, NCAP7
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NCAP7), K(NCAP7)
C     .. Local Scalars ..
      DOUBLE PRECISION  FOUR, SUM, ZERO
      INTEGER           IERROR, J, NCAP, NCAP3, NCAP4
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      ZERO = 0.0D+0
      FOUR = 4.0D+0
C
C     *********  DATA VALIDATION  *********
C
C     CHECK WHETHER AT LEAST ONE INTERVAL HAS BEEN SPECIFIED -
C
      IERROR = 1
      IF (NCAP7.LT.8) GO TO 60
C
C     CHECK WHETHER THE KNOTS SUPPLIED ARE REASONABLE -
C
      IERROR = 2
      NCAP = NCAP7 - 7
      NCAP4 = NCAP + 4
      IF (K(4).GE.K(NCAP4)) GO TO 60
      IF (K(1).NE.K(2) .OR. K(2).NE.K(3) .OR. K(3).NE.K(4)) GO TO 60
      IF (K(NCAP4).NE.K(NCAP+5) .OR. K(NCAP+5).NE.K(NCAP+6)
     *    .OR. K(NCAP+6).NE.K(NCAP+7)) GO TO 60
      DO 20 J = 5, NCAP4
         IF (K(J).LT.K(J-1)) GO TO 60
   20 CONTINUE
C
C     *********  COMPUTATION  *********
C
      IERROR = 0
      NCAP3 = NCAP + 3
      SUM = ZERO
      DO 40 J = 1, NCAP3
         SUM = SUM + (K(J+4)-K(J))*C(J)
   40 CONTINUE
      DEFINT = SUM/FOUR
C
C     *********  ERROR DIAGNOSTICS  *********
C
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
C
C     END OF SUBROUTINE  E02BDF
C
      END
