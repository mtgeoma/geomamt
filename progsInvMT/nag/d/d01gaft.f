      SUBROUTINE D01GAF(X,Y,N,ANS,ER,IFAIL)
C
C     THIS SUBROUTINE INTEGRATES A FUNCTION (Y) SPECIFIED
C     NUMERICALLY AT N POINTS (X), WHERE N IS AT LEAST 4,
C     OVER THE RANGE X(1) TO X(N).  THE POINTS NEED NOT BE
C     EQUALLY SPACED, BUT SHOULD BE DISTINCT AND IN ASCENDING
C     OR DESCENDING ORDER.  AN ERROR ESTIMATE IS RETURNED.
C     THE METHOD IS DUE TO GILL AND MILLER.
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 7 REVISED IER-154 (DEC 1978)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01GAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANS, ER
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D1, D2, D3, H1, H2, H3, H4, R1, R2, R3, R4, S
      INTEGER           I, NN
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      ANS = 0.0D0
      ER = 0.0D0
      IF (N.GE.4) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
C
C     CHECK POINTS ARE STRICTLY INCREASING OR DECREASING
C
   20 H2 = X(2) - X(1)
      DO 80 I = 3, N
         H3 = X(I) - X(I-1)
         IF (H2*H3) 40, 60, 80
   40    IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
         RETURN
   60    IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
         RETURN
   80 CONTINUE
C
C     INTEGRATE OVER INITIAL INTERVAL
C
      D3 = (Y(2)-Y(1))/H2
      H3 = X(3) - X(2)
      D1 = (Y(3)-Y(2))/H3
      H1 = H2 + H3
      D2 = (D1-D3)/H1
      H4 = X(4) - X(3)
      R1 = (Y(4)-Y(3))/H4
      R2 = (R1-D1)/(H4+H3)
      H1 = H1 + H4
      R3 = (R2-D2)/H1
      ANS = H2*(Y(1)+H2*(D3/2.0D0-H2*(D2/6.0D0-(H2+2.0D0*H3)*R3/12.0D0))
     *      )
      S = -(H2**3)*(H2*(3.0D0*H2+5.0D0*H4)+10.0D0*H3*H1)/60.0D0
      R4 = 0.0D0
C
C     INTEGRATE OVER CENTRAL PORTION OF RANGE
C
      NN = N - 1
      DO 120 I = 3, NN
         ANS = ANS + H3*((Y(I)+Y(I-1))/2.0D0-H3*H3*(D2+R2+(H2-H4)*R3)
     *         /12.0D0)
         C = H3**3*(2.0D0*H3*H3+5.0D0*(H3*(H4+H2)+2.0D0*H4*H2))/120.0D0
         ER = ER + (C+S)*R4
         IF (I.NE.3) S = C
         IF (I.EQ.3) S = S + 2.0D0*C
         IF (I-N+1) 100, 140, 100
  100    H1 = H2
         H2 = H3
         H3 = H4
         D1 = R1
         D2 = R2
         D3 = R3
         H4 = X(I+2) - X(I+1)
         R1 = (Y(I+2)-Y(I+1))/H4
         R4 = H4 + H3
         R2 = (R1-D1)/R4
         R4 = R4 + H2
         R3 = (R2-D2)/R4
         R4 = R4 + H1
         R4 = (R3-D3)/R4
  120 CONTINUE
C
C     INTEGRATE OVER FINAL INTERVAL
C
  140 CONTINUE
      ANS = ANS + H4*(Y(N)-H4*(R1/2.0D0+H4*(R2/6.0D0+(2.0D0*H3+H4)
     *      *R3/12.0D0)))
      ER = ER - H4**3*R4*(H4*(3.0D0*H4+5.0D0*H2)+10.0D0*H3*(H2+H3+H4))
     *     /60.0D0 + S*R4
      ANS = ANS + ER
      IFAIL = 0
      RETURN
      END
