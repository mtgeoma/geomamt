      DOUBLE PRECISION FUNCTION G01CEF(P,IFAIL)
C     MARK 10 RE-ISSUE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12A REVISED. IER-511 (AUG 1986).
C
C     G01CEF RETURNS THE DEVIATE ASSOCIATED WITH THE LOWER
C     TAIL PROBABILITY P FROM THE STANDARD NORMAL
C     DISTRIBUTION.
C
C     WRITTEN BY N.M.MACLAREN AND COLLEAGUES
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     RELATIVE ACCURACY IS 5.0E-13
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01CEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 P
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 P1, P2, P3, P4, P5, P6, P7, P8,
     *                                 Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8,
     *                                 R2PI, S, T2B1, T2B2, T2B3, T2B4,
     *                                 T2B5, T2B6, T2T1, T2T2, T2T3,
     *                                 T2T4, T2T5, T2T6, T3B1, T3B2,
     *                                 T3B3, T3B4, T3B5, T3B6, T3B7,
     *                                 T3T1, T3T2, T3T3, T3T4, T3T5,
     *                                 T3T6, T3T7, X, Y, YD
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 X01AAF
      INTEGER                          P01ABF
      EXTERNAL                         X01AAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MIN, SIGN, LOG, SQRT
C     .. Data statements ..
      DATA                             P1, P2, P3, P4, P5, P6, P7,
     *                                 P8/0.1000000000000000D+1,
     *                                 -0.4163362762616374D+2,
     *                                 0.6494128979404664D+2,
     *                                 -0.3689355386573687D+2,
     *                                 0.8984374887949291D+1,
     *                                 -0.8486205099916682D+0,
     *                                 0.1799227452322391D-1,
     *                                 0.1243484778425483D-3/
      DATA                             Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *                                 Q8/0.1000000000000000D+1,
     *                                 -0.4180029429283041D+2,
     *                                 0.7184967217618504D+2,
     *                                 -0.4645534714071778D+2,
     *                                 0.1357696314152835D+2,
     *                                 -0.1714887182301168D+1,
     *                                 0.6929108687616727D-1, 0.15D-8/
      DATA                             T2T1, T2T2, T2T3, T2T4, T2T5,
     *                                 T2T6/8.416212335733377D-1,
     *                                 -6.002063030711376D0,
     *                                 -4.722617278101355D-1,
     *                                 68.25030402764662D0,
     *                                 -96.23966784668191D0,
     *                                 -10.47655255702662D0/
      DATA                             T2B1, T2B2, T2B3, T2B4, T2B5,
     *                                 T2B6/1.0D0, -11.37563668160269D0,
     *                                 41.33878042998575D0,
     *                                 -43.59278882675467D0,
     *                                 -21.25589975172773D0,
     *                                 25.08610603077810D0/
      DATA                             T3T1, T3T2, T3T3, T3T4, T3T5,
     *                                 T3T6, T3T7/-3.012472607913056D0,
     *                                 2.982419254309647D0,
     *                                 22.70676420727861D0,
     *                                 7.670763790822552D0,
     *                                 5.519186629029667D-1,
     *                                 7.985443046076538D-3,
     *                                 4.150184151574655D-6/
      DATA                             T3B1, T3B2, T3B3, T3B4, T3B5,
     *                                 T3B6, T3B7/1.0D0,
     *                                 1.239758125817922D0,
     *                                 -1.205350980555889D1,
     *                                 -12.02359058219926D0,
     *                                 -2.373035796200643D0,
     *                                 -1.193440282031508D-1,
     *                                 -1.216250189896074D-3/
C     .. Executable Statements ..
      IF (P.LE.0.0D0 .OR. P.GE.1.0D0) GO TO 60
      IFAIL = 0
      X = P - 0.5D0
      IF (ABS(X).GT.0.3D0) GO TO 20
C
C     BREAK UP RANGE
C     FOR 0.2 LE P LE 0.8 WE USE A RATIONAL TCHEBYCHEV (P/Q)
C
      R2PI = SQRT(2.0D0*X01AAF(0.0D0))
      S = X*R2PI
      X = S*S
      Y = P1 + X*(P2+X*(P3+X*(P4+X*(P5+X*(P6+X*(P7+X*P8))))))
      YD = Q1 + X*(Q2+X*(Q3+X*(Q4+X*(Q5+X*(Q6+X*(Q7+X*Q8))))))
      Y = Y/YD
      G01CEF = Y*S
      RETURN
C
C     FOR 0.08 LE P LE 0.2 OR 0.8 LE P LE 0.92
C     WE USE RATIONAL TCHEBYCHEV (T2T/T2B)
C
   20 S = SIGN(1.0D0,X)
      IF (ABS(X).GT.0.42D0) GO TO 40
C
C     BREAK UP RANGE SOME MORE
C
      X = ABS(X) - 0.3D0
      Y = T2T1 + X*(T2T2+X*(T2T3+X*(T2T4+X*(T2T5+X*T2T6))))
      YD = T2B1 + X*(T2B2+X*(T2B3+X*(T2B4+X*(T2B5+X*T2B6))))
      Y = Y/YD
      G01CEF = Y*S
      RETURN
C
C     THE FOLLOWING CASE HANDLES ASYMPTOTIC BEHAVIOUR.
C     UNFORTUNATELY WE MUST MAKE A TRANSFORMATION.
C
   40 X = SQRT((-2.0D0)*LOG(MIN(P,1.0D0-P)))
      Y = T3T1 + X*(T3T2+X*(T3T3+X*(T3T4+X*(T3T5+X*(T3T6+X*T3T7)))))
      YD = T3B1 + X*(T3B2+X*(T3B3+X*(T3B4+X*(T3B5+X*(T3B6+X*T3B7)))))
      Y = Y/YD + X
      G01CEF = S*Y
      RETURN
C
C     ERROR EXITS
C
   60 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      G01CEF = 0.0D0
      RETURN
      END
