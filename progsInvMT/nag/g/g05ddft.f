      DOUBLE PRECISION FUNCTION G05DDF(A,B)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11 REVISED. IER-441 (FEB 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER, NORMALLY DISTRIBUTED (GAUSSIAN)
C     WITH MEAN A AND STANDARD DEVIATION B.
C     THE METHOD USED IS A MODIFICATION OF CACM ALGORITHM 488
C     (BY R.P.BRENT).
C     THE FOLLOWING CAN BE EXTENDED IF A TRUNCATION PROBABILITY
C     OF 1.0E-12 IS REGARDED AS UNSATISFACTORY.
C     THE CODE ASSUMES THAT THE CONTENTS OF COMMON BLOCK /BG05CA/
C     ARE SAVED BETWEEN CALLS OF THIS ROUTINE AND CALLS OF OTHER
C     G05 ROUTINES THAT REFERENCE IT. A SAVE STATEMENT
C     ENSURES THAT THIS IS SO.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
C     .. Scalars in Common ..
      DOUBLE PRECISION                 STORE1, STORE2
C     .. Local Scalars ..
      DOUBLE PRECISION                 HALF, ONE, T, U, V, W
      INTEGER                          N
C     .. Local Arrays ..
      DOUBLE PRECISION                 D(41)
C     .. External Functions ..
      DOUBLE PRECISION                 G05CAF
      EXTERNAL                         G05CAF
C     .. Common blocks ..
      COMMON                           /BG05CA/STORE1, STORE2
C     .. Save statement ..
      SAVE                             /BG05CA/
C     .. Data statements ..
      DATA                             ONE/1.0D0/, HALF/0.5D0/
      DATA                             D(1), D(2), D(3), D(4), D(5),
     *                                 D(6), D(7), D(8), D(9), D(10),
     *                                 D(11), D(12), D(13), D(14)/0.0D0,
     *                                 0.674489750196082D0,
     *                                 1.150349380376008D0,
     *                                 1.534120544352546D0,
     *                                 1.862731867421652D0,
     *                                 2.153874694061456D0,
     *                                 2.417559016236505D0,
     *                                 2.660067468617460D0,
     *                                 2.885634912426757D0,
     *                                 3.097269078198785D0,
     *                                 3.297193345691964D0,
     *                                 3.487104104114431D0,
     *                                 3.668329285121323D0,
     *                                 3.841930685501911D0/
      DATA                             D(15), D(16), D(17), D(18),
     *                                 D(19), D(20), D(21), D(22),
     *                                 D(23), D(24), D(25), D(26),
     *                                 D(27)/4.008772594168585D0,
     *                                 4.169569323349106D0,
     *                                 4.324919040826046D0,
     *                                 4.475328424654204D0,
     *                                 4.621231001499247D0,
     *                                 4.763001034267814D0,
     *                                 4.900964207963193D0,
     *                                 5.035405969463927D0,
     *                                 5.166578119728753D0,
     *                                 5.294704084854598D0,
     *                                 5.419983174916868D0,
     *                                 5.542594057802940D0,
     *                                 5.662697617459439D0/
      DATA                             D(28), D(29), D(30), D(31),
     *                                 D(32), D(33), D(34), D(35),
     *                                 D(36), D(37), D(38), D(39),
     *                                 D(40)/5.780439324478935D0,
     *                                 5.895951216739571D0,
     *                                 6.009353565530745D0,
     *                                 6.120756285971941D0,
     *                                 6.230260137989044D0,
     *                                 6.337957754553790D0,
     *                                 6.443934526538564D0,
     *                                 6.548269367831731D0,
     *                                 6.651035379893011D0,
     *                                 6.752300431407015D0,
     *                                 6.852127665896068D0,
     *                                 6.950575947916750D0,
     *                                 7.047700256664409D0/
      DATA                             D(41)/7.143552034352190D0/
C     .. Executable Statements ..
C     FIRST CALL TO G05CAF COMES BEFORE FIRST REFERENCE TO STORE1
C     TO MAKE SURE THAT STORE1 IS INITIALIZED
      V = G05CAF(0.0D0)
      U = STORE1
      DO 20 N = 1, 39
         IF (U.GT.HALF) GO TO 40
         U = U + U
   20 CONTINUE
      N = 40
   40 T = D(N)
      U = V
   60 W = (D(N+1)-T)*U
      V = W*(W*HALF+T)
   80 U = G05CAF(0.0D0)
      IF (V.LE.U) GO TO 100
      V = G05CAF(0.0D0)
      IF (U.GT.V) GO TO 80
      U = (V-U)/(ONE-U)
      GO TO 60
  100 U = (U-V)/(ONE-V)
      IF (U.GT.HALF) GO TO 120
      STORE1 = U + U
      G05DDF = A + B*(W+T)
      RETURN
  120 STORE1 = U + U - ONE
      G05DDF = A - B*(W+T)
      RETURN
      END
