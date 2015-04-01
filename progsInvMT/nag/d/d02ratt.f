      SUBROUTINE D02RAT(A,C,DEL,M,N,P,R,IR,IC,SING,ICA,AUX,MTNMAX,NMAX,
     *                  MMAX2)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ------------------------------------------------------
C
C     LU  DECOMPOSITION OF THE BLOCK TRIDIAGONAL MATRIX
C
C
C     .            I  A1(1)   C1(1)     0   ............   0  I
C     .            I                                          I
C     .            I  B1(2)   A1(2)   C1(2)     0    ...   0  I
C     .            I                                          I
C     .       AA = I    .     .     .     .      .      .     I
C     .            I                                          I
C     .            I    0       .   B1(N-1)  A1(N-1)  C1(N-1) I
C     .            I                                          I
C     .            I  D(1)      0      0      B1(N)    A1(N)  I
C
C     WHERE THE BLOCKS  B1,A1,C1,D  ARE  M*M,  B1(J)  HAS
C     NONZERO ELEMENTS ONLY IN ITS FIRST  P  ROWS, C1(J)
C     HAS NONZERO ELEMENTS ONLY IN ITS LAST  (M-P)  ROWS,
C     WHILE  D(1)  HAS NONZERO ELEMENTS IN THE ROWS  (P+1),
C     . . .,(P+R) .
C     THIS NONZERO INFORMATION IS PACKED IN THE TWO -DIMENSIONAL
C     ARRAYS  A,C,DEL, IN THE FOLLOWING WAY
C
C     .            I   A1(1)   I                 I I B1(2) I I
C     .            I           I                 I I C1(1)
C     .            I           I                 I           I
C     .            I   A1(2)   I                 I I B1(3) I I
C     .            I           I                 I I C1(2) I I
C     .        A = I     .     I    ,        C = I     .     I
C     .            I           I                 I           I
C     .            I   A1(N)   I                 I I B1(N) I I
C     .            I           I (M*N)*M         I IC1(N-1)I I
C
C
C     DEL = (   D(1) , 0 , ... , 0 )    (R*(M*N))
C
C
C     THE  LU  DECOMPOSITION HAS THE FORM
C
C
C     .            I    I      0      .     .     .     0  I
C     .            I                                       I
C     .            I  BE(2)    I     0      .     .     0  I
C     .            I                                       I
C     .        L = I                                       I
C     .            I                                       I
C     .            I    0      0       BE(N-1)   I      0  I
C     .            I                                       I
C     .            I  DE(1)  DE(2)    .    .   BE(N)    I  I
C
C
C     .            I  AL(1)  GA(1)     0      .         0  I
C     .            I                                       I
C     .            I     0    AL(2)  GA(2)    0         0  I
C     .            I                                       I
C     .        U = I                                       I
C     .            I                                       I
C     .            I                       AL(N-1)  GA(N-1)I
C     .            I                                       I
C     .            I                          0      AL(N) I
C
C
C     WHERE  BE(J),GA(J),DE(1), HAVE THE SAME ZERO PATTERN
C     AS  B1(J),C1(J),D(1), WITH THE EXCEPTION OF BE(N)
C     WHICH HAS ITS ROWS  (P+1),...,(P+R)  FILLED IN.
C     THOSE SAME ROWS ARE ALSO NONZERO IN  DE(J).
C     UPON OUTPUT, THIS INFORMATION OVERWRITES THE DATA ARRAYS
C     A,C,DEL.
C     THE  AL(J)  ARE STORED ALSO IN  LU  FACTORED FORM,AND
C     I  IS THE  M*M  IDENTITY MATRIX.
C     WE SHOULD ALSO POINT OUT THAT THE FIRST  P  ROWS OF
C     BE(N)  ARE STORED IN THE FIRST  P  ROWS OF  C(N-1),
C     WHILE ITS FOLLOWING  R  ROWS ARE STORED IN DEL(N-1).
C     DEL(N)  IS USED AS AN AUXILIARY WORKING AREA.
C     IN REFERENCE  (1)  IT IS SHOWN THAT IF  AA  IS NONSINGULAR
C     THEN BY USING A PARTIAL PIVOTING STRATEGY THAT DOES
C     NOT DESTROY THE SPARSENESS, THE  LU  DECOMPOSITION
C     IS FEASIBLE. IN THE  KTH  STEP OF THIS
C     PROCESS,COLUMN PIVOTING IS USED IN ORDER TO
C     TRIANGULARIZE THE FIRST  P  ROWS OF A1(K). THEN, ROW
C     INTERCHANGES ARE USED BETWEEN THE LAST (M-P) ROW OF
C     A1(K) AND THE  P  NONZERO ROWS OF  B1(K+1) IN ORDER
C     TO COMPLETE THE TRIANGULARIZATION OF  A1(K). THE
C     LAST STEP PRODUCES AN  R  ROW FILL IN IN THE LAST
C     ROW BLOCK OF  L.
C     IF EITHER  P  OR  R  ARE ZERO THEN THERE ARE SIMPLIFICATIONS.
C     THE INTERCHANGES ARE PERFORMED PHYSICALLY ON  AA ,
C     AND THE INFORMATION ABOUT THEM IS STORED IN THE INTEGER
C     ARRAYS  IR,IC.  UPON OUTPUT, THE CONTENTS OF
C     IC(I,J)  TELLS US WHICH ORIGINAL COLUMN IN BLOCK  I
C     OCCUPIES NOW POSITION J. WHEN SOLVING A SYSTEM OF
C     EQUATIONS WITH THE MATRIX AA USING THIS LU
C     DECOMPOSITION  THE FINAL ANSWERS SHOULD BE PERMUTED
C     ACCORDINGLY.
C     DUE TO THE SPECIAL WAY IN WHICH THE ROW INTERCHANGES
C     ARE PERFORMED THE ORIGIN OF THE BLOCK COUNTING IN  IR
C     IS SHIFTED IN  P  ROWS (I.E., THE FIRST  P  ROWS OF  AA
C     ARE IGNORED, AND THE FIRST BLOCK CONTAINS THE  (M-P)  LAST
C     ROWS OF A1(1)  AND THE  P  NONZERO ROWS OF  B1(2)).
C     OTHERWISE  THE INFORMATION IS STORED ,AND HAS THE SAME
C     MEANING AS FOR  IC, CHANGING COLUMN BY ROWS.  IN SOLVING
C     A SYSTEM OF EQUATIONS, THE RIGHT HAND SIDE SHOULD BE
C     PERMUTED ACCORDINGLY BEFORE PROCESSING.
C     THE COMPANION SUBROUTINE SOLVE , TAKES THE INFORMATION
C     PROVIDED BY  DECOMP  AND SOLVES THE SYSTEM OF LINEAR
C     EQUATIONS
C
C     .           AA * X = B  .
C
C     IF AT ANY TIME A PIVOTAL ELEMENT TURNS OUT TO BE ZERO
C     THE PROCESS IS INTERRUPTED AND THE LOGICAL VARIABLE
C     SING  IS SET EQUAL TO .TRUE. .
C
C     (1)  H.B.KELLER, ACCURATE DIFFERENCE METHODS FOR NON-
C     .            LINEAR TWO-POINT BOUNDARY VALUE PROBLEMS.
C     .            SIAM J. NUMER. ANAL. 11, PP. 305-320 (1974).
C
C     ------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           M, MMAX2, MTNMAX, N, NMAX, P, R
      LOGICAL           SING
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MTNMAX,M), AUX(MMAX2,M), C(MTNMAX,M),
     *                  DEL(M,MTNMAX)
      INTEGER           IC(NMAX,M), ICA(M), IR(NMAX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  AK, TE, XMU, Y
      INTEGER           I, I0, I01, I1, I1S, I2, I2PJS, I3, I4, I7, I8,
     *                  II, II1, II2, INI, IP, IPIV, IPM, IPP, IS, IS1,
     *                  ISM, ISP, ITE, IX, IXI, IXM, IXMPJS, IXP, IXPI,
     *                  IXPJS, IXPKS, J, J1, JMIX, JP, JS, JSX, JX, KS,
     *                  M1, MP, N1, P1
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      SING = .FALSE.
      P1 = P + 1
      MP = M + P
      N1 = N - 1
C
C     MAIN LOOP
C
      DO 1100 II = 1, N1
         IX = (II-1)*M
         IXM = IX + M
         IXI = IX - M
         DO 40 I = 1, M
            IC(II,I) = I
            IR(II,I) = I
            IXPI = IX + I
            DO 20 J = 1, M
               AUX(I,J) = A(IXPI,J)
   20       CONTINUE
   40    CONTINUE
         IF (P) 260, 260, 60
   60    DO 100 I = 1, P
            IXPI = IX + I
            IPM = I + M
            DO 80 J = 1, M
               AUX(IPM,J) = C(IXPI,J)
   80       CONTINUE
  100    CONTINUE
C
C        REDUCTION OF FIRST P ROWS OF A1(II)
C
         DO 240 I = 1, P
            TE = 0.D0
            DO 120 J = I, M
               Y = ABS(AUX(I,J))
               IF (Y.LE.TE) GO TO 120
               TE = Y
               IPIV = J
  120       CONTINUE
            IF (TE.EQ.0.D0) GO TO 1340
C
C           COLUMN INTERCHANGES
C
            IF (IPIV.EQ.I) GO TO 180
            ITE = IC(II,I)
            IC(II,I) = IC(II,IPIV)
            IC(II,IPIV) = ITE
            DO 140 IS = 1, MP
               TE = AUX(IS,I)
               AUX(IS,I) = AUX(IS,IPIV)
               AUX(IS,IPIV) = TE
  140       CONTINUE
            IF (II.EQ.1) GO TO 180
            I1 = IXI + P1
            I2 = IXI + M
            DO 160 IS = I1, I2
               TE = C(IS,I)
               C(IS,I) = C(IS,IPIV)
               C(IS,IPIV) = TE
  160       CONTINUE
C
C           FACTORIZATION
C
  180       AK = 1.D0/AUX(I,I)
            I1 = I + 1
            DO 220 IS = I1, MP
               XMU = AUX(IS,I)*AK
               AUX(IS,I) = XMU
               DO 200 JS = I1, M
                  AUX(IS,JS) = AUX(IS,JS) - XMU*AUX(I,JS)
  200          CONTINUE
  220       CONTINUE
  240    CONTINUE
C
C        REDUCTION OF REMAINING PART OF A1(II)
C
  260    DO 440 J = P1, M
            TE = 0.D0
            DO 280 I = J, MP
               Y = ABS(AUX(I,J))
               IF (Y.LE.TE) GO TO 280
               TE = Y
               IPP = I
  280       CONTINUE
            IF (TE.EQ.0.D0) GO TO 1340
C
C           ROW INTERCHANGES
C
            IF (IPP.EQ.J) GO TO 380
            IPIV = IPP - P
            JP = J - P
            ITE = IR(II,JP)
            IR(II,JP) = IR(II,IPIV)
            IR(II,IPIV) = ITE
            DO 300 JS = 1, M
               TE = AUX(J,JS)
               AUX(J,JS) = AUX(IPP,JS)
               AUX(IPP,JS) = TE
  300       CONTINUE
            II1 = IX + J
            II2 = IX + IPP
            IF (IPP.GT.M) GO TO 340
            DO 320 JS = 1, M
               TE = C(II1,JS)
               C(II1,JS) = C(II2,JS)
               C(II2,JS) = TE
  320       CONTINUE
            GO TO 380
  340       DO 360 JS = 1, M
               TE = C(II1,JS)
               C(II1,JS) = A(II2,JS)
               A(II2,JS) = TE
  360       CONTINUE
C
C           FACTORIZATION
C
  380       IF (J.EQ.M) GO TO 440
            AK = 1.D0/AUX(J,J)
            J1 = J + 1
            DO 420 IS = J1, MP
               XMU = AUX(IS,J)*AK
               AUX(IS,J) = XMU
               DO 400 JS = J1, M
                  AUX(IS,JS) = AUX(IS,JS) - XMU*AUX(J,JS)
  400          CONTINUE
  420       CONTINUE
  440    CONTINUE
C
C        RESTORING B1(II+1) AND AL(II)
C
         IF (P.EQ.0) GO TO 500
         I1 = M - P + 1
         DO 480 IS = I1, M
            IP = IR(II,IS)
            IF (IP.EQ.IS) GO TO 480
            ISP = IS + P + IXI
            IPP = IP + P + IX
            DO 460 JS = 1, M
               C(ISP,JS) = A(IPP,JS)
  460       CONTINUE
  480    CONTINUE
  500    DO 540 I = 1, M
            IXPI = IX + I
            DO 520 J = 1, M
               A(IXPI,J) = AUX(I,J)
  520       CONTINUE
  540    CONTINUE
C
         IS = IX + 1
         ISP = IX + P
         DO 560 JS = 1, M
            ICA(JS) = IC(II,JS)
  560    CONTINUE
         JS = 1
  580    IF (JS.GE.M) GO TO 720
         IP = ICA(JS)
         IF (IP.NE.JS) GO TO 600
         JS = JS + 1
         GO TO 580
C
  600    INI = JS
         ICA(INI) = INI
C
C        COLUMN INTERCHANGES FOR B1(II+1)
C
         IF (P.EQ.0) GO TO 660
  620    DO 640 I1 = IS, ISP
            TE = C(I1,INI)
            C(I1,INI) = C(I1,IP)
            C(I1,IP) = TE
  640    CONTINUE
  660    IF (R.EQ.0) GO TO 700
C
C        COLUMN INTERCHANGES OF D(II)
C
         IPP = IX + IP
         JX = IX + INI
         DO 680 IS1 = 1, R
            TE = DEL(IS1,IPP)
            DEL(IS1,IPP) = DEL(IS1,JX)
            DEL(IS1,JX) = TE
  680    CONTINUE
  700    INI = IP
         IP = ICA(IP)
         ICA(INI) = INI
         IF (IP.NE.JS) GO TO 620
         JS = JS + 1
         GO TO 580
C
  720    IF (R.EQ.0) GO TO 920
C
C        SOLUTION OF DE(II)*AL(II)=- DELTA(II)
C
         I1 = IX + 1
         DO 840 IS = 1, R
            DO 780 JS = I1, IXM
               TE = DEL(IS,JS)
               JSX = JS - IX
               IF (JS.EQ.I1) GO TO 760
               J1 = JS - 1
               DO 740 KS = I1, J1
                  TE = TE - DEL(IS,KS)*A(KS,JSX)
  740          CONTINUE
  760          DEL(IS,JS) = TE/A(JS,JSX)
  780       CONTINUE
            M1 = M - 1
            DO 820 JS = 1, M1
               J = IXM - JS
               JMIX = J - IX
               TE = DEL(IS,J)
               J1 = J + 1
               DO 800 KS = J1, IXM
                  TE = TE - DEL(IS,KS)*A(KS,JMIX)
  800          CONTINUE
               DEL(IS,J) = TE
  820       CONTINUE
  840    CONTINUE
C
C        COMPUTATION OF DE(II+1)=-DE(II)*GA(II)
C
         IXP = IX + P1
         DO 900 IS = 1, R
            DO 880 JS = 1, M
               TE = 0.D0
               DO 860 KS = IXP, IXM
                  TE = TE - DEL(IS,KS)*C(KS,JS)
  860          CONTINUE
               IXMPJS = IXM + JS
               DEL(IS,IXMPJS) = TE
  880       CONTINUE
  900    CONTINUE
C
C        SOLUTION OF BE(II+1)*AL(II)=B(II+1)
C
  920    I1 = IX + 1
         IF (P.EQ.0) GO TO 1100
         I2 = IX + P
         DO 1080 IS = I1, I2
            DO 980 JS = 1, M
               TE = C(IS,JS)
               IXPJS = IX + JS
               IF (JS.EQ.1) GO TO 960
               J1 = JS - 1
               DO 940 KS = 1, J1
                  IXPKS = IX + KS
                  TE = TE - C(IS,KS)*A(IXPKS,JS)
  940          CONTINUE
  960          C(IS,JS) = TE/A(IXPJS,JS)
  980       CONTINUE
            M1 = M - 1
            DO 1020 JS = 1, M1
               J = M - JS
               TE = C(IS,J)
               J1 = J + 1
               DO 1000 KS = J1, M
                  IXPKS = IX + KS
                  TE = TE - C(IS,KS)*A(IXPKS,J)
 1000          CONTINUE
               C(IS,J) = TE
 1020       CONTINUE
C
C           COMPUTATION OF AL(II+1)=A(II+1)-BE(II+1)*GA(II)
C
            ISM = IS + M
            DO 1060 JS = 1, M
               TE = A(ISM,JS)
               DO 1040 KS = P1, M
                  IXPKS = IX + KS
                  TE = TE - C(IS,KS)*C(IXPKS,JS)
 1040          CONTINUE
               A(ISM,JS) = TE
 1060       CONTINUE
 1080    CONTINUE
 1100 CONTINUE
C
C     COMPLETE COMPUTATION OF ALFA(N)
C
      I2 = N1*M
      I1 = I2 + P
      IF (R.EQ.0) GO TO 1160
      DO 1140 IS = 1, R
         I1S = I1 + IS
         DO 1120 JS = 1, M
            I2PJS = I2 + JS
            A(I1S,JS) = A(I1S,JS) + DEL(IS,I2PJS)
 1120    CONTINUE
 1140 CONTINUE
C
C     L U DECOMPOSITION OF ALFA(N) (COLUMN PIVOTING)
C
 1160 I4 = I1 - M + 1
      I1 = I2 + 1
      I3 = N*M - 1
      I8 = I3 + 1
      DO 1180 I = 1, M
         IC(N,I) = I
 1180 CONTINUE
      DO 1320 I = I1, I3
         I0 = I - I2
         TE = 0.D0
         DO 1200 J = I0, M
            Y = ABS(A(I,J))
            IF (Y.LE.TE) GO TO 1200
            TE = Y
            IPIV = J
 1200    CONTINUE
         IF (TE.EQ.0.D0) GO TO 1340
C
C        COLUMN INTERCHANGES
C
         IF (IPIV.EQ.I0) GO TO 1260
         ITE = IC(N,I0)
         IC(N,I0) = IC(N,IPIV)
         IC(N,IPIV) = ITE
         DO 1220 IS = I1, I8
            TE = A(IS,I0)
            A(IS,I0) = A(IS,IPIV)
            A(IS,IPIV) = TE
 1220    CONTINUE
         DO 1240 IS = I4, I2
            TE = C(IS,I0)
            C(IS,I0) = C(IS,IPIV)
            C(IS,IPIV) = TE
 1240    CONTINUE
C
C        FACTORIZATION
C
 1260    AK = 1.D0/A(I,I0)
         I01 = I0 + 1
         I7 = I + 1
         DO 1300 IS = I7, I8
            XMU = A(IS,I0)*AK
            A(IS,I0) = XMU
            DO 1280 JS = I01, M
               A(IS,JS) = A(IS,JS) - XMU*A(I,JS)
 1280       CONTINUE
 1300    CONTINUE
 1320 CONTINUE
C
      RETURN
 1340 SING = .TRUE.
      RETURN
      END
