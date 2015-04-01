      SUBROUTINE D02RAU(M,N,P,R,X,Y,IY,FCN,G,FCNEP,JACOBG,JACOBE,B20,
     *                  FCNA,A1,B1,A,C,DEL,CASI,SING,IR,IC,UU,RES,LIN,
     *                  MTNMAX,NMAX,MMAX2,HX,GRADF,AUX,ICA,XAU,EPSNU,
     *                  HMAX,IG,NIG,H,IRN,IP1,NIP,F,ALPHA,IFLAG,LP)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9C REVISED. IER-367 (MAY 1982).
C     MARK 10B REVISED. IER-401 (JAN 1983).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     FCN, FCNA, FCNEP, G, JACOBE, JACOBG
C
C     CASI = .TRUE.  NO DECOMPOSITION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPSNU
      INTEGER           IFLAG, IY, LIN, LP, M, MMAX2, MTNMAX, N, NIG,
     *                  NIP, NMAX, P, R
      LOGICAL           CASI, SING
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MTNMAX,M), A1(M,M), ALPHA(M), AUX(MMAX2,M),
     *                  B1(M,M), B20(M,2), C(MTNMAX,M), DEL(M,MTNMAX),
     *                  F(MTNMAX), GRADF(MTNMAX), H(M), HMAX(M),
     *                  HX(NMAX), RES(MTNMAX), UU(MTNMAX), X(NMAX),
     *                  XAU(MTNMAX), Y(IY,N)
      INTEGER           IC(NMAX,M), ICA(M), IG(NIG), IP1(NIP),
     *                  IR(NMAX,M), IRN(M,M)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, FCNA, FCNEP, G, JACOBE, JACOBG
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSMCH, H1, TE
      INTEGER           I, I1, I1J, I1PJ, I1PJJ, I2, I2J, I2PIT, I3,
     *                  I3P, I3PIT, I3PJJ, I3PPI, I4, I5, IFAIL, II, IM,
     *                  IND, IP, IPP, IT, J, J1, JJ, JJPJ, K, M1, MMP,
     *                  N1, NERR, P1, P2, P3, PM
      CHARACTER*40      REC
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          D02RAR, D02RAS, D02RAT, X04AAF, X04BAF
C     .. Executable Statements ..
      IF (CASI) GO TO 1200
      IF (LP.NE.0) CALL X04AAF(0,NERR)
      EPSMCH = X02AJF()
      P1 = P - 1
      N1 = N - 1
      M1 = M - P
      MMP = M1 - 1
      PM = P + 1
C
C     CONSTRUCTION OF A,C,DEL
C
      DO 660 I = 1, N
         JJ = (I-1)*M
         GO TO (580,420,20,40) LIN
   20    CALL JACOBE(X(I),EPSNU,Y(1,I),DEL,M)
         IF (I.NE.1) GO TO 600
         CALL JACOBG(EPSNU,Y(1,1),Y(1,N),A1,B1,M)
         GO TO 220
   40    IND = 1
         IFAIL = 1
   60    CALL D02RAR(M,M,IRN,M*M,IP1,NIP,H,Y(1,I),F(JJ+1)
     *               ,HMAX,10.0D0,100.0D0,1000.0D0,EPSMCH,EPSMCH,.TRUE.,
     *               DEL,IG,NIG,A,M+M,XAU,IND,IFAIL)
         IF (IND.EQ.0) GO TO 80
         CALL FCNEP(X(I),EPSNU,Y(1,I),A,M)
         GO TO 60
   80    IF (IFAIL.EQ.0) GO TO 100
C        IMPOSS EXIT FROM D02RAR
         IFLAG = 6
         RETURN
  100    IF (I.NE.1) GO TO 600
         IND = 1
         IFAIL = 1
  120    CALL D02RAR(M,M,IRN,M*M,IP1,NIP,H,Y(1,1)
     *               ,ALPHA,HMAX,10.0D0,100.0D0,1000.0D0,EPSMCH,EPSMCH,
     *               .TRUE.,A1,IG,NIG,A,M+M,XAU,IND,IFAIL)
         IF (IND.EQ.0) GO TO 140
         CALL G(EPSNU,Y(1,1),Y(1,N),A,M)
         GO TO 120
  140    IF (IFAIL.EQ.0) GO TO 160
         IFLAG = 6
         RETURN
  160    IND = 1
         IFAIL = 1
  180    CALL D02RAR(M,M,IRN,M*M,IP1,NIP,H,Y(1,N)
     *               ,ALPHA,HMAX,10.0D0,100.0D0,1000.0D0,EPSMCH,EPSMCH,
     *               .TRUE.,B1,IG,NIG,A,M+M,XAU,IND,IFAIL)
         IF (IND.EQ.0) GO TO 200
         CALL G(EPSNU,Y(1,1),Y(1,N),A,M)
         GO TO 180
  200    IF (IFAIL.EQ.0) GO TO 220
C        IMPOSSIBLE EXIT FROM D02RAR
         IFLAG = 6
         RETURN
  220    IF (P.EQ.0) GO TO 300
         DO 260 I5 = 1, P
            K = 0
            DO 240 J = 1, M
               IF (A1(I5,J).NE.0.0D0) K = 1
  240       CONTINUE
            IF (K.EQ.0) GO TO 280
  260    CONTINUE
         GO TO 300
C        JACOBIAN SINGULAR
  280    IFLAG = 8
         IF (LP.NE.0) THEN
            CALL X04BAF(NERR,' ')
            WRITE (REC,FMT=99999)
            CALL X04BAF(NERR,REC)
         END IF
         RETURN
  300    P3 = P + R
         IF (R.EQ.0) GO TO 360
         P2 = P + 1
         DO 340 I5 = P2, P3
            K = 0
            DO 320 J = 1, M
               IF (A1(I5,J).NE.0.0D0) K = 1
               IF (B1(I5,J).NE.0.0D0) K = 1
  320       CONTINUE
            IF (K.EQ.0) GO TO 280
  340    CONTINUE
         IF (P3.EQ.M) GO TO 600
  360    P3 = P3 + 1
         DO 400 I5 = P3, M
            K = 0
            DO 380 J = 1, M
               IF (B1(I5,J).NE.0.0D0) K = 1
  380       CONTINUE
            IF (K.EQ.0) GO TO 280
  400    CONTINUE
         GO TO 600
  420    IND = 1
         IFAIL = 1
  440    CALL D02RAR(M,M,IRN,M*M,IP1,NIP,H,Y(1,I),F(JJ+1)
     *               ,HMAX,10.0D0,100.0D0,1000.0D0,EPSMCH,EPSMCH,.TRUE.,
     *               DEL,IG,NIG,A,M+M,XAU,IND,IFAIL)
         IF (IND.EQ.0) GO TO 460
         CALL FCN(X(I),Y(1,I),A)
         GO TO 440
  460    IF (IFAIL.EQ.0) GO TO 480
         IFLAG = 6
         RETURN
  480    IF (I.NE.1) GO TO 600
         DO 520 J = 1, M
            DO 500 J1 = 1, M
               A1(J,J1) = 0.0D0
               B1(J,J1) = 0.0D0
  500       CONTINUE
  520    CONTINUE
         J = 0
         DO 540 J1 = 1, M
            IF (B20(J1,1).NE.0.0D0) GO TO 540
            J = J + 1
            A1(J,J1) = 1.0D0
  540    CONTINUE
         DO 560 J1 = 1, M
            IF (B20(J1,2).NE.0.0D0) GO TO 560
            J = J + 1
            B1(J,J1) = 1.0D0
  560    CONTINUE
         GO TO 600
  580    CALL FCNA(X(I),DEL)
  600    CONTINUE
         DO 640 J = 1, M
            JJPJ = JJ + J
            DO 620 K = 1, M
               C(JJPJ,K) = DEL(J,K)
  620       CONTINUE
  640    CONTINUE
  660 CONTINUE
      IF (P.EQ.0) GO TO 720
      DO 700 I = 1, P
         DO 680 J = 1, M
            A(I,J) = A1(I,J)
  680    CONTINUE
  700 CONTINUE
  720 DO 880 II = 1, N1
         H1 = .5D0*HX(II)
         I1 = (II-1)*M + 1
         I2 = I1 + MMP
         IT = 0
         DO 760 I = I1, I2
            IP = I + P
            IT = IT + 1
            DO 740 J = 1, M
               A(IP,J) = -H1*C(I,J)
               IF (IT.EQ.J) A(IP,J) = A(IP,J) - 1.D0
  740       CONTINUE
  760    CONTINUE
         I3 = I2 + M
         I4 = I1 + P1
         IT = 0
         IF (P.EQ.0) GO TO 820
         DO 800 I = I1, I4
            IT = IT + 1
            I3PIT = I3 + IT
            I2PIT = I2 + IT
            DO 780 J = 1, M
               IM = I + M
               A(IM,J) = -H1*C(I3PIT,J)
               C(I,J) = -H1*C(I2PIT,J)
               IF (IT+M1.NE.J) GO TO 780
               A(IM,J) = A(IM,J) + 1.D0
               C(I,J) = C(I,J) - 1.D0
  780       CONTINUE
  800    CONTINUE
  820    IT = 0
         DO 860 I = I1, I2
            IP = I + P
            IT = IT + 1
            IM = I + M
            DO 840 J = 1, M
               C(IP,J) = -H1*C(IM,J)
               IF (IT.EQ.J) C(IP,J) = C(IP,J) + 1.D0
  840       CONTINUE
  860    CONTINUE
  880 CONTINUE
C
      I1 = N1*M + PM
      I2 = N*M
      IT = P
      DO 920 I = I1, I2
         IT = IT + 1
         DO 900 J = 1, M
            A(I,J) = B1(IT,J)
  900    CONTINUE
  920 CONTINUE
      IF (R.EQ.0) GO TO 980
      DO 960 I = 1, R
         IPP = I + P
         DO 940 J = 1, M
            DEL(I,J) = A1(IPP,J)
  940    CONTINUE
  960 CONTINUE
  980 IF (LIN.EQ.1) GO TO 1180
C     COMPUTATION  OF  GRADIENT
      IT = 0
      DO 1120 I = 1, N
         I1 = (I-1)*M
         I2 = I1 - M1
         I3 = I*M
         DO 1100 J = 1, M
            IT = IT + 1
            TE = 0.D0
            I1PJ = I1 + J
            DO 1000 JJ = 1, M
               I1J = I1 + JJ
               TE = TE + A(I1J,J)*RES(I1J)
 1000       CONTINUE
            IF (I.EQ.N) GO TO 1040
            DO 1020 JJ = 1, P
               I1PJJ = I1 + JJ
               I3PJJ = I3 + JJ
               TE = TE + C(I1PJJ,J)*RES(I3PJJ)
 1020       CONTINUE
 1040       IF (I.EQ.1) GO TO 1080
            DO 1060 JJ = 1, M1
               I2J = I2 + JJ
               TE = TE + C(I2J,J)*RES(I2J)
 1060       CONTINUE
 1080       GRADF(I1PJ) = TE
 1100    CONTINUE
 1120 CONTINUE
C
      IF (R.EQ.0) GO TO 1180
      I3P = N1*M + P
      DO 1160 J = 1, M
         TE = 0.D0
         DO 1140 I = 1, R
            I3PPI = I3P + I
            TE = TE + DEL(I,J)*RES(I3PPI)
 1140    CONTINUE
         GRADF(J) = GRADF(J) + TE
 1160 CONTINUE
 1180 CALL D02RAT(A,C,DEL,M,N,P,R,IR,IC,SING,ICA,AUX,MTNMAX,NMAX,MMAX2)
      IF (SING) GO TO 1220
 1200 CALL D02RAS(A,C,DEL,RES,M,N,P,R,IR,IC,UU,MTNMAX,NMAX,XAU,SING)
 1220 RETURN
C
99999 FORMAT (' SINGULAR BOUNDARY CONDITION JACOBIAN')
      END
