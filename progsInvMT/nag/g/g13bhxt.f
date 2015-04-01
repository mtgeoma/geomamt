      SUBROUTINE G13BHX(STTF,NSTTF,KPST,NWD,XC,LXC,ZC,LZC,PARA,NPARA,
     *                  KPPA,XN,ZN,NNV,KUSTF)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ROUTINE G13BHX TAKES AN ARRAY XC (CONSISTING OF THE
C     VALUES OF X RELATING TO THIS INPUT SERIES IN STTF,
C     PLUS A SINGLE NEW VALUE OF X) AND A CORRESPONDING
C     ARRAY ZC, AND DERIVES UPDATED VALUES OF X AND Z.
C     AS AN OPTION IT UPDATES THE X,Z COMPONENT OF STTF
C
C     .. Scalar Arguments ..
      INTEGER           KPPA, KPST, KUSTF, LXC, LZC, NNV, NPARA, NSTTF,
     *                  NWD
C     .. Array Arguments ..
      DOUBLE PRECISION  PARA(NPARA), STTF(NSTTF), XC(LXC), XN(NNV),
     *                  ZC(LZC), ZN(NNV)
C     .. Local Scalars ..
      DOUBLE PRECISION  Q, ZERO, ZNA
      INTEGER           I, IREV, J, K, L, NGD, NGW
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C
C     NGD AND NGW HOLD THE NUMBERS OF DELTAS AND OMEGAS
C
      NGD = LZC - 1
      NGW = NWD - NGD
C
C     PROCESS A LOOP IN WHICH ONE NEW VALUE OF X IS
C     INTRODUCED EACH TIME
C
      DO 260 L = 1, NNV
C
C        TRANSFER X AND Z BLOCKS FROM STTF TO XC AND ZC
C        ONLY NEEDED ON FIRST CYCLE OF LOOP
C
         K = KPST
         IF (L.NE.1) GO TO 80
         IF (LXC.LE.1) GO TO 40
         J = LXC - 1
         DO 20 I = 1, J
            K = K + 1
            XC(I) = STTF(K)
   20    CONTINUE
   40    IF (LZC.LE.1) GO TO 80
         DO 60 I = 1, NGD
            K = K + 1
            ZC(I) = STTF(K)
   60    CONTINUE
C
C        ADD NEW VALUE OF X TO XC AND CALCULATE CORRESPONDING
C        NEW VALUE OF Z
C
   80    ZNA = ZERO
         XC(LXC) = XN(L)
         K = KPPA + 1
         IREV = K + NGW
         DO 120 I = 1, NGW
            IREV = IREV - 1
            Q = PARA(IREV)
            IF (IREV.EQ.K) GO TO 100
            Q = -Q
  100       ZNA = ZNA + Q*XC(I)
  120    CONTINUE
         IF (NGD.LE.0) GO TO 160
         IREV = K + NWD
         DO 140 I = 1, NGD
            IREV = IREV - 1
            ZNA = ZNA + PARA(IREV)*ZC(I)
  140    CONTINUE
  160    ZC(LZC) = ZNA
C
C        USE NEW VALUES OF X AND Z TO UPDATE XC AND ZC.
C        UPDATE STTF IF KUSTF = 1.
C
         J = KPST
         IF (LXC.LE.1) GO TO 200
         DO 180 I = 2, LXC
            XC(I-1) = XC(I)
            IF (KUSTF.EQ.0) GO TO 180
            J = J + 1
            STTF(J) = XC(I)
  180    CONTINUE
  200    IF (LZC.LE.1) GO TO 240
         DO 220 I = 2, LZC
            ZC(I-1) = ZC(I)
            IF (KUSTF.EQ.0) GO TO 220
            J = J + 1
            STTF(J) = ZC(I)
  220    CONTINUE
  240    ZN(L) = ZNA
  260 CONTINUE
      RETURN
      END
