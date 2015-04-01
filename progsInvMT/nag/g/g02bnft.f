      SUBROUTINE G02BNF(N,M,X,IX,ITYPE,RR,IRR,KWORKA,KWORKB,WORK1,WORK2,
     *                  IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     NAG SUBROUTINE G02BNF
C     WRITTEN 27. 7.73 BY PAUL GRIFFITHS (OXFORD UNIVERSITY)
C
C     COMPUTES KENDALL AND/OR SPEARMAN RANK CORRELATION
C     COEFFICIENTS
C     FOR A SET OF DATA IN THE ARRAY X -- THE ARRAY X IS
C     OVERWRITTEN
C     BY THE ROUTINE.
C
C     USES NAG ERROR ROUTINE P01AAF
C     AND AUXILIARY ROUTINES G02BNZ
C     G02BNX
C
C
C     ABOVE DATA STATEMENT MAY BE MACHINE-DEPENDENT -- DEPENDS ON
C     NUMBER OF CHARACTERS WHICH CAN BE STORED IN A REAL VARIABLE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BNF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IRR, ITYPE, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RR(IRR,M), WORK1(M), WORK2(M), X(IX,M)
      INTEGER           KWORKA(N), KWORKB(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D, U
      INTEGER           I, IERROR, IP, IV, IVP, J, JV, M1, N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G02BNX, G02BNZ
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN, SQRT
C     .. Executable Statements ..
      IERROR = 0
      IF (ITYPE.LT.-1 .OR. ITYPE.GT.1) IERROR = 4
      IF (IX.LT.N .OR. IRR.LT.M) IERROR = 3
      IF (M.LT.2) IERROR = 2
      IF (N.LT.2) IERROR = 1
      IF (IERROR) 20, 20, 540
   20 DO 60 IV = 1, M
         DO 40 I = 1, N
            KWORKA(I) = I
            KWORKB(I) = I
   40    CONTINUE
         CALL G02BNZ(X(1,IV),N,N,KWORKA,KWORKB)
         CALL G02BNX(X(1,IV),N,N,KWORKA,KWORKB(1),WORK1(IV),WORK2(IV))
   60 CONTINUE
      N1 = N - 1
      M1 = M - 1
      IF (ITYPE) 80, 80, 220
   80 DO 200 IV = 1, M1
         IVP = IV + 1
         DO 180 JV = IVP, M
            U = 0.0D0
            DO 120 I = 1, N1
               IP = I + 1
               DO 100 J = IP, N
                  C = X(I,IV) - X(J,IV)
                  D = X(I,JV) - X(J,JV)
                  IF (C.NE.0.0D0 .AND. D.NE.0.0D0) U = U + SIGN(1.0D0,C)
     *                *SIGN(1.0D0,D)
  100          CONTINUE
  120       CONTINUE
            D = WORK1(IV)*WORK1(JV)
            IF (D) 140, 140, 160
  140       RR(JV,IV) = 0.0D0
            GO TO 180
  160       RR(JV,IV) = U/SQRT(D)
  180    CONTINUE
  200 CONTINUE
  220 IF (ITYPE) 360, 240, 240
  240 DO 340 IV = 2, M
         IVP = IV - 1
         DO 320 JV = 1, IVP
            U = 0.0D0
            DO 260 I = 1, N
               D = X(I,IV) - X(I,JV)
               U = U + D*D
  260       CONTINUE
            D = WORK2(IV)*WORK2(JV)
            IF (D) 280, 280, 300
  280       RR(JV,IV) = 0.0D0
            GO TO 320
  300       RR(JV,IV) = 0.5D0*(WORK2(IV)+WORK2(JV)-6.0D0*U)/SQRT(D)
  320    CONTINUE
  340 CONTINUE
  360 IF (ITYPE) 380, 500, 440
  380 DO 420 IV = 2, M
         IVP = IV - 1
         DO 400 JV = 1, IVP
            RR(JV,IV) = RR(IV,JV)
  400    CONTINUE
  420 CONTINUE
      GO TO 500
  440 DO 480 IV = 1, M1
         IVP = IV + 1
         DO 460 JV = IVP, M
            RR(JV,IV) = RR(IV,JV)
  460    CONTINUE
  480 CONTINUE
  500 DO 520 I = 1, M
         RR(I,I) = 1.0D0
  520 CONTINUE
      IFAIL = 0
      RETURN
  540 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
