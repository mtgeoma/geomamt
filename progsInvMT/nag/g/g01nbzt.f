      DOUBLE PRECISION FUNCTION G01NBZ(X,ICODE,IMU,N,IS,ISROW,ISPRTN,
     *                                 RLANDA,SMU,SLA,AA,CC,VLA,VMU,
     *                                 GAMMA,R,WORK1,WORK2,POW)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Based on routine F of Magnus and Pesaran
C
C     .. Parameters ..
      INTEGER                          ISDIM, ISPAR
      PARAMETER                        (ISDIM=12,ISPAR=77)
      INTEGER                          ISSYM
      PARAMETER                        (ISSYM=ISDIM*(ISDIM+1)/2)
      DOUBLE PRECISION                 ONE, TWO, ZERO, HALF
      PARAMETER                        (ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0,
     *                                 HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          ICODE, IMU, IS, ISROW, N
C     .. Array Arguments ..
      DOUBLE PRECISION                 AA(*), CC(*), GAMMA(*), POW(N,*),
     *                                 R(*), RLANDA(*), SLA(N), SMU(N),
     *                                 VLA(N), VMU(N), WORK1(*),
     *                                 WORK2(*)
      INTEGER                          ISPRTN(ISPAR,ISDIM)
C     .. Local Scalars ..
      DOUBLE PRECISION                 F2, F3, F4, FACT, RB, SP, SPI,
     *                                 SPJ, SUM, SUM1, SUM2, TRACE
      INTEGER                          I, II, IJ, J, K, L, LCODE, NN
C     .. Local Arrays ..
      DOUBLE PRECISION                 PP(ISSYM), TR1(ISDIM), TR2(ISDIM)
      INTEGER                          NK(ISDIM)
C     .. External Functions ..
      DOUBLE PRECISION                 DDOT, G01NBX, G02AAR
      EXTERNAL                         DDOT, G01NBX, G02AAR
C     .. External Subroutines ..
      EXTERNAL                         DCOPY, DGEMV, DSPMV, G02AAS
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE, EXP, SQRT
C     .. Executable Statements ..
      NN = N*(N+1)/2
      F2 = ONE
      TRACE = ZERO
      IJ = 0
      DO 40 I = 1, N
         SPI = ONE + TWO*X*RLANDA(I)
         SPI = ONE/SQRT(SPI)
         F2 = F2*SPI
         IF (I.LT.IS) F2 = F2*X
         IF (IMU.NE.0) VMU(I) = SMU(I)*SPI
         IF (ICODE.EQ.2) VLA(I) = SLA(I)*SPI
         DO 20 J = 1, I
            SPJ = ONE + TWO*X*RLANDA(J)
            SPJ = ONE/SQRT(SPJ)
            SP = SPI*SPJ
            IJ = IJ + 1
            R(IJ) = AA(IJ)*SP
            IF (ICODE.EQ.3) THEN
               GAMMA(IJ) = CC(IJ)*SP
               IF (I.EQ.J) TRACE = TRACE + GAMMA(IJ)
            END IF
   20    CONTINUE
   40 CONTINUE
      IF (IS.GT.N+1) F2 = F2*(X**(IS-N-1))
      F3 = ONE
      IF (IMU.NE.0) THEN
         F3 = ZERO
         DO 60 J = 1, N
            SP = SMU(J)*SMU(J) - VMU(J)*VMU(J)
            F3 = F3 + SP
   60    CONTINUE
         F3 = EXP(-HALF*F3)
      END IF
      DO 120 L = 1, IS
         IF (L.EQ.1) THEN
            CALL DCOPY(NN,R,1,WORK1,1)
         ELSE
            CALL G02AAS('U',N,R,WORK1,WORK2)
            CALL DCOPY(NN,WORK2,1,WORK1,1)
         END IF
         TR1(L) = G02AAR('U',N,WORK1)
         IF (IMU.NE.0) CALL DSPMV('U',N,1.0D0,WORK1,VMU,1,0.0D0,POW(1,L)
     *                            ,1)
         IF (ICODE.EQ.3) THEN
            SUM = ZERO
            IJ = 0
            DO 100 I = 1, N
               DO 80 J = 1, I - 1
                  IJ = IJ + 1
                  SUM = SUM + TWO*WORK1(IJ)*GAMMA(IJ)
   80          CONTINUE
               IJ = IJ + 1
               SUM = SUM + WORK1(IJ)*GAMMA(IJ)
  100       CONTINUE
            TR2(L) = SUM
         END IF
  120 CONTINUE
      IF (IMU.NE.0) THEN
         IF (ICODE.EQ.2) THEN
            SUM = DDOT(N,VLA,1,VMU,1)
            TRACE = TRACE + SUM
         ELSE IF (ICODE.EQ.3) THEN
            SUM = ZERO
            IJ = 0
            DO 160 I = 1, N
               DO 140 J = 1, I - 1
                  IJ = IJ + 1
                  SUM = SUM + TWO*VMU(I)*GAMMA(IJ)*VMU(J)
  140          CONTINUE
               IJ = IJ + 1
               SUM = SUM + VMU(I)*GAMMA(IJ)*VMU(I)
  160       CONTINUE
            TRACE = TRACE + SUM
         END IF
         DO 180 L = 1, IS
            SUM1 = DDOT(N,VMU,1,POW(1,L),1)
            IF (ICODE.EQ.3) CALL DSPMV('U',N,1.0D0,GAMMA,VMU,1,0.0D0,
     *                                 VLA,1)
            IF (ICODE.GE.2) SUM2 = DDOT(N,VLA,1,POW(1,L),1)
            TR1(L) = TR1(L) + DBLE(L)*SUM1
            IF (ICODE.EQ.2) THEN
               TR2(L) = SUM2
            ELSE IF (ICODE.EQ.3) THEN
               TR2(L) = TR2(L) + TWO*SUM2
            END IF
  180    CONTINUE
         IF (ICODE.EQ.3) THEN
            IJ = 1
            DO 200 I = 1, IS
               CALL DSPMV('U',N,1.0D0,GAMMA,POW(1,I),1,0.0D0,VLA,1)
               CALL DGEMV('T',N,I,1.0D0,POW,N,VLA,1,0.0D0,PP(IJ),1)
               IJ = IJ + I
  200       CONTINUE
         END IF
      END IF
      F4 = ZERO
      LCODE = ICODE
      IF (ICODE.EQ.2 .AND. IMU.EQ.0) LCODE = 0
      IF (ICODE.EQ.3 .AND. IMU.EQ.0) LCODE = 2
      DO 260 I = 1, ISROW
         RB = 1
         FACT = 1
         DO 240 L = 1, IS
            FACT = FACT*2*L
            K = ISPRTN(I,L)
            NK(L) = K
            DO 220 II = 1, K
               RB = RB*II*2*L
  220       CONTINUE
  240    CONTINUE
         RB = FACT/RB
         F4 = F4 + RB*G01NBX(LCODE,IS,NK,TR1,TR2,PP,TRACE)
  260 CONTINUE
      G01NBZ = F2*F3*F4
      RETURN
      END
