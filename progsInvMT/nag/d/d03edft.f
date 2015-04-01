      SUBROUTINE D03EDF(NGX,NGY,LDA,A,RHS,UB,MAXIT,ACC,US,U,IOUT,NUMIT,
     *                  IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-534 (FEB 1987).
C     MARK 13 REVISED. IER-624 (APR 1988).
C     MARK 14 REVISED. IER-811 (DEC 1989).
C
C     MGD1 author -- P Wesseling, Delft University, Holland
C     Modified by -- C P Thompson and G J McCarthy, CSSD, AERE Harwell
C     Modified for NAG usage -- NAG Central Office
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03EDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC
      INTEGER           IFAIL, IOUT, LDA, MAXIT, NGX, NGY, NUMIT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,7), RHS(LDA), U(LDA), UB(NGX*NGY), US(LDA)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACCX, EPS, RESN2, RESN3, RESN32
      INTEGER           I, IERROR, J, JOUT, K, KB, KBET, KBX, KBY, KDEL,
     *                  KDELX2, KDELY2, KEND, KEPSZ, KEPX1, KEPX2,
     *                  KEPY1, KEPY2, KGAM, KGAMX, KGAMY, KK, LEV,
     *                  LEVELS, M, MUB, MUBG, MUG, NADV, NG, NG1, NG2,
     *                  NG3, NG4, NGXY, NIT, NKBY, NPC, NPCC, NPF, NPFF,
     *                  NPFO, NREC, NXC, NXF, NYC, NYF
C     .. Local Arrays ..
      INTEGER           IATX(7), IATY(7), NGP(12), NGRIDX(12),
     *                  NGRIDY(12)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DNRM2, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D03EDT, D03EDU, D03EDV, D03EDW, D03EDX, D03EDY,
     *                  D03EDZ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              IATX/0, 1, -1, 0, 1, -1, 0/
      DATA              IATY/-1, -1, 0, 0, 0, 1, 1/
C     .. Executable Statements ..
C
      CALL X04ABF(0,NADV)
C
      IERROR = 0
      JOUT = IOUT
      NIT = 0
      IF (NGX.LT.3) THEN
         WRITE (REC,FMT=99999) NGX
         GO TO 540
      ELSE IF (NGY.LT.3) THEN
         WRITE (REC,FMT=99998) NGY
         GO TO 540
      ELSE IF (ACC.LT.0.0D0) THEN
         WRITE (REC,FMT=99996) ACC
         GO TO 540
      ELSE IF (MAXIT.LT.0) THEN
         WRITE (REC,FMT=99995) MAXIT
         GO TO 540
      ELSE IF (IOUT.LT.0 .OR. IOUT.GT.8) THEN
         WRITE (REC,FMT=99997) IOUT
         GO TO 540
      END IF
      EPS = X02AJF()
      IF (ACC.EQ.0.0D0) THEN
         ACCX = EPS
      ELSE IF (ACC.LT.EPS) THEN
         IERROR = 4
         ACCX = EPS
      ELSE
         ACCX = ACC
      END IF
      LEVELS = 0
      NG3 = NGX
      NG4 = NGY
      DO 20 I = 1, 12
         LEVELS = LEVELS + 1
         NG1 = (NG3/2)*2
         NG2 = (NG4/2)*2
         IF (NG3.EQ.NG1 .OR. NG4.EQ.NG2) GO TO 40
         IF (NG3.LE.3 .OR. NG4.LE.3) GO TO 40
         NG3 = (NG3+1)/2
         NG4 = (NG4+1)/2
   20 CONTINUE
   40 NGRIDX(LEVELS) = NGX
      NGRIDY(LEVELS) = NGY
      DO 60 LEV = LEVELS - 1, 1, -1
         NGRIDX(LEV) = (NGRIDX(LEV+1)+1)/2
         NGRIDY(LEV) = (NGRIDY(LEV+1)+1)/2
   60 CONTINUE
      NGP(LEVELS) = NGRIDX(LEVELS)*NGRIDY(LEVELS)
      DO 80 LEV = LEVELS - 1, 1, -1
         NGP(LEV) = NGP(LEV+1) + NGRIDX(LEV)*NGRIDY(LEV)
   80 CONTINUE
      NG = NGP(1)
      IF (LDA.LT.NG) THEN
         WRITE (REC,FMT=99994) NG, LDA
         GO TO 540
      END IF
C
      IF (JOUT.GE.5) THEN
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,' The Finite Difference Operator')
         CALL D03EDU(A,LEVELS,LDA,NGP,NGRIDX,NGRIDY)
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,' The Right Hand Side')
         CALL D03EDT(RHS,LEVELS,NGP,NGRIDX,NGRIDY)
      END IF
C
C     Construction of coarse grid matrices with
C     Galerkin approximation.
C
      DO 120 J = 1, 7
         DO 100 K = NGP(LEVELS) + 1, NGP(1)
            A(K,J) = 0.0D0
  100    CONTINUE
  120 CONTINUE
      DO 280 LEV = LEVELS - 1, 1, -1
         NXF = NGRIDX(LEV+1)
         NYF = NGRIDY(LEV+1)
         NXC = NGRIDX(LEV)
         NYC = NGRIDY(LEV)
         NPCC = NGP(LEV) - NXC*NYC - NXC
         NPFF = NGP(LEV+1) - NXF*NYF - 2*NXF - 1
         DO 220 KDEL = 1, 7
            KDELX2 = 2*IATX(KDEL)
            KDELY2 = 2*IATY(KDEL)
            DO 200 KBET = 1, 7
               KBX = IATX(KBET)
               KBY = IATY(KBET)
               MUB = 1
               IF (KBET.EQ.4) MUB = 2
               NKBY = KBY*NXF
               KEPX1 = KDELX2 - KBX
               KEPY1 = KDELY2 - KBY
               DO 180 KGAM = 1, 7
                  KGAMX = IATX(KGAM)
                  KGAMY = IATY(KGAM)
                  MUG = 1
                  IF (KGAM.EQ.4) MUG = 2
                  KEPX2 = KEPX1 + KGAMX
                  IF (ABS(KEPX2).GT.1) GO TO 180
                  KEPY2 = KEPY1 + KGAMY
                  IF (ABS(KEPY2).GT.1) GO TO 180
                  IF (ABS(KEPX2+KEPY2).GT.1) GO TO 180
                  KEPSZ = KEPX2 + 3*KEPY2 + 4
                  MUBG = MUB*MUG
                  NPFO = NPFF + NKBY + KBX
                  DO 160 K = MAX(1,1-KBY), MIN(NYC,NYC-KBY)
                     NPC = NPCC + NXC*K
                     NPF = NPFO + 2*NXF*K
                     DO 140 J = MAX(1,1-KBX), MIN(NXC,NXC-KBX)
                        A(NPC+J,KDEL) = A(NPC+J,KDEL) + MUBG*A(NPF+2*J,
     *                                  KEPSZ)
  140                CONTINUE
  160             CONTINUE
  180          CONTINUE
  200       CONTINUE
  220    CONTINUE
         KEND = NGP(LEV)
         KB = KEND - NXC*NYC + 1
         DO 260 J = 1, 7
            DO 240 K = KB, KEND
               A(K,J) = A(K,J)/4.0D0
  240       CONTINUE
  260    CONTINUE
  280 CONTINUE
C
C     End of computation of coarse grid matrices
C
      IF (JOUT.GE.6 .AND. LEVELS.GT.1) THEN
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,' Galerkin Coarse Grid Approximations')
         DO 300 LEV = LEVELS - 1, 1, -1
            CALL D03EDU(A,LEV,LDA,NGP,NGRIDX,NGRIDY)
  300    CONTINUE
      END IF
C
C     Computation of Incomplete Crout-Decompositions
C
      DO 380 LEV = 1, LEVELS
         M = NGRIDX(LEV)
         KK = NGP(LEV) - M*NGRIDY(LEV)
         DO 320 K = KK + 2, KK + M
            A(K,3) = A(K,3)/A(K-1,4)
            A(K,4) = A(K,4) - A(K,3)*A(K-1,5)
            A(K,6) = A(K,6) - A(K,3)*A(K-1,7)
  320    CONTINUE
         K = KK + M
         A(K,2) = A(K,2)/A(K-M+1,4)
         A(K,4) = A(K,4) - A(K,2)*A(K-M+1,6)
         A(K,5) = A(K,5) - A(K,2)*A(K-M+1,7)
         DO 340 K = KK + M + 1, NGP(LEV) - M + 1
            A(K,1) = A(K,1)/A(K-M,4)
            A(K,2) = (A(K,2)-A(K,1)*A(K-M,5))/A(K-M+1,4)
            A(K,3) = (A(K,3)-A(K,1)*A(K-M,6))/A(K-1,4)
            A(K,4) = A(K,4) - A(K,3)*A(K-1,5) - A(K,2)*A(K-M+1,6) - A(K,
     *               1)*A(K-M,7)
            A(K,5) = A(K,5) - A(K,2)*A(K-M+1,7)
            A(K,6) = A(K,6) - A(K,3)*A(K-1,7)
  340    CONTINUE
         DO 360 K = NGP(LEV) - M + 2, NGP(LEV)
            A(K,1) = A(K,1)/A(K-M,4)
            A(K,2) = (A(K,2)-A(K,1)*A(K-M,5))/A(K-M+1,4)
            A(K,3) = (A(K,3)-A(K,1)*A(K-M,6))/A(K-1,4)
            A(K,4) = A(K,4) - A(K,3)*A(K-1,5) - A(K,2)*A(K-M+1,6) - A(K,
     *               1)*A(K-M,7)
            A(K,5) = A(K,5) - A(K,2)*A(K-M+1,7)
  360    CONTINUE
         A(NGP(LEV),5) = 0.0D0
  380 CONTINUE
C
      IF (JOUT.GE.7) THEN
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,' Incomplete Crout Decompositions')
         DO 400 LEV = 1, LEVELS
            CALL D03EDU(A,LEV,LDA,NGP,NGRIDX,NGRIDY)
  400    CONTINUE
      END IF
C
C     Find initial U and UB
C
      NGXY = NGP(LEVELS)
      CALL D03EDY(A,US,UB,RHS,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
      CALL D03EDV(A,U,US,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
      NIT = 0
C
C     Compute initial residual and test
C
      CALL D03EDZ(A,US,U,UB,NGRIDX(LEVELS),NGXY,LDA)
      RESN3 = DNRM2(NGXY,UB,1)
      IF (JOUT.GE.2) THEN
         WRITE (REC,FMT=99992) NIT, RESN3
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,REC(1))
      END IF
      IF (RESN3.LE.ACCX .AND. MAXIT.EQ.0) GO TO 520
C
C     Start of Multigrid iterations
C
      DO 500 NIT = 1, MAXIT
         IF (JOUT.GE.8) THEN
            WRITE (REC,FMT=99993) NIT - 1
            CALL X04BAF(NADV,' ')
            CALL X04BAF(NADV,REC(1))
            CALL D03EDT(UB,LEVELS,NGP,NGRIDX,NGRIDY)
         END IF
         RESN2 = RESN3 + ACCX*1.0D-06
         IF (LEVELS.GT.1) THEN
            CALL D03EDW(RHS,UB,LEVELS-1,NGP(LEVELS-1),NGXY,NGP,NGRIDX,
     *                  NGRIDY)
            DO 420 LEV = LEVELS - 2, 1, -1
               CALL D03EDW(RHS,RHS,LEV,NGP(LEV),NGP(LEV+1),NGP,NGRIDX,
     *                     NGRIDY)
  420       CONTINUE
            CALL D03EDV(A,U,RHS,1,LDA,NG,NGP,NGRIDX,NGRIDY)
            DO 440 LEV = 2, LEVELS - 1
               CALL D03EDX(U,U,LEV,NGP(LEV),NGP(LEV-1),NGP,NGRIDX,
     *                     NGRIDY)
               CALL D03EDY(A,US,U,RHS,LEV,LDA,NGP(LEV),NGP,NGRIDX,
     *                     NGRIDY)
               CALL D03EDV(A,U,US,LEV,LDA,NGP(LEV),NGP,NGRIDX,NGRIDY)
  440       CONTINUE
            CALL D03EDX(UB,U,LEVELS,NGXY,NGP(LEVELS-1),NGP,NGRIDX,
     *                  NGRIDY)
            DO 460 K = 1, NGP(LEVELS)
               UB(K) = UB(K) + U(K)
  460       CONTINUE
            CALL D03EDY(A,U,UB,RHS,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
            CALL D03EDV(A,U,U,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
         ELSE IF (LEVELS.EQ.1) THEN
            CALL D03EDV(A,US,UB,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
            DO 480 K = 1, NGP(LEVELS)
               UB(K) = US(K) + U(K)
  480       CONTINUE
            CALL D03EDY(A,U,UB,RHS,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
            CALL D03EDV(A,U,U,LEVELS,LDA,NGXY,NGP,NGRIDX,NGRIDY)
         END IF
         CALL D03EDZ(A,US,U,UB,NGRIDX(LEVELS),NGXY,LDA)
         RESN3 = DNRM2(NGXY,UB,1)
         RESN32 = RESN3/RESN2
         IF (RESN32.GE.1.0D0 .AND. NIT.GT.1) THEN
            IERROR = 3
            IF (JOUT.GT.0) THEN
               JOUT = MAX(JOUT,2)
               WRITE (REC,FMT=99988)
               CALL X04BAF(NADV,' ')
               CALL X04BAF(NADV,REC(1))
            END IF
         END IF
         IF (JOUT.GE.2) THEN
            WRITE (REC,FMT=99991) NIT, RESN3, RESN32
            CALL X04BAF(NADV,' ')
            CALL X04BAF(NADV,REC(1))
         END IF
C        Test for convergence
         IF (RESN3.LE.ACCX) THEN
            IF (IERROR.EQ.3) IERROR = 0
            GO TO 520
         END IF
  500 CONTINUE
      NIT = NIT - 1
      IF (IERROR.NE.3) IERROR = 2
C
C     Termination
C
  520 US(1) = RESN3
      IF (JOUT.GE.4) THEN
         WRITE (REC,FMT=99993) NIT
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,REC(1))
         CALL D03EDT(UB,LEVELS,NGP,NGRIDX,NGRIDY)
      END IF
      IF (JOUT.EQ.1 .OR. JOUT.GE.3) THEN
         CALL X04BAF(NADV,' ')
         CALL X04BAF(NADV,' Numerical Solution')
         CALL D03EDT(U,LEVELS,NGP,NGRIDX,NGRIDY)
      END IF
C
      IF (IERROR.EQ.2 .OR. IERROR.EQ.3) THEN
         WRITE (REC,FMT=99990) MAXIT, RESN3, ACCX
         IF (IERROR.EQ.2) WRITE (REC(3),FMT=99987)
         IF (IERROR.EQ.3) WRITE (REC(3),FMT=99986)
         NREC = 3
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC,FMT=99989) RESN3, EPS, ACC
         NREC = 3
      END IF
      GO TO 560
  540 IERROR = 1
      NREC = 1
  560 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      NUMIT = NIT
      RETURN
C
99999 FORMAT (' ** On entry, NGX is less than 3: NGX = ',I16)
99998 FORMAT (' ** On entry, NGY is less than 3: NGY = ',I16)
99997 FORMAT (' ** On entry, IOUT is not in the range 0 to 8: IOUT = ',
     *       I16)
99996 FORMAT (' ** On entry, ACC is less than 0.0: ACC =',1P,D13.5)
99995 FORMAT (' ** On entry, MAXIT is less than 0: MAXIT =',I16)
99994 FORMAT (' ** On entry, LDA must be at least',I8,': LDA =',I16)
99993 FORMAT (' Residual For Iteration',I4)
99992 FORMAT (' Iteration',I4,'  Residual norm =',1P,D10.2)
99991 FORMAT (' Iteration',I4,'  Residual norm =',1P,D10.2,'  Reductio',
     *       'n factor =',D10.2)
99990 FORMAT (' ** After MAXIT iterations the residual norm is not les',
     *       's than the tolerance',/'    MAXIT =',I8,'  residual norm',
     *       ' =',1P,D10.2,'  tolerance =',D10.2)
99989 FORMAT (' ** On entry, ACC is less than machine precision. The r',
     *       'outine terminated',/'    because the residual norm is le',
     *       'ss th','an machine precision.',/'    residual norm =',1P,
     *       D10.2,'  mach','ine precision =',D10.2,'  ACC =',D10.2)
99988 FORMAT (' ** WARNING: the residual norm has increased **')
99987 FORMAT ('    The residual norm has decreased at each iteration a',
     *       'fter the first')
99986 FORMAT ('    The residual norm increased at one or more iteratio',
     *       'ns after the first')
      END
