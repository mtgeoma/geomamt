      SUBROUTINE D03FAF(XS,XF,L,LBDCND,BDXS,BDXF,YS,YF,M,MBDCND,BDYS,
     *                  BDYF,ZS,ZF,N,NBDCND,BDZS,BDZF,LAMBDA,LDIMF,
     *                  MDIMF,F,PERTRB,W,LWRK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03FAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  LAMBDA, PERTRB, XF, XS, YF, YS, ZF, ZS
      INTEGER           IFAIL, L, LBDCND, LDIMF, LWRK, M, MBDCND, MDIMF,
     *                  N, NBDCND
C     .. Array Arguments ..
      DOUBLE PRECISION  BDXF(MDIMF,N+1), BDXS(MDIMF,N+1),
     *                  BDYF(LDIMF,N+1), BDYS(LDIMF,N+1),
     *                  BDZF(LDIMF,M+1), BDZS(LDIMF,M+1),
     *                  F(LDIMF,MDIMF,N+1), W(LWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, C3, DX, DY, DZ, S, S1, S2, TWBYDX,
     *                  TWBYDY, TWBYDZ, XLP, YLP, ZLP
      INTEGER           I, IERROR, IR, IWB, IWC, IWW, J, JTRIGL, JTRIGM,
     *                  JWORK1, JWORK2, K, LSTART, LSTOP, LUNK, LWKACT,
     *                  MSTART, MSTOP, MUNK, NPEROD, NREC, NSTART,
     *                  NSTOP, NUNK
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D03FAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE
C     .. Executable Statements ..
C
C     Check for invalid input
C
      IERROR = 0
      NREC = 0
      IF (XF.LE.XS) THEN
         WRITE (REC,FMT=99999) XS, XF
         IERROR = 1
         NREC = 2
      ELSE IF (L.LT.5) THEN
         WRITE (REC,FMT=99998) L
         IERROR = 1
         NREC = 1
      ELSE IF (LBDCND.LT.0 .OR. LBDCND.GT.4) THEN
         WRITE (REC,FMT=99997) LBDCND
         NREC = 2
         IERROR = 1
      ELSE IF (YF.LE.YS) THEN
         WRITE (REC,FMT=99996) YS, YF
         NREC = 2
         IERROR = 1
      ELSE IF (M.LT.5) THEN
         WRITE (REC,FMT=99995) M
         NREC = 1
         IERROR = 1
      ELSE IF (MBDCND.LT.0 .OR. MBDCND.GT.4) THEN
         WRITE (REC,FMT=99994) MBDCND
         NREC = 2
         IERROR = 1
      ELSE IF (ZF.LE.ZS) THEN
         WRITE (REC,FMT=99993) ZS, ZF
         NREC = 2
         IERROR = 1
      ELSE IF (N.LT.5) THEN
         WRITE (REC,FMT=99992) N
         NREC = 1
         IERROR = 1
      ELSE IF (NBDCND.LT.0 .OR. NBDCND.GT.4) THEN
         WRITE (REC,FMT=99991) NBDCND
         NREC = 2
         IERROR = 1
      ELSE IF (LDIMF.LT.L+1) THEN
         WRITE (REC,FMT=99990) LDIMF, L
         NREC = 2
         IERROR = 1
      ELSE IF (MDIMF.LT.M+1) THEN
         WRITE (REC,FMT=99989) MDIMF, M
         NREC = 2
         IERROR = 1
      END IF
      IF (IERROR.EQ.1) THEN
         IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
         RETURN
      END IF
      IF (LAMBDA.GT.0.0D0) THEN
C
C        A solution may not exist, but carry on anyway
C
         IERROR = 3
         WRITE (REC,FMT=99987)
         NREC = 1
      END IF
C
      DY = (YF-YS)/DBLE(M)
      TWBYDY = 2.D0/DY
      C2 = 1.D0/(DY**2)
      MSTART = 1
      MSTOP = M
      IF (MBDCND.EQ.1 .OR. MBDCND.EQ.2) MSTART = 2
      IF (MBDCND.EQ.2 .OR. MBDCND.EQ.3) MSTOP = M + 1
      MUNK = MSTOP - MSTART + 1
C
      DZ = (ZF-ZS)/DBLE(N)
      TWBYDZ = 2.D0/DZ
      C3 = 1.D0/(DZ**2)
      NSTART = 1
      NSTOP = N
      IF (NBDCND.EQ.1 .OR. NBDCND.EQ.2) NSTART = 2
      IF (NBDCND.EQ.2 .OR. NBDCND.EQ.3) NSTOP = N + 1
      NUNK = NSTOP - NSTART + 1
C
      DX = (XF-XS)/DBLE(L)
      C1 = 1.D0/(DX**2)
      TWBYDX = 2.D0/DX
      LSTART = 1
      LSTOP = L
      IF (LBDCND.EQ.1 .OR. LBDCND.EQ.2) LSTART = 2
      IF (LBDCND.EQ.2 .OR. LBDCND.EQ.3) LSTOP = L + 1
      LUNK = LSTOP - LSTART + 1
C
C     Calculate the exact workspace requirements
C
C     First the storage for A, B, C, XRT and YRT
C
      LWKACT = 3*NUNK + LUNK + MUNK
C
C     Now the workspace for WORK1, WORK2 and TRIGL
C
      IF (LBDCND.EQ.0) THEN
         JWORK1 = NUNK*LUNK
         JWORK2 = NUNK*LUNK
         JTRIGL = 2*LUNK
      ELSE IF (LBDCND.EQ.1) THEN
         JWORK1 = NUNK*(LUNK+1)
         JWORK2 = NUNK*(LUNK+1)
         JTRIGL = 2*(LUNK+1)
      ELSE IF (LBDCND.EQ.2) THEN
         JWORK1 = NUNK*LUNK
         JWORK2 = NUNK*LUNK
         JTRIGL = 2*LUNK
      ELSE IF (LBDCND.EQ.3) THEN
         JWORK1 = NUNK*LUNK
         JWORK2 = NUNK*(LUNK-1)
         JTRIGL = 2*(LUNK-1)
      ELSE IF (LBDCND.EQ.4) THEN
         JWORK1 = NUNK*LUNK
         JWORK2 = NUNK*LUNK
         JTRIGL = 2*LUNK
      END IF
C
C     WORK1 and WORK2 need to be big enough for transformations in both
C     the x and y directions. TRIGM has to be kept separately
C     from TRIGL
C
      IF (MBDCND.EQ.0) THEN
         JWORK1 = MAX(JWORK1,NUNK*MUNK)
         JWORK2 = MAX(JWORK2,NUNK*MUNK)
         JTRIGM = 2*MUNK
      ELSE IF (MBDCND.EQ.1) THEN
         JWORK1 = MAX(JWORK1,NUNK*(MUNK+1))
         JWORK2 = MAX(JWORK2,NUNK*(MUNK+1))
         JTRIGM = 2*(MUNK+1)
      ELSE IF (MBDCND.EQ.2) THEN
         JWORK1 = MAX(JWORK1,NUNK*MUNK)
         JWORK2 = MAX(JWORK2,NUNK*MUNK)
         JTRIGM = 2*MUNK
      ELSE IF (MBDCND.EQ.3) THEN
         JWORK1 = MAX(JWORK1,NUNK*MUNK)
         JWORK2 = MAX(JWORK2,NUNK*(MUNK-1))
         JTRIGM = 2*(MUNK-1)
      ELSE IF (MBDCND.EQ.4) THEN
         JWORK1 = MAX(JWORK1,NUNK*MUNK)
         JWORK2 = MAX(JWORK2,NUNK*MUNK)
         JTRIGM = 2*MUNK
      END IF
      LWKACT = LWKACT + JWORK1 + JWORK2 + JTRIGL + JTRIGM
C
C     Exit if not enough workspace
C
      IF (LWRK.LT.LWKACT) THEN
         WRITE (REC,FMT=99988) LWKACT, LWRK
         NREC = 2
         IERROR = 2
         IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
         RETURN
      END IF
C
C     Enter boundary data for X-boundaries.
C
      IF (LBDCND.EQ.1 .OR. LBDCND.EQ.2) THEN
C
C        Solution specified at x = XS
C
         DO 40 K = NSTART, NSTOP
            DO 20 J = MSTART, MSTOP
               F(2,J,K) = F(2,J,K) - C1*F(1,J,K)
   20       CONTINUE
   40    CONTINUE
      END IF
C
      IF (LBDCND.EQ.3 .OR. LBDCND.EQ.4) THEN
C
C        Derivative specified at x = XS
C
         DO 80 K = NSTART, NSTOP
            DO 60 J = MSTART, MSTOP
               F(1,J,K) = F(1,J,K) + TWBYDX*BDXS(J,K)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      IF (LBDCND.EQ.1 .OR. LBDCND.EQ.4) THEN
C
C        Solution specified at x = XF
C
         DO 120 K = NSTART, NSTOP
            DO 100 J = MSTART, MSTOP
               F(L,J,K) = F(L,J,K) - C1*F(L+1,J,K)
  100       CONTINUE
  120    CONTINUE
      END IF
C
      IF (LBDCND.EQ.2 .OR. LBDCND.EQ.3) THEN
C
C        Derivative specified at x = XF
C
         DO 160 K = NSTART, NSTOP
            DO 140 J = MSTART, MSTOP
               F(L+1,J,K) = F(L+1,J,K) - TWBYDX*BDXF(J,K)
  140       CONTINUE
  160    CONTINUE
      END IF
C
C     Enter boundary data for Y-boundaries.
C
      IF (MBDCND.EQ.1 .OR. MBDCND.EQ.2) THEN
C
C        Solution specified at y = YS
C
         DO 200 K = NSTART, NSTOP
            DO 180 I = LSTART, LSTOP
               F(I,2,K) = F(I,2,K) - C2*F(I,1,K)
  180       CONTINUE
  200    CONTINUE
      END IF
C
      IF (MBDCND.EQ.3 .OR. MBDCND.EQ.4) THEN
C
C        Derivative specified at y = YS
C
         DO 240 K = NSTART, NSTOP
            DO 220 I = LSTART, LSTOP
               F(I,1,K) = F(I,1,K) + TWBYDY*BDYS(I,K)
  220       CONTINUE
  240    CONTINUE
      END IF
C
      IF (MBDCND.EQ.1 .OR. MBDCND.EQ.4) THEN
C
C        Solution specified at y = YF
C
         DO 280 K = NSTART, NSTOP
            DO 260 I = LSTART, LSTOP
               F(I,M,K) = F(I,M,K) - C2*F(I,M+1,K)
  260       CONTINUE
  280    CONTINUE
      END IF
C
      IF (MBDCND.EQ.2 .OR. MBDCND.EQ.3) THEN
C
C        Derivative specified at y = YF
C
         DO 320 K = NSTART, NSTOP
            DO 300 I = LSTART, LSTOP
               F(I,M+1,K) = F(I,M+1,K) - TWBYDY*BDYF(I,K)
  300       CONTINUE
  320    CONTINUE
      END IF
C
C     Enter boundary data for Z-boundaries.
C
      IF (NBDCND.EQ.1 .OR. NBDCND.EQ.2) THEN
C
C        Solution specified at z = ZS
C
         DO 360 J = MSTART, MSTOP
            DO 340 I = LSTART, LSTOP
               F(I,J,2) = F(I,J,2) - C3*F(I,J,1)
  340       CONTINUE
  360    CONTINUE
      END IF
C
      IF (NBDCND.EQ.3 .OR. NBDCND.EQ.4) THEN
C
C        Derivative specified at z = ZS
C
         DO 400 J = MSTART, MSTOP
            DO 380 I = LSTART, LSTOP
               F(I,J,1) = F(I,J,1) + TWBYDZ*BDZS(I,J)
  380       CONTINUE
  400    CONTINUE
      END IF
C
      IF (NBDCND.EQ.1 .OR. NBDCND.EQ.4) THEN
C
C        Solution specified at z = ZF
C
         DO 440 J = MSTART, MSTOP
            DO 420 I = LSTART, LSTOP
               F(I,J,N) = F(I,J,N) - C3*F(I,J,N+1)
  420       CONTINUE
  440    CONTINUE
      END IF
C
      IF (NBDCND.EQ.2 .OR. NBDCND.EQ.3) THEN
C
C        Derivative specified at z = ZF
C
         DO 480 J = MSTART, MSTOP
            DO 460 I = LSTART, LSTOP
               F(I,J,N+1) = F(I,J,N+1) - TWBYDZ*BDZF(I,J)
  460       CONTINUE
  480    CONTINUE
      END IF
C
C     Define A, B, C coefficients in W-array.
C
      IWB = NUNK + 1
      IWC = IWB + NUNK
      IWW = IWC + NUNK
      DO 500 K = 1, NUNK
         I = IWC + K - 1
         W(K) = C3
         W(I) = C3
         I = IWB + K - 1
         W(I) = -2.0D0*C3 + LAMBDA
  500 CONTINUE
      IF (NBDCND.EQ.2 .OR. NBDCND.EQ.3) W(IWB-1) = 2.0D0*C3
      IF (NBDCND.EQ.3 .OR. NBDCND.EQ.4) W(IWC) = 2.0D0*C3
      PERTRB = 0.0D0
C
C     For singular problems adjust data to ensure a
C     solution will exist.
C
      IF ((LBDCND.EQ.0 .OR. LBDCND.EQ.3)
     *    .AND. (MBDCND.EQ.0 .OR. MBDCND.EQ.3)
     *    .AND. (NBDCND.EQ.0 .OR. NBDCND.EQ.3) .AND. LAMBDA.EQ.0.0D0)
     *    THEN
C
C           Boundary conditions are either completely Neumann or
C           completely Periodic.
C
C
C           Poisson's equation
C
         XLP = (3+LBDCND)/3
         YLP = (3+MBDCND)/3
         ZLP = (3+NBDCND)/3
         S1 = 0.0D0
         DO 580 K = 2, NSTOP - 1
            DO 540 J = 2, MSTOP - 1
               DO 520 I = 2, LSTOP - 1
                  S1 = S1 + F(I,J,K)
  520          CONTINUE
               S1 = S1 + (F(1,J,K)+F(LSTOP,J,K))/XLP
  540       CONTINUE
            S2 = 0.0D0
            DO 560 I = 2, LSTOP - 1
               S2 = S2 + F(I,1,K) + F(I,MSTOP,K)
  560       CONTINUE
            S2 = (S2+(F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)
     *           +F(LSTOP,MSTOP,K))/XLP)/YLP
            S1 = S1 + S2
  580    CONTINUE
         S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)
     *       +F(1,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)
     *       +F(LSTOP,MSTOP,NSTOP))/(XLP*YLP)
         DO 620 J = 2, MSTOP - 1
            DO 600 I = 2, LSTOP - 1
               S = S + F(I,J,1) + F(I,J,NSTOP)
  600       CONTINUE
  620    CONTINUE
         S2 = 0.0D0
         DO 640 I = 2, LSTOP - 1
            S2 = S2 + F(I,1,1) + F(I,1,NSTOP) + F(I,MSTOP,1) + F(I,
     *           MSTOP,NSTOP)
  640    CONTINUE
         S = S2/YLP + S
         S2 = 0.0D0
         DO 660 J = 2, MSTOP - 1
            S2 = S2 + F(1,J,1) + F(1,J,NSTOP) + F(LSTOP,J,1) + F(LSTOP,
     *           J,NSTOP)
  660    CONTINUE
         S = S2/XLP + S
         PERTRB = (S/ZLP+S1)/((DBLE(LUNK+1)-XLP)*(DBLE(MUNK+1)-YLP)
     *            *(DBLE(NUNK+1)-ZLP))
         DO 720 K = 1, NUNK
            DO 700 J = 1, MUNK
               DO 680 I = 1, LUNK
                  F(I,J,K) = F(I,J,K) - PERTRB
  680          CONTINUE
  700       CONTINUE
  720    CONTINUE
      END IF
      NPEROD = 0
      IF (NBDCND.NE.0) THEN
         NPEROD = 1
         W(1) = 0.0D0
         W(IWW-1) = 0.0D0
      END IF
      CALL D03FAZ(LBDCND,LUNK,C1,MBDCND,MUNK,C2,NPEROD,NUNK,W,W(IWB),
     *            W(IWC),LDIMF,MDIMF,F(LSTART,MSTART,NSTART),IR,W(IWW),
     *            JWORK1,JWORK2,JTRIGL,JTRIGM)
C
C     Fill in sides for periodic boundary conditions.
C
      IF (LBDCND.EQ.0) THEN
         IF (MBDCND.EQ.0) THEN
            DO 740 K = NSTART, NSTOP
               F(1,M+1,K) = F(1,1,K)
  740       CONTINUE
            MSTOP = M + 1
         END IF
         IF (NBDCND.EQ.0) THEN
            DO 760 J = MSTART, MSTOP
               F(1,J,N+1) = F(1,J,1)
  760       CONTINUE
            NSTOP = N + 1
         END IF
         DO 800 J = MSTART, MSTOP
            DO 780 K = NSTART, NSTOP
               F(L+1,J,K) = F(1,J,K)
  780       CONTINUE
  800    CONTINUE
      END IF
      IF (MBDCND.EQ.0) THEN
         IF (NBDCND.EQ.0) THEN
            DO 820 I = LSTART, LSTOP
               F(I,1,N+1) = F(I,1,1)
  820       CONTINUE
            NSTOP = N + 1
         END IF
         DO 860 I = LSTART, LSTOP
            DO 840 K = NSTART, NSTOP
               F(I,M+1,K) = F(I,1,K)
  840       CONTINUE
  860    CONTINUE
      END IF
      IF (NBDCND.EQ.0) THEN
         DO 900 J = MSTART, MSTOP
            DO 880 I = LSTART, LSTOP
               F(I,J,N+1) = F(I,J,1)
  880       CONTINUE
  900    CONTINUE
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** XS must be less than XF',/' ** XS = ',F13.5,
     *       ', XF = ',F13.5)
99998 FORMAT (' ** L must be at least 5, L = ',I10)
99997 FORMAT (' ** LBDCND must satisfy 0 .LE. LBDCND .LE. 4           ',
     *       '  ',/' ** LBDCND = ',I10)
99996 FORMAT (' ** YS must be less than YF ',/' ** YS = ',F13.5,', YF ',
     *       '= ',F13.5)
99995 FORMAT (' ** M must be at least 5, M = ',I10)
99994 FORMAT (' ** MBDCND must satisfy 0 .LE. MBDCND .LE. 4           ',
     *       '  ',/' ** MBDCND = ',I10)
99993 FORMAT (' ** ZS must be less than ZF',/' ** ZS = ',F13.5,
     *       ', ZF = ',F13.5)
99992 FORMAT (' ** N must be at least 5, N = ',I10)
99991 FORMAT (' ** NBDCND must satisfy 0 .LE. NBDCND .LE. 4           ',
     *       '  ',/' ** NBDCND = ',I10)
99990 FORMAT (' ** LDIMF must be at least L+1',/'  ** LDIMF = ',I10,
     *       ', L = ',I10)
99989 FORMAT (' ** MDIMF must be at least M+1',/' ** MDIMF = ',I10,', ',
     *       'M = ',I10)
99988 FORMAT (' ** Workspace array W is too small',/' ** Routine requi',
     *       'res ',I10,' elements. You have supplied ',I10,' elements')
99987 FORMAT (' ** LAMBDA .GT. 0 in Helmholtz''s equation - a s','olut',
     *       'ion may not exist')
      END
