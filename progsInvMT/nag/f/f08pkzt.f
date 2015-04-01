      SUBROUTINE F08PKZ(RIGHTV,NOINIT,N,H,LDH,WR,WI,VR,VI,B,LDB,WORK,
     *                  EPS3,SMLNUM,BIGNUM,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLAEIN(RIGHTV,NOINIT,N,H,LDH,WR,WI,VR,VI,B,LDB,
C    *                  WORK,EPS3,SMLNUM,BIGNUM,INFO)
C
C  Purpose
C  =======
C
C  DLAEIN uses inverse iteration to find a right or left eigenvector
C  corresponding to the eigenvalue (WR,WI) of a real upper Hessenberg
C  matrix H.
C
C  Arguments
C  =========
C
C  RIGHTV   (input) LOGICAL
C          = .TRUE. : compute right eigenvector;
C          = .FALSE.: compute left eigenvector.
C
C  NOINIT   (input) LOGICAL
C          = .TRUE. : no initial vector supplied in (VR,VI).
C          = .FALSE.: initial vector supplied in (VR,VI).
C
C  N       (input) INTEGER
C          The order of the matrix H.  N >= 0.
C
C  H       (input) DOUBLE PRECISION array, dimension (LDH,N)
C          The upper Hessenberg matrix H.
C
C  LDH     (input) INTEGER
C          The leading dimension of the array H.  LDH >= max(1,N).
C
C  WR      (input) DOUBLE PRECISION
C  WI      (input) DOUBLE PRECISION
C          The real and imaginary parts of the eigenvalue of H whose
C          corresponding right or left eigenvector is to be computed.
C
C  VR      (input/output) DOUBLE PRECISION array, dimension (N)
C  VI      (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, if NOINIT = .FALSE. and WI = 0.0, VR must contain
C          a real starting vector for inverse iteration using the real
C          eigenvalue WR; if NOINIT = .FALSE. and WI.ne.0.0, VR and VI
C          must contain the real and imaginary parts of a complex
C          starting vector for inverse iteration using the complex
C          eigenvalue (WR,WI); otherwise VR and VI need not be set.
C          On exit, if WI = 0.0 (real eigenvalue), VR contains the
C          computed real eigenvector; if WI.ne.0.0 (complex eigenvalue),
C          VR and VI contain the real and imaginary parts of the
C          computed complex eigenvector. The eigenvector is normalized
C          so that the component of largest magnitude has magnitude 1;
C          here the magnitude of a complex number (x,y) is taken to be
C          |x| + |y|.
C          VI is not referenced if WI = 0.0.
C
C  B       (workspace) DOUBLE PRECISION array, dimension (LDB,N)
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= N+1.
C
C  WORK   (workspace) DOUBLE PRECISION array, dimension (N)
C
C  EPS3    (input) DOUBLE PRECISION
C          A small machine-dependent value which is used to perturb
C          close eigenvalues, and to replace zero pivots.
C
C  SMLNUM  (input) DOUBLE PRECISION
C          A machine-dependent value close to the underflow threshold.
C
C  BIGNUM  (input) DOUBLE PRECISION
C          A machine-dependent value close to the overflow threshold.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          = 1:  inverse iteration did not converge; VR is set to the
C                last iterate, and so is VI if WI.ne.0.0.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TENTH
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,TENTH=1.0D-1)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIGNUM, EPS3, SMLNUM, WI, WR
      INTEGER           INFO, LDB, LDH, N
      LOGICAL           NOINIT, RIGHTV
C     .. Array Arguments ..
      DOUBLE PRECISION  B(LDB,*), H(LDH,*), VI(*), VR(*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSBII, ABSBJJ, EI, EJ, GROWTO, NORM, NRMSML,
     *                  REC, ROOTN, SCALE, TEMP, VCRIT, VMAX, VNORM, W,
     *                  W1, X, XI, XR, Y
      INTEGER           I, I1, I2, I3, IERR, ITS, J
      CHARACTER         NORMIN, TRANS
C     .. External Functions ..
      DOUBLE PRECISION  DASUM, DNRM2, F06BNF
      INTEGER           IDAMAX
      EXTERNAL          DASUM, DNRM2, F06BNF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          A02ACF, DROT, DSCAL, F07TGZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, SQRT
C     .. Executable Statements ..
C
      INFO = 0
C
C     GROWTO is the threshold used in the acceptance test for an
C     eigenvector.
C
      ROOTN = SQRT(DBLE(N))
      GROWTO = TENTH/ROOTN
      NRMSML = MAX(ONE,EPS3*ROOTN)*SMLNUM
C
C     Form B = H - (WR,WI)*I (except that the subdiagonal elements and
C     the imaginary parts of the diagonal elements are not stored).
C
      DO 40 J = 1, N
         DO 20 I = 1, J - 1
            B(I,J) = H(I,J)
   20    CONTINUE
         B(J,J) = H(J,J) - WR
   40 CONTINUE
C
      IF (WI.EQ.ZERO) THEN
C
C        Real eigenvalue.
C
         IF (NOINIT) THEN
C
C           Set initial vector.
C
            DO 60 I = 1, N
               VR(I) = EPS3
   60       CONTINUE
         ELSE
C
C           Scale supplied initial vector.
C
            VNORM = DNRM2(N,VR,1)
            CALL DSCAL(N,(EPS3*ROOTN)/MAX(VNORM,NRMSML),VR,1)
         END IF
C
         IF (RIGHTV) THEN
C
C           LU decomposition with partial pivoting of B, replacing zero
C           pivots by EPS3.
C
            DO 120 I = 1, N - 1
               EI = H(I+1,I)
               IF (ABS(B(I,I)).LT.ABS(EI)) THEN
C
C                 Interchange rows and eliminate.
C
                  X = B(I,I)/EI
                  B(I,I) = EI
                  DO 80 J = I + 1, N
                     TEMP = B(I+1,J)
                     B(I+1,J) = B(I,J) - X*TEMP
                     B(I,J) = TEMP
   80             CONTINUE
               ELSE
C
C                 Eliminate without interchange.
C
                  IF (B(I,I).EQ.ZERO) B(I,I) = EPS3
                  X = EI/B(I,I)
                  IF (X.NE.ZERO) THEN
                     DO 100 J = I + 1, N
                        B(I+1,J) = B(I+1,J) - X*B(I,J)
  100                CONTINUE
                  END IF
               END IF
  120       CONTINUE
            IF (B(N,N).EQ.ZERO) B(N,N) = EPS3
C
            TRANS = 'N'
C
         ELSE
C
C           UL decomposition with partial pivoting of B, replacing zero
C           pivots by EPS3.
C
            DO 180 J = N, 2, -1
               EJ = H(J,J-1)
               IF (ABS(B(J,J)).LT.ABS(EJ)) THEN
C
C                 Interchange columns and eliminate.
C
                  X = B(J,J)/EJ
                  B(J,J) = EJ
                  DO 140 I = 1, J - 1
                     TEMP = B(I,J-1)
                     B(I,J-1) = B(I,J) - X*TEMP
                     B(I,J) = TEMP
  140             CONTINUE
               ELSE
C
C                 Eliminate without interchange.
C
                  IF (B(J,J).EQ.ZERO) B(J,J) = EPS3
                  X = EJ/B(J,J)
                  IF (X.NE.ZERO) THEN
                     DO 160 I = 1, J - 1
                        B(I,J-1) = B(I,J-1) - X*B(I,J)
  160                CONTINUE
                  END IF
               END IF
  180       CONTINUE
            IF (B(1,1).EQ.ZERO) B(1,1) = EPS3
C
            TRANS = 'T'
C
         END IF
C
         NORMIN = 'N'
         DO 220 ITS = 1, N
C
C           Solve U*x = scale*v for a right eigenvector
C             or U'*x = scale*v for a left eigenvector,
C           overwriting x on v.
C
            CALL F07TGZ('Upper',TRANS,'Nonunit',NORMIN,N,B,LDB,VR,SCALE,
     *                  WORK,IERR)
            NORMIN = 'Y'
C
C           Test for sufficient growth in the norm of v.
C
            VNORM = DASUM(N,VR,1)
            IF (VNORM.GE.GROWTO*SCALE) GO TO 240
C
C           Choose new orthogonal starting vector and try again.
C
            TEMP = EPS3/(ROOTN+ONE)
            VR(1) = EPS3
            DO 200 I = 2, N
               VR(I) = TEMP
  200       CONTINUE
            VR(N-ITS+1) = VR(N-ITS+1) - EPS3*ROOTN
  220    CONTINUE
C
C        Failure to find eigenvector in N iterations.
C
         INFO = 1
C
  240    CONTINUE
C
C        Normalize eigenvector.
C
         I = IDAMAX(N,VR,1)
         CALL DSCAL(N,ONE/VR(I),VR,1)
      ELSE
C
C        Complex eigenvalue.
C
         IF (NOINIT) THEN
C
C           Set initial vector.
C
            DO 260 I = 1, N
               VR(I) = EPS3
               VI(I) = ZERO
  260       CONTINUE
         ELSE
C
C           Scale supplied initial vector.
C
            NORM = F06BNF(DNRM2(N,VR,1),DNRM2(N,VI,1))
            REC = (EPS3*ROOTN)/MAX(NORM,NRMSML)
            CALL DSCAL(N,REC,VR,1)
            CALL DSCAL(N,REC,VI,1)
         END IF
C
         IF (RIGHTV) THEN
C
C           LU decomposition with partial pivoting of B, replacing zero
C           pivots by EPS3.
C
C           The imaginary part of the (i,j)-th element of U is stored in
C           B(j+1,i).
C
            B(2,1) = -WI
            DO 280 I = 2, N
               B(I+1,1) = ZERO
  280       CONTINUE
C
            DO 340 I = 1, N - 1
               ABSBII = F06BNF(B(I,I),B(I+1,I))
               EI = H(I+1,I)
               IF (ABSBII.LT.ABS(EI)) THEN
C
C                 Interchange rows and eliminate.
C
                  XR = B(I,I)/EI
                  XI = B(I+1,I)/EI
                  B(I,I) = EI
                  B(I+1,I) = ZERO
                  DO 300 J = I + 1, N
                     TEMP = B(I+1,J)
                     B(I+1,J) = B(I,J) - XR*TEMP
                     B(J+1,I+1) = B(J+1,I) - XI*TEMP
                     B(I,J) = TEMP
                     B(J+1,I) = ZERO
  300             CONTINUE
                  B(I+2,I) = -WI
                  B(I+1,I+1) = B(I+1,I+1) - XI*WI
                  B(I+2,I+1) = B(I+2,I+1) + XR*WI
               ELSE
C
C                 Eliminate without interchanging rows.
C
                  IF (ABSBII.EQ.ZERO) THEN
                     B(I,I) = EPS3
                     B(I+1,I) = ZERO
                     ABSBII = EPS3
                  END IF
                  EI = (EI/ABSBII)/ABSBII
                  XR = B(I,I)*EI
                  XI = -B(I+1,I)*EI
                  DO 320 J = I + 1, N
                     B(I+1,J) = B(I+1,J) - XR*B(I,J) + XI*B(J+1,I)
                     B(J+1,I+1) = -XR*B(J+1,I) - XI*B(I,J)
  320             CONTINUE
                  B(I+2,I+1) = B(I+2,I+1) - WI
               END IF
C
C              Compute 1-norm of offdiagonal elements of i-th row.
C
               WORK(I) = DASUM(N-I,B(I,I+1),LDB) + DASUM(N-I,B(I+2,I),1)
  340       CONTINUE
            IF (B(N,N).EQ.ZERO .AND. B(N+1,N).EQ.ZERO) B(N,N) = EPS3
            WORK(N) = ZERO
C
            I1 = N
            I2 = 1
            I3 = -1
         ELSE
C
C           UL decomposition with partial pivoting of conjg(B),
C           replacing zero pivots by EPS3.
C
C           The imaginary part of the (i,j)-th element of U is stored in
C           B(j+1,i).
C
            B(N+1,N) = WI
            DO 360 J = 1, N - 1
               B(N+1,J) = ZERO
  360       CONTINUE
C
            DO 420 J = N, 2, -1
               EJ = H(J,J-1)
               ABSBJJ = F06BNF(B(J,J),B(J+1,J))
               IF (ABSBJJ.LT.ABS(EJ)) THEN
C
C                 Interchange columns and eliminate
C
                  XR = B(J,J)/EJ
                  XI = B(J+1,J)/EJ
                  B(J,J) = EJ
                  B(J+1,J) = ZERO
                  DO 380 I = 1, J - 1
                     TEMP = B(I,J-1)
                     B(I,J-1) = B(I,J) - XR*TEMP
                     B(J,I) = B(J+1,I) - XI*TEMP
                     B(I,J) = TEMP
                     B(J+1,I) = ZERO
  380             CONTINUE
                  B(J+1,J-1) = WI
                  B(J-1,J-1) = B(J-1,J-1) + XI*WI
                  B(J,J-1) = B(J,J-1) - XR*WI
               ELSE
C
C                 Eliminate without interchange.
C
                  IF (ABSBJJ.EQ.ZERO) THEN
                     B(J,J) = EPS3
                     B(J+1,J) = ZERO
                     ABSBJJ = EPS3
                  END IF
                  EJ = (EJ/ABSBJJ)/ABSBJJ
                  XR = B(J,J)*EJ
                  XI = -B(J+1,J)*EJ
                  DO 400 I = 1, J - 1
                     B(I,J-1) = B(I,J-1) - XR*B(I,J) + XI*B(J+1,I)
                     B(J,I) = -XR*B(J+1,I) - XI*B(I,J)
  400             CONTINUE
                  B(J,J-1) = B(J,J-1) + WI
               END IF
C
C              Compute 1-norm of offdiagonal elements of j-th column.
C
               WORK(J) = DASUM(J-1,B(1,J),1) + DASUM(J-1,B(J+1,1),LDB)
  420       CONTINUE
            IF (B(1,1).EQ.ZERO .AND. B(2,1).EQ.ZERO) B(1,1) = EPS3
            WORK(1) = ZERO
C
            I1 = 1
            I2 = N
            I3 = 1
         END IF
C
         DO 540 ITS = 1, N
            SCALE = ONE
            VMAX = ONE
            VCRIT = BIGNUM
C
C           Solve U*(xr,xi) = scale*(vr,vi) for a right eigenvector,
C             or U'*(xr,xi) = scale*(vr,vi) for a left eigenvector,
C           overwriting (xr,xi) on (vr,vi).
C
            DO 500 I = I1, I2, I3
C
               IF (WORK(I).GT.VCRIT) THEN
                  REC = ONE/VMAX
                  CALL DSCAL(N,REC,VR,1)
                  CALL DSCAL(N,REC,VI,1)
                  SCALE = SCALE*REC
                  VMAX = ONE
                  VCRIT = BIGNUM
               END IF
C
               XR = VR(I)
               XI = VI(I)
               IF (RIGHTV) THEN
                  DO 440 J = I + 1, N
                     XR = XR - B(I,J)*VR(J) + B(J+1,I)*VI(J)
                     XI = XI - B(I,J)*VI(J) - B(J+1,I)*VR(J)
  440             CONTINUE
               ELSE
                  DO 460 J = 1, I - 1
                     XR = XR - B(J,I)*VR(J) + B(I+1,J)*VI(J)
                     XI = XI - B(J,I)*VI(J) - B(I+1,J)*VR(J)
  460             CONTINUE
               END IF
C
               W = ABS(B(I,I)) + ABS(B(I+1,I))
               IF (W.GT.SMLNUM) THEN
                  IF (W.LT.ONE) THEN
                     W1 = ABS(XR) + ABS(XI)
                     IF (W1.GT.W*BIGNUM) THEN
                        REC = ONE/W1
                        CALL DSCAL(N,REC,VR,1)
                        CALL DSCAL(N,REC,VI,1)
                        XR = VR(I)
                        XI = VI(I)
                        SCALE = SCALE*REC
                        VMAX = VMAX*REC
                     END IF
                  END IF
C
C                 Divide by diagonal element of B.
C
                  CALL A02ACF(XR,XI,B(I,I),B(I+1,I),VR(I),VI(I))
                  VMAX = MAX(ABS(VR(I))+ABS(VI(I)),VMAX)
                  VCRIT = BIGNUM/VMAX
               ELSE
                  DO 480 J = 1, N
                     VR(J) = ZERO
                     VI(J) = ZERO
  480             CONTINUE
                  VR(I) = ONE
                  VI(I) = ONE
                  SCALE = ZERO
                  VMAX = ONE
                  VCRIT = BIGNUM
               END IF
  500       CONTINUE
C
C           Test for sufficient growth in the norm of (VR,VI).
C
            VNORM = DASUM(N,VR,1) + DASUM(N,VI,1)
            IF (VNORM.GE.GROWTO*SCALE) GO TO 560
C
C           Choose a new orthogonal starting vector and try again.
C
            Y = EPS3/(ROOTN+ONE)
            VR(1) = EPS3
            VI(1) = ZERO
C
            DO 520 I = 2, N
               VR(I) = Y
               VI(I) = ZERO
  520       CONTINUE
            VR(N-ITS+1) = VR(N-ITS+1) - EPS3*ROOTN
  540    CONTINUE
C
C        Failure to find eigenvector in N iterations
C
         INFO = 1
C
  560    CONTINUE
C
C        Normalize eigenvector.
C
         DO 580 I = 1, N
            WORK(I) = F06BNF(VR(I),VI(I))
  580    CONTINUE
         I = IDAMAX(N,WORK,1)
         XR = VR(I)/WORK(I)
         XI = VI(I)/WORK(I)
         CALL DROT(N,VR,1,VI,1,XR,XI)
         DO 600 I = 1, N
            WORK(I) = ABS(VR(I)) + ABS(VI(I))
  600    CONTINUE
         I = IDAMAX(N,WORK,1)
         VNORM = WORK(I)
         CALL DSCAL(N,ONE/VNORM,VR,1)
         CALL DSCAL(N,ONE/VNORM,VI,1)
C
      END IF
C
      RETURN
C
C     End of F08PKZ (DLAEIN)
C
      END
