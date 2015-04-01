      SUBROUTINE G01NBF(CASE,MEAN,N,A,LDA,B,LDB,C,LDC,ELA,EMU,SIGMA,
     *                  LDSIG,L1,L2,LMAX,RMOM,ABSERR,EPS,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C       Works out L1-th to L2-th moments relating to
C       ratios of quadratic forms in normal variables.
C       Based on routine QRMOM by Jan Magnus and Bahran Pesaran
C
C       Maximum number of moments set by ISDIM (and ISPAR) is 12.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01NBF')
      INTEGER           ISDIM, ISPAR
      PARAMETER         (ISDIM=12,ISPAR=77)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ABSERR, EPS
      INTEGER           IFAIL, L1, L2, LDA, LDB, LDC, LDSIG, LMAX, N
      CHARACTER         CASE, MEAN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), B(LDA,N), C(LDC,*), ELA(*), EMU(*),
     *                  RMOM(L2-L1+1), SIGMA(LDSIG,N),
     *                  WK(3*N*N+(8+L2)*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSE, EPSM, FACT, RES, RL1
      INTEGER           I, IERROR, IFAULT, II, IRANK, IS, ISROW, ITEM,
     *                  JJ, LU, LWK1, LWK2, LWKA, LWKC, LWKR, LWKSLA,
     *                  LWKSMU, NI, NREC
C     .. Local Arrays ..
      INTEGER           ISPRTN(ISPAR*ISDIM+10)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01NBV, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G01NBV, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DGEMV, DPPTRF, DSYMV, DTPMV, DTPSV,
     *                  F02ABF, G01NAZ, G01NBT, G01NBY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, MAX, SQRT
C     .. Executable Statements ..
      EPSM = X02AJF()
      NREC = 1
C
C       Check the inputs
C
      IERROR = 1
      IF (N.LE.1) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (LDA.LT.N) THEN
         WRITE (P01REC(1),FMT=99997) LDA, N
      ELSE IF (LDB.LT.N) THEN
         WRITE (P01REC(1),FMT=99996) LDB, N
      ELSE IF (LDSIG.LT.N) THEN
         WRITE (P01REC(1),FMT=99994) LDSIG, N
      ELSE IF ((CASE.EQ.'Q' .OR. CASE.EQ.'q') .AND. LDC.LT.N) THEN
         WRITE (P01REC(1),FMT=99995) LDC, N
      ELSE IF (L1.LT.1) THEN
         WRITE (P01REC(1),FMT=99993) L1
      ELSE IF (L2.LT.L1) THEN
         WRITE (P01REC(1),FMT=99992) L1, L2
      ELSE IF (L2.GT.ISDIM) THEN
         WRITE (P01REC(1),FMT=99991) L2
      ELSE IF (CASE.NE.'R' .AND. CASE.NE.'r' .AND. CASE.NE.'L' .AND.
     *         CASE.NE.'l' .AND. CASE.NE.'Q' .AND. CASE.NE.'q') THEN
         WRITE (P01REC(1),FMT=99990) CASE
      ELSE IF (MEAN.NE.'M' .AND. MEAN.NE.'m' .AND. MEAN.NE.'Z' .AND.
     *         MEAN.NE.'z') THEN
         WRITE (P01REC(1),FMT=99989) MEAN
      ELSE IF (EPS.NE.0.0D0 .AND. EPS.LT.EPSM) THEN
         WRITE (P01REC(1),FMT=99998) EPS
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         IF (EPS.EQ.0.0D0) THEN
            EPSM = SQRT(EPSM)
         ELSE
            EPSM = EPS
         END IF
C
C        Partition workspace
C
         LWK1 = N*N + 1
         LWK2 = LWK1 + N*N
         LWKA = LWK2 + (4+L2)*N
         LWKC = LWKA + N*(N+1)/2
         LWKR = LWKC + N*(N+1)/2
         LWKSLA = LWKR + N
         LWKSMU = LWKSLA + N
C
C        CHECK THE EXISTENCE OF EXPECTATIONS AND
C        INITIALIZE ALL MATRICES AND VECTORS
C
         II = 0
         DO 20 I = 1, N
            NI = N - I + 1
            CALL DCOPY(NI,SIGMA(I,I),1,WK(LWKC+II),1)
            II = II + NI
   20    CONTINUE
C
C        Computes Choleski factor of SIGMA, L
C
         IFAULT = 1
         CALL DPPTRF('L',N,WK(LWKC),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99986)
            GO TO 240
         END IF
C
C        Check B is positive semi-definite
C
         IFAULT = 1
         CALL F02ABF(B,LDB,N,WK(LWKR),WK,N,WK(LWK2),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 5
            WRITE (P01REC(1),FMT=99988)
            GO TO 240
         END IF
         IRANK = N
         RL1 = WK(LWKR)
         DO 40 I = 1, N
            IF (WK(LWKR+I-1).GT.EPSM) GO TO 60
            WK(LWKR+I-1) = ZERO
            IRANK = IRANK - 1
   40    CONTINUE
   60    CONTINUE
         IF (IRANK.EQ.0 .OR. RL1.LT.-EPSM) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99987)
            GO TO 240
         END IF
C
C        Check existence of moments
C
         CALL G01NBY(CASE,N,A,LDA,C,LDC,ELA,WK,WK(LWK1),IRANK,ITEM,LMAX)
         IF (LMAX.LT.0) THEN
            LU = L2
            LMAX = L2
         ELSE IF (L1.GT.LMAX) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99983) LMAX, L1
            GO TO 240
         ELSE IF (L2.GT.LMAX) THEN
            IERROR = 6
            WRITE (P01REC(1),FMT=99982) LMAX, L2
            LU = LMAX
         ELSE
            LU = L2
            LMAX = L2
         END IF
C
C        Check L'BL = PDP' is positive semi-definite
C
         II = 0
         JJ = LWK1
         DO 80 I = 1, N
            NI = N - I + 1
            CALL DSYMV('L',NI,1.0D0,B(I,I),LDB,WK(LWKC+II),1,0.0D0,
     *                 WK(JJ),1)
            CALL DTPMV('L','T','N',NI,WK(LWKC+II),WK(JJ),1)
            II = II + NI
            JJ = JJ + N + 1
   80    CONTINUE
         IFAULT = 1
         CALL F02ABF(WK(LWK1),N,N,WK(LWKR),WK,N,WK(LWK2),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 5
            WRITE (P01REC(1),FMT=99988)
            GO TO 240
         END IF
         IRANK = N
         RL1 = WK(LWKR)
         DO 100 I = 1, N
            IF (WK(LWKR+I-1).GT.EPSM) GO TO 120
            WK(LWKR+I-1) = ZERO
            IRANK = IRANK - 1
  100    CONTINUE
  120    CONTINUE
         IF (IRANK.EQ.0 .OR. RL1.LT.-EPSM) THEN
            IERROR = 4
            WRITE (P01REC(1),FMT=99985)
            GO TO 240
         END IF
C
C         Compute: A* = P'L'ALP
C                  C* = P'L'CLP
C                  a* = P'L'a
C                 mu* = P'INV(L)mu
C
         DO 140 I = 1, N
            CALL DCOPY(N,WK((I-1)*N+1),1,WK(LWK1+(I-1)*N),1)
            CALL DTPMV('L','N','N',N,WK(LWKC),WK(LWK1+(I-1)*N),1)
  140    CONTINUE
         II = 0
         DO 160 I = 1, N
            CALL DSYMV('L',N,1.0D0,A,LDA,WK(LWK1+(I-1)*N),1,0.0D0,
     *                 WK(LWK2),1)
            CALL DGEMV('T',N,I,1.0D0,WK(LWK1),N,WK(LWK2),1,0.0D0,
     *                 WK(LWKA+II),1)
            II = II + I
  160    CONTINUE
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            CALL DCOPY(N,EMU,1,WK(LWK2),1)
            CALL DTPSV('L','N','N',N,WK(LWKC),WK(LWK2),1)
            CALL DGEMV('T',N,N,1.0D0,WK,N,WK(LWK2),1,0.0D0,WK(LWKSMU),1)
         END IF
         IF (CASE.EQ.'L' .OR. CASE.EQ.'l') THEN
            ISPRTN(1) = 2
            CALL DGEMV('T',N,N,1.0D0,WK(LWK1),N,ELA,1,0.0D0,WK(LWKSLA),
     *                 1)
         ELSE IF (CASE.EQ.'Q' .OR. CASE.EQ.'q') THEN
            ISPRTN(1) = 3
            II = 0
            DO 180 I = 1, N
               CALL DSYMV('L',N,1.0D0,C,LDC,WK(LWK1+(I-1)*N),1,0.0D0,
     *                    WK(LWK2),1)
               CALL DGEMV('T',N,I,1.0D0,WK(LWK1),N,WK(LWK2),1,0.0D0,
     *                    WK(LWKC+II),1)
               II = II + I
  180       CONTINUE
         ELSE
            ISPRTN(1) = 1
         END IF
C
C        Initialize FACT
C
         FACT = 1.0D0
         DO 200 I = 2, L1 - 1
            FACT = FACT*DBLE(I)
  200    CONTINUE
C
C        Put information into array ISPRTN
C
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            ISPRTN(2) = 1
         ELSE
            ISPRTN(2) = 0
         END IF
         ISPRTN(3) = N
         ISPRTN(6) = LWKR
         ISPRTN(7) = LWKSMU
         ISPRTN(8) = LWKSLA
         ISPRTN(9) = LWKA
         ISPRTN(10) = LWKC
C
C        WORK OUT L1-TH TO LU-TH MOMENTS
C
         ABSERR = 0.0D0
         DO 220 IS = L1, LU
            ISPRTN(4) = IS
C
C           WORK OUT ALL PARTITIONS OF IS
C
            CALL G01NAZ(IS,ISPRTN(11),ISPAR,ISROW,IFAULT)
            ISPRTN(5) = ISROW
C
C           CALCULATE THE IS-TH MOMENT AS AN INTEGRAL
C
            CALL G01NBT(G01NBV,0.0D0,1,0.0D0,EPSM,RES,ABSE,IFAULT,
     *                  ISPRTN,WK)
C
C           STORE RESULTS
C
            RMOM(IS-L1+1) = RES/FACT
            FACT = FACT*DBLE(IS)
            ABSE = ABSE/FACT
            ABSERR = MAX(ABSERR,ABSE)
  220    CONTINUE
         IF (IFAULT.NE.0) THEN
            IERROR = 7
            WRITE (P01REC,FMT=99984)
         END IF
      END IF
  240 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.le.1: N = ',I16)
99998 FORMAT (' ** On entry, EPS.lt.machine precision: EPS = ',D13.5)
99997 FORMAT (' ** On entry, LDA.lt.N: LDA = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, LDB.lt.N: LDB = ',I16,' N = ',I16)
99995 FORMAT (' ** On entry, LDC.lt.N: LDC = ',I16,' N = ',I16)
99994 FORMAT (' ** On entry, LDSIG.lt.N: LDSIG = ',I16,' N = ',I16)
99993 FORMAT (' ** On entry, L1.lt.1: L1 = ',I16)
99992 FORMAT (' ** On entry, L2.lt.L1: L1 = ',I16,' L2 = ',I16)
99991 FORMAT (' ** On entry, L2.gt.12: L2 = ',I16)
99990 FORMAT (' ** On entry, CASE is not valid: CASE = ',A1)
99989 FORMAT (' ** On entry, MEAN is not valid: MEAN = ',A1)
99988 FORMAT (' ** Failure in computing eigenvalues')
99987 FORMAT (' ** On entry, B is not positive semi-definite or is null'
     *       )
99986 FORMAT (' ** On entry, SIGMA is not positive-definite')
99985 FORMAT (' ** The matrix L''BL is not positive semi-definite or i',
     *       's null')
99984 FORMAT (' ** Full accuracy not achieved in integration')
99983 FORMAT (' ** Only ',I16,' moments exist, less than L1 = ',I16)
99982 FORMAT (' ** Only ',I16,' moments exist, less than L2 = ',I16)
      END
