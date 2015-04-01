      SUBROUTINE G01NAF(MOM,MEAN,N,A,LDA,EMU,SIGMA,LDSIG,L,RKUM,RMOM,WK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C       Calculates cumulants and moments of a
C       quadratic form in normal variables
C       Based on routine CUM by Jan Magnus and Bahram Pesaran
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01NAF')
      INTEGER           ISDIM, ISPAR
      PARAMETER         (ISDIM=12,ISPAR=77)
      DOUBLE PRECISION  ONE, TWO, ZERO
      PARAMETER         (ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, L, LDA, LDSIG, N
      CHARACTER         MEAN, MOM
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), EMU(*), RKUM(L), RMOM(*),
     *                  SIGMA(LDSIG,N), WK(3*N*(N+1)/2+N)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, F, FACTIS, R2L, RA, RB, RS, SUM, TP
      INTEGER           I, IERROR, IFAULT, II, IS, ISROW, J, K, LWK2,
     *                  LWK3, LWK4, NI, NN, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  TR(ISDIM)
      INTEGER           ISPRTN(ISPAR,ISDIM)
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G02AAR, X02AJF
      INTEGER           P01ABF
      EXTERNAL          G02AAR, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, DPPTRF, DSYMV, DTPMV, DTPSV, G01NAZ,
     *                  G02AAS, G02AAX
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (N.LE.1) THEN
         WRITE (P01REC(1),FMT=99999) N
      ELSE IF (L.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) L
      ELSE IF (L.GT.ISDIM) THEN
         WRITE (P01REC(1),FMT=99995) L
      ELSE IF (LDA.LT.N) THEN
         WRITE (P01REC(1),FMT=99994) LDA, N
      ELSE IF (LDSIG.LT.N) THEN
         WRITE (P01REC(1),FMT=99993) LDSIG, N
      ELSE IF (MOM.NE.'M' .AND. MOM.NE.'m' .AND. MOM.NE.'C' .AND.
     *         MOM.NE.'c') THEN
         WRITE (P01REC(1),FMT=99997) MOM
      ELSE IF (MEAN.NE.'M' .AND. MEAN.NE.'m' .AND. MEAN.NE.'Z' .AND.
     *         MEAN.NE.'z') THEN
         WRITE (P01REC(1),FMT=99996) MEAN
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.NE.1) THEN
         NN = N*(N+1)/2
         LWK2 = NN
         LWK3 = LWK2 + NN
         LWK4 = LWK3 + NN
         EPS = X02AJF()
C
C        SET WK2 = L WHERE SIGMA = LL'
C
         II = LWK2 + 1
         DO 20 I = 1, N
            NI = N - I + 1
            CALL DCOPY(NI,SIGMA(I,I),1,WK(II),1)
            II = II + NI
   20    CONTINUE
         CALL DPPTRF('L',N,WK(LWK2+1),IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 2
            WRITE (P01REC(1),FMT=99992)
            GO TO 180
         END IF
C
C        Set up L'AL in work
C
         II = 1
         DO 40 I = 1, N
            NI = N - I + 1
            CALL DSYMV('L',NI,1.0D0,A(I,I),LDA,WK(LWK2+II),1,0.0D0,
     *                 WK(II),1)
            CALL DTPMV('L','T','N',NI,WK(LWK2+II),WK(II),1)
            II = II + NI
   40    CONTINUE
C
C        Calculate LV*EMU if MEAN = 'M'
C
         IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            CALL DCOPY(N,EMU,1,WK(LWK4+1),1)
            CALL DTPSV('L','N','N',N,WK(LWK2+1),WK(LWK4+1),1)
         END IF
C
C        Calculate cumulants of quadratic form in RKUM
C
         FACTIS = 1.0D0
         TP = 1.0D0
         DO 160 IS = 1, L
            IF (IS.EQ.1) THEN
               CALL G02AAX('L',N,WK,WK(LWK2+1))
               CALL DCOPY(NN,WK(LWK2+1),1,WK,1)
            ELSE
               CALL G02AAS('U',N,WK,WK(LWK2+1),WK(LWK3+1))
               CALL DCOPY(NN,WK(LWK3+1),1,WK(LWK2+1),1)
            END IF
            TR(IS) = G02AAR('U',N,WK(LWK2+1))
            IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
               SUM = ZERO
               II = LWK2
               DO 80 I = 1, N
                  DO 60 J = 1, I - 1
                     II = II + 1
                     SUM = SUM + TWO*WK(LWK4+I)*WK(II)*WK(LWK4+J)
   60             CONTINUE
                  II = II + 1
                  SUM = SUM + WK(LWK4+I)*WK(II)*WK(LWK4+I)
   80          CONTINUE
               RS = IS
               TR(IS) = TR(IS) + RS*SUM
            END IF
            RKUM(IS) = TP*FACTIS*TR(IS)
            TP = TP*2.0D0
            FACTIS = FACTIS*DBLE(IS)
C
C           Calculate moments of quadratic form (if required) in RMOM
C
            IF (MOM.EQ.'M' .OR. MOM.EQ.'m') THEN
               CALL G01NAZ(IS,ISPRTN,ISPAR,ISROW,IFAULT)
               SUM = ZERO
               DO 140 I = 1, ISROW
                  RB = ONE
                  RA = ONE
                  DO 120 J = 1, IS
                     R2L = 2*J
                     K = ISPRTN(I,J)
                     IF (K.NE.0) THEN
                        F = R2L
                        DO 100 II = 2, K
                           F = F*R2L*DBLE(II)
  100                   CONTINUE
                        RB = RB*F
                        RA = RA*(TR(J)**K)
                     END IF
  120             CONTINUE
                  RB = (FACTIS*(TWO**IS))/RB
                  SUM = SUM + RA*RB
  140          CONTINUE
               RMOM(IS) = SUM
            END IF
  160    CONTINUE
      END IF
  180 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N.le.1: N = ',I16)
99998 FORMAT (' ** On entry, L.lt.1: L = ',I16)
99997 FORMAT (' ** On entry, MOM is not valid: MOM = ',A1)
99996 FORMAT (' ** On entry, MEAN is not valid: MEAN = ',A1)
99995 FORMAT (' ** On entry, L.gt.12: L = ',I16)
99994 FORMAT (' ** On entry, LDA.lt.N: LDA = ',I16,' N = ',I16)
99993 FORMAT (' ** On entry, LDSIG.lt.N: LDSIG = ',I16,' N = ',I16)
99992 FORMAT (' ** On entry, SIGMA is not positive-definite')
      END
