      SUBROUTINE G03ADF(WEIGHT,N,M,Z,LDZ,ISZ,NX,NY,WT,E,LDE,NCV,CVX,
     *                  LDCVX,MCV,CVY,LDCVY,TOL,WK,IWK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03ADF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IFAIL, IWK, LDCVX, LDCVY, LDE, LDZ, M, MCV, N,
     *                  NCV, NX, NY
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  CVX(LDCVX,MCV), CVY(LDCVY,MCV), E(LDE,6),
     *                  WK(IWK), WT(*), Z(LDZ,M)
      INTEGER           ISZ(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, RDF, WSUM
      INTEGER           I, IERROR, IFAULT, IRANKX, IRANKY, IWKMIN, IXY,
     *                  K, LQY, LWK, MXY, NREC, NWT, NXN
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G03AAZ, G03ACZ, G03ADZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, DBLE, SQRT
C     .. Executable Statements ..
      NREC = 1
      IERROR = 1
      IF (NX.GE.NY) THEN
         LQY = MAX(NY*N,NX*NX+5*(NX-1))
         IWKMIN = N*NX + NX + NY + LQY
         MXY = NY
      ELSE
         LQY = MAX(NX*N,NY*NY+5*(NY-1))
         IWKMIN = N*NY + NY + NX + LQY
         MXY = NX
      END IF
      IXY = NX + NY
      IF (NX.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) NX
      ELSE IF (NY.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) NY
      ELSE IF (M.LT.IXY) THEN
         WRITE (P01REC(1),FMT=99997) M, IXY
      ELSE IF (N.LE.IXY) THEN
         WRITE (P01REC(1),FMT=99996) N, IXY
      ELSE IF (LDZ.LT.N) THEN
         WRITE (P01REC(1),FMT=99995) LDZ, N
      ELSE IF (LDCVX.LT.NX) THEN
         WRITE (P01REC(1),FMT=99994) LDCVX, NX
      ELSE IF (LDCVY.LT.NY) THEN
         WRITE (P01REC(1),FMT=99993) LDCVY, NY
      ELSE IF (MCV.LT.MXY) THEN
         NREC = 2
         WRITE (P01REC,FMT=99980) MCV, MXY
      ELSE IF (LDE.LT.MXY) THEN
         NREC = 2
         WRITE (P01REC,FMT=99992) LDE, MXY
      ELSE IF (IWK.LT.IWKMIN) THEN
         NREC = 2
         WRITE (P01REC,FMT=99991) IWKMIN, IWK
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         WRITE (P01REC(1),FMT=99988) WEIGHT
      ELSE IF (TOL.LT.0.0D0) THEN
         WRITE (P01REC(1),FMT=99985) TOL
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.0) THEN
         EPS = X02AJF()
         IF (TOL.LT.EPS) THEN
            EPS = SQRT(EPS)
         ELSE
            EPS = TOL
         END IF
         K = 0
         DO 20 I = 1, M
            IF (ISZ(I).GT.0) K = K + 1
   20    CONTINUE
         IF (K.NE.NX) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99982) K, NX
            GO TO 120
         END IF
         K = 0
         DO 40 I = 1, M
            IF (ISZ(I).LT.0) K = K + 1
   40    CONTINUE
         IF (K.NE.NY) THEN
            IERROR = 3
            WRITE (P01REC(1),FMT=99981) K, NY
            GO TO 120
         END IF
C
C        CHECK WEIGHTS
C
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            WSUM = 0.0D0
            NWT = (IXY-1)*N
            DO 60 I = 1, N
               IF (WT(I).LT.0.0D0) GO TO 80
               IF (WT(I).GT.0.0D0) THEN
                  WSUM = WSUM + WT(I)
                  WK(NWT+I) = SQRT(WT(I))
               ELSE
                  WK(NWT+I) = 0.0D0
               END IF
   60       CONTINUE
            GO TO 100
   80       IERROR = 2
            WRITE (P01REC(1),FMT=99989) I
            GO TO 120
  100       CONTINUE
            IF (WSUM.LT.DBLE(IXY+1)) THEN
               IERROR = 4
               WRITE (P01REC(1),FMT=99984)
               GO TO 120
            END IF
         ELSE
            WSUM = DBLE(N)
         END IF
         RDF = SQRT(WSUM-1.0D0)
         IF (NX.GE.NY) THEN
            NXN = N*NX + 1
            LWK = NXN + LQY
            CALL G03AAZ('U',WEIGHT,N,Z,LDZ,M,ISZ,IXY,WT,WSUM,WK,N,E,
     *                  WK(LWK))
            CALL G03AAZ('N',WEIGHT,N,Z,LDZ,M,ISZ,NY,WT,WSUM,WK(NXN),N,E,
     *                  WK(LWK))
            CALL G03ADZ(N,NX,NY,CVX,LDCVX,CVY,LDCVY,WK,WK(NXN),LQY,RDF,
     *                  EPS,IRANKX,IRANKY,NCV,E,WK(LWK),IFAULT)
            IF (IFAULT.EQ.1) THEN
               IERROR = 5
               WRITE (P01REC(1),FMT=99990)
               GO TO 120
            ELSE IF (IFAULT.EQ.2) THEN
               IERROR = 7
               WRITE (P01REC(1),FMT=99987)
               GO TO 120
            ELSE IF (IFAULT.EQ.3) THEN
               IERROR = 7
               WRITE (P01REC(1),FMT=99986)
               GO TO 120
            END IF
         ELSE
            NXN = N*NY + 1
            LWK = NXN + LQY
            CALL G03AAZ('U',WEIGHT,N,Z,LDZ,M,ISZ,NX,WT,WSUM,WK(NXN),N,E,
     *                  WK(LWK))
            CALL G03AAZ('N',WEIGHT,N,Z,LDZ,M,ISZ,IXY,WT,WSUM,WK,N,E,
     *                  WK(LWK))
            CALL G03ADZ(N,NY,NX,CVY,LDCVY,CVX,LDCVX,WK,WK(NXN),LQY,RDF,
     *                  EPS,IRANKY,IRANKX,NCV,E,WK(LWK),IFAULT)
            IF (IFAULT.EQ.1) THEN
               IERROR = 5
               WRITE (P01REC(1),FMT=99990)
               GO TO 120
            ELSE IF (IFAULT.EQ.3) THEN
               IERROR = 7
               WRITE (P01REC(1),FMT=99987)
               GO TO 120
            ELSE IF (IFAULT.EQ.2) THEN
               IERROR = 7
               WRITE (P01REC(1),FMT=99986)
               GO TO 120
            END IF
         END IF
C
C        CALCULATE TEST STATISTICS
C
         CALL G03ACZ(E,LDE,WSUM,NCV,IRANKX,IRANKY,IFAULT)
         IF (IFAULT.EQ.1) THEN
            IERROR = 6
            WRITE (P01REC(1),FMT=99983)
         END IF
         NX = IRANKX
         NY = IRANKY
      END IF
  120 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (' ** On entry, NX.lt.1 : NX = ',I16)
99998 FORMAT (' ** On entry, NY.lt.1 : NY = ',I16)
99997 FORMAT (' ** On entry, M.lt.NX+NY : M = ',I16,' NX+NY = ',I16)
99996 FORMAT (' ** On entry, N.le.NX+NY : N = ',I16,' NX+NY = ',I16)
99995 FORMAT (' ** On entry, LDZ.lt.N : LDZ = ',I16,' N = ',I16)
99994 FORMAT (' ** On entry, LDCVX.lt.NX : LDCVX = ',I16,' NX = ',I16)
99993 FORMAT (' ** On entry, LDCVY.lt.NY : LDCVY = ',I16,' NY = ',I16)
99992 FORMAT (' ** On entry, LDE.lt.MIN(NX,NY) : LDE = ',I16,/'       ',
     *       '                    MIN(NX,NY) = ',I16)
99991 FORMAT (' ** On entry, IWK is too small, minimum value = ',I16,
     *       /'                                          IWK = ',I16)
99990 FORMAT (' ** An SVD has failed to converge')
99989 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
99988 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99987 FORMAT (' ** Rank of X matrix is 0')
99986 FORMAT (' ** Rank of Y matrix is 0')
99985 FORMAT (' ** On entry, TOL.lt.0.0 : TOL = ',D13.5)
99984 FORMAT (' ** Effective number of observations is less than NX+NY',
     *       '+1')
99983 FORMAT (' ** Canonical correlation equal to 1.0')
99982 FORMAT (' ** On entry, there are ',I16,' X vars instead of ',I16)
99981 FORMAT (' ** On entry, there are ',I16,' Y vars instead of ',I16)
99980 FORMAT (' ** On entry, MCV.lt.MIN(NX,NY) : MCV = ',I16,/'       ',
     *       '                    MIN(NX,NY) = ',I16)
      END
