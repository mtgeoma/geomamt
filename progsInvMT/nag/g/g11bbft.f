      SUBROUTINE G11BBF(TYPE,WEIGHT,N,NFAC,ISF,LFAC,IFAC,LDF,PERCNT,Y,
     *                  WT,TABLE,MAXT,NCELLS,NDIM,IDIM,ICOUNT,IWK,WK,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes table from data and factors using quantiles for precent
C     given by PERCNT.
C
C     TYPE indicates if the distribution is considered to be discrete,
C     i.e. CDF is a step function, or continuous with Y representing the
C     class mid-points, i.e. CDF is a polygon.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G11BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PERCNT
      INTEGER           IFAIL, LDF, MAXT, N, NCELLS, NDIM, NFAC
      CHARACTER         TYPE, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  TABLE(MAXT), WK(2*N), WT(*), Y(N)
      INTEGER           ICOUNT(MAXT), IDIM(NFAC), IFAC(LDF,NFAC),
     *                  ISF(NFAC), IWK(2*NFAC+N), LFAC(NFAC)
C     .. Local Scalars ..
      DOUBLE PRECISION  APT, CPR, CPRMAX, CPRP, F, SUM, WSUM
      INTEGER           I, ICS, IERROR, IFAULT, IL, IRANK, IU, J, K, M,
     *                  NOBS, NREC
      LOGICAL           DDN, EMPTY, WTL
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G04CAU, M01CAF, M01DAF, M01EAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT
C     .. Executable Statements ..
C
C     Carry out initial error checks
C
      NREC = 1
      IERROR = 1
      IF (N.LT.2) THEN
         WRITE (P01REC,FMT=99999) N
      ELSE IF (NFAC.LT.1) THEN
         WRITE (P01REC,FMT=99998) NFAC
      ELSE IF (LDF.LT.N) THEN
         WRITE (P01REC,FMT=99997) LDF, N
      ELSE IF (PERCNT.LE.0.0D0 .OR. PERCNT.GE.100) THEN
         WRITE (P01REC,FMT=99995) PERCNT
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u') THEN
         WRITE (P01REC,FMT=99991) WEIGHT
      ELSE IF (TYPE.NE.'D' .AND. TYPE.NE.'d' .AND. TYPE.NE.'C' .AND.
     *         TYPE.NE.'c') THEN
         WRITE (P01REC,FMT=99996) TYPE
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 380
      IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
         WTL = .TRUE.
      ELSE
         WTL = .FALSE.
      END IF
      IF (TYPE.EQ.'D' .OR. TYPE.EQ.'d') THEN
         DDN = .TRUE.
      ELSE
         DDN = .FALSE.
      END IF
C
C     Compute number of cells in the table.
C
      NDIM = 0
      DO 20 I = 1, NFAC
         IF (ISF(I).GT.0) NDIM = NDIM + 1
   20 CONTINUE
C
C     If no factors compute for single group
C
      IF (NDIM.EQ.0) THEN
         NDIM = 1
         IDIM(1) = 1
         NCELLS = 1
         IWK(1) = 0
         IF (WTL) THEN
            K = 0
            DO 40 I = 1, N
               IF (WT(I).GT.0.D0) THEN
                  K = K + 1
                  WK(K) = Y(I)
                  WK(N+K) = WT(I)
               END IF
   40       CONTINUE
            IWK(2) = K
         ELSE
            DO 60 I = 1, N
               WK(I) = Y(I)
   60       CONTINUE
            IWK(2) = N
         END IF
         ICOUNT(1) = 0
         NOBS = IWK(2)
      ELSE
         J = NDIM
         IWK(J) = 1
         K = NDIM
         DO 80 I = NFAC, 1, -1
            IF (ISF(I).GT.0) THEN
               IF (LFAC(I).LE.1) THEN
                  NREC = 1
                  WRITE (P01REC,FMT=99994) I
                  IERROR = 2
                  GO TO 380
               ELSE
                  IDIM(K) = LFAC(I)
                  K = K - 1
                  J = J + 1
                  IWK(J) = IWK(J-1)*LFAC(I)
               END IF
            END IF
   80    CONTINUE
         NCELLS = IWK(J)
         IF (NCELLS.GT.MAXT) THEN
            IERROR = 2
            NREC = 2
            WRITE (P01REC,FMT=99989) MAXT, NCELLS
            GO TO 380
         END IF
         ICS = NDIM + 1
         DO 100 I = 1, NCELLS
            ICOUNT(I) = 0
  100    CONTINUE
C
C        Find cell counts
C
         DO 160 I = 1, N
            IF (WTL) THEN
               IF (WT(I).LT.0.0D0) THEN
                  IERROR = 2
                  WRITE (P01REC,FMT=99990) I
                  GO TO 380
               ELSE IF (WT(I).LE.0.0D0) THEN
                  GO TO 140
               END IF
            END IF
            K = 0
            DO 120 J = NFAC, 1, -1
               IF (ISF(J).GT.0) THEN
                  IF (IFAC(I,J).LE.0) THEN
                     IERROR = 2
                     WRITE (P01REC,FMT=99993) I, J
                     GO TO 380
                  ELSE IF (IFAC(I,J).GT.LFAC(J)) THEN
                     IERROR = 2
                     NREC = 2
                     WRITE (P01REC,FMT=99992) I, J, J, LFAC(J)
                     GO TO 380
                  ELSE
                     K = K + 1
                     IWK(K) = IFAC(I,J)
                  END IF
               END IF
  120       CONTINUE
            CALL G04CAU(.FALSE.,.TRUE.,IWK(ICS),NDIM,K,IWK,IFAULT)
            ICOUNT(K) = ICOUNT(K) + 1
  140       CONTINUE
  160    CONTINUE
C
C        Cummulate cell counts
C
         DO 180 I = 2, NCELLS
            ICOUNT(I) = ICOUNT(I) + ICOUNT(I-1)
  180    CONTINUE
         NOBS = ICOUNT(NCELLS)
C
C        Sort observations into cells
C
         DO 240 I = 1, N
            IF (WTL) THEN
               IF (WT(I).LE.0.0D0) GO TO 220
            END IF
            K = 0
            DO 200 J = NFAC, 1, -1
               IF (ISF(J).GT.0) THEN
                  K = K + 1
                  IWK(K) = IFAC(I,J)
               END IF
  200       CONTINUE
            CALL G04CAU(.FALSE.,.TRUE.,IWK(ICS),NDIM,K,IWK,IFAULT)
            WK(ICOUNT(K)) = Y(I)
            IF (WTL) WK(N+ICOUNT(K)) = WT(I)
            ICOUNT(K) = ICOUNT(K) - 1
  220       CONTINUE
  240    CONTINUE
      END IF
C
C       Compute quantiles
C
      EMPTY = .FALSE.
      IF (WTL) THEN
         IRANK = 2*NDIM + 1
         IFAULT = 0
         DO 320 I = 1, NCELLS
            IL = ICOUNT(I)
            IF (I.EQ.NCELLS) THEN
               IU = NOBS
            ELSE
               IU = ICOUNT(I+1)
            END IF
            IF (IU.GT.IL) THEN
               CALL M01DAF(WK,IL+1,IU,'A',IWK(IRANK),IFAULT)
               CALL M01EAF(WK,IL+1,IU,IWK(IRANK),IFAULT)
               CALL M01EAF(WK(N+1),IL+1,IU,IWK(IRANK),IFAULT)
               WSUM = 0.0D0
               DO 260 J = N + IL + 1, N + IU
                  WSUM = WSUM + WK(J)
  260          CONTINUE
               IF (DDN) THEN
                  APT = (PERCNT/100.00D0)*WSUM
                  IF (APT.GE.WSUM) THEN
                     TABLE(I) = WK(IU)
                     GO TO 320
                  END IF
                  SUM = 0.0D0
                  DO 280 J = IL + 1, IU
                     SUM = SUM + WK(N+J)
                     IF (SUM.EQ.APT) THEN
                        TABLE(I) = 0.5D0*(WK(J)+WK(J+1))
                        GO TO 320
                     ELSE IF (SUM.GT.APT) THEN
                        TABLE(I) = WK(J)
                        GO TO 320
                     END IF
  280             CONTINUE
               ELSE
                  CPRMAX = WSUM - 0.5D0*WK(N+IU)
                  APT = (PERCNT/100.00D0)*CPRMAX
                  CPRP = 0.5D0*WK(N+IL+1)
                  IF (APT.GE.CPRMAX) THEN
                     TABLE(I) = WK(IU)
                  ELSE IF (APT.LE.CPRP) THEN
                     TABLE(I) = WK(IL+1)
                  ELSE
                     SUM = WK(N+IL+1)
                     DO 300 J = IL + 2, IU
                        CPR = SUM + 0.5D0*WK(N+J)
                        IF (APT.LE.CPR) THEN
                           F = (APT-CPRP)/(CPR-CPRP)
                           TABLE(I) = (1.0D0-F)*WK(J-1) + F*WK(J)
                           GO TO 320
                        END IF
                        SUM = SUM + WK(N+J)
                        CPRP = CPR
  300                CONTINUE
                  END IF
               END IF
            ELSE
               TABLE(I) = 0.0D0
               EMPTY = .TRUE.
            END IF
  320    CONTINUE
      ELSE
         DO 340 I = 1, NCELLS
            IFAULT = 0
            IL = ICOUNT(I)
            IF (I.EQ.NCELLS) THEN
               IU = NOBS
            ELSE
               IU = ICOUNT(I+1)
            END IF
            IF (IU.GT.IL) THEN
               M = IU - IL
               CALL M01CAF(WK,IL+1,IU,'A',IFAULT)
               IF (DDN) THEN
                  APT = (PERCNT/100.00D0)*DBLE(M)
                  K = INT(APT)
                  IF (K.EQ.0) THEN
                     TABLE(I) = WK(IL+1)
                  ELSE IF (K.EQ.M) THEN
                     TABLE(I) = WK(IL+K)
                  ELSE IF (DBLE(K).EQ.APT) THEN
                     TABLE(I) = 0.5D0*(WK(IL+K)+WK(IL+K+1))
                  ELSE
                     TABLE(I) = WK(IL+K+1)
                  END IF
               ELSE
                  APT = (PERCNT/100.00D0)*(DBLE(M)-0.5D0) + 0.5D0
                  K = INT(APT)
                  IF (K.EQ.0) THEN
                     TABLE(I) = WK(IL+1)
                  ELSE IF (K.EQ.M) THEN
                     TABLE(I) = WK(IL+K)
                  ELSE
                     F = APT - DBLE(K)
                     TABLE(I) = (1.0D0-F)*WK(IL+K) + F*WK(IL+K+1)
                  END IF
               END IF
            ELSE
               TABLE(I) = 0.0D0
               EMPTY = .TRUE.
            END IF
  340    CONTINUE
      END IF
      DO 360 I = 1, NCELLS - 1
         ICOUNT(I) = ICOUNT(I+1) - ICOUNT(I)
  360 CONTINUE
      ICOUNT(NCELLS) = NOBS - ICOUNT(NCELLS)
      IF (EMPTY) THEN
         IERROR = 3
         WRITE (P01REC,FMT=99988)
      END IF
  380 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N .lt. 2: N = ',I16)
99998 FORMAT (' ** On entry, NFAC .lt. 1: NFAC = ',I16)
99997 FORMAT (' ** On entry, LDF .lt. N: LDF = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, TYPE is not valid: TYPE = ',A1)
99995 FORMAT (' ** On entry, PERCNT is not valid: PERCNT = ',D13.5)
99994 FORMAT (' ** On entry, the ',I16,'th element of LFAC .le. 1')
99993 FORMAT (' ** On entry, the ',I16,',',I16,'th element of IFAC .le',
     *       '. 0')
99992 FORMAT (' ** On entry, the ',I16,',',I16,'th element of',/'     ',
     *       '  ','        IFAC .gt. LFAC(',I16,') = ',I16)
99991 FORMAT (' ** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99990 FORMAT (' ** On entry, the ',I16,'th element of WT .lt. 0')
99989 FORMAT (' ** On entry, MAXT too small: MAXT = ',I16,/'          ',
     *       '         minimum value   = ',I16)
99988 FORMAT (' ** Some cells of the table are empty')
      END
