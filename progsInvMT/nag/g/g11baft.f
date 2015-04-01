      SUBROUTINE G11BAF(STAT,UPDATE,WEIGHT,N,NFAC,ISF,LFAC,IFAC,LDF,Y,
     *                  WT,TABLE,MAXT,NCELLS,NDIM,IDIM,ICOUNT,AUXT,IWK,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes table from data and factors using selected statistics.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G11BAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDF, MAXT, N, NCELLS, NDIM, NFAC
      CHARACTER         STAT, UPDATE, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  AUXT(*), TABLE(MAXT), WT(*), Y(N)
      INTEGER           ICOUNT(MAXT), IDIM(NFAC), IFAC(LDF,NFAC),
     *                  ISF(NFAC), IWK(2*NFAC), LFAC(NFAC)
C     .. Local Scalars ..
      DOUBLE PRECISION  DEV, SW
      INTEGER           I, ICS, IERROR, IFAULT, IND, J, K, NREC
      LOGICAL           NZV, SINGLE, WTL
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G04CAU
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
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
      ELSE IF (UPDATE.NE.'U' .AND. UPDATE.NE.'u' .AND. UPDATE.NE.
     *         'I' .AND. UPDATE.NE.'i') THEN
         WRITE (P01REC,FMT=99995) UPDATE
      ELSE IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.
     *         'U' .AND. WEIGHT.NE.'u' .AND. WEIGHT.NE.'V' .AND.
     *         WEIGHT.NE.'v') THEN
         WRITE (P01REC,FMT=99991) WEIGHT
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 220
      IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w' .OR. WEIGHT.EQ.'V' .OR.
     *    WEIGHT.EQ.'v') THEN
         WTL = .TRUE.
         IF (WEIGHT.EQ.'V' .OR. WEIGHT.EQ.'v') THEN
            NZV = .TRUE.
         ELSE
            NZV = .FALSE.
         END IF
      ELSE
         WTL = .FALSE.
      END IF
C
C     Compute number of cells in the table.
C
      NDIM = 0
      DO 20 I = 1, NFAC
         IF (ISF(I).GT.0) NDIM = NDIM + 1
   20 CONTINUE
      IF (NDIM.EQ.0) THEN
         NDIM = 1
         IWK(2) = 1
         IDIM(1) = 1
         SINGLE = .TRUE.
      ELSE
         SINGLE = .FALSE.
         J = NDIM
         IWK(J) = 1
         K = NDIM
         DO 40 I = NFAC, 1, -1
            IF (ISF(I).GT.0) THEN
               IF (LFAC(I).LE.1) THEN
                  WRITE (P01REC,FMT=99994) I
                  IERROR = 2
                  GO TO 220
               ELSE
                  IDIM(K) = LFAC(I)
                  K = K - 1
                  J = J + 1
                  IWK(J) = IWK(J-1)*LFAC(I)
               END IF
            END IF
   40    CONTINUE
      END IF
      IF (MAXT.LT.IWK(2*NDIM)) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99985) MAXT, IWK(2*NDIM)
         GO TO 220
      END IF
      ICS = NDIM + 1
      IF (STAT.EQ.'N' .OR. STAT.EQ.'n') THEN
         IND = 1
      ELSE IF (STAT.EQ.'T' .OR. STAT.EQ.'t') THEN
         IND = 2
      ELSE IF (STAT.EQ.'A' .OR. STAT.EQ.'a') THEN
         IND = 3
      ELSE IF (STAT.EQ.'V' .OR. STAT.EQ.'v') THEN
         IND = 5
      ELSE IF (STAT.EQ.'L' .OR. STAT.EQ.'l') THEN
         IND = 6
      ELSE IF (STAT.EQ.'S' .OR. STAT.EQ.'s') THEN
         IND = 7
      ELSE
         IERROR = 1
         WRITE (P01REC,FMT=99996) STAT
         GO TO 220
      END IF
C
C     Initialize TABLE if UPDATE = I
C
      IF (UPDATE.EQ.'I' .OR. UPDATE.EQ.'i') THEN
         NCELLS = IWK(2*NDIM)
         DO 60 I = 1, NCELLS
            ICOUNT(I) = 0
            TABLE(I) = 0.0D0
   60    CONTINUE
         IF (IND.EQ.3) THEN
            DO 80 I = 1, NCELLS
               AUXT(I) = 0.0D0
   80       CONTINUE
         ELSE IF (IND.EQ.5) THEN
            DO 100 I = 1, 2*NCELLS
               AUXT(I) = 0.0D0
  100       CONTINUE
         END IF
      ELSE
C
C        Check re-input values if UPDATE = U
C
         IF (NCELLS.NE.IWK(2*NDIM)) THEN
            IERROR = 4
            NREC = 2
            WRITE (P01REC,FMT=99988)
            GO TO 220
         END IF
         IF (IND.EQ.3 .OR. IND.EQ.5) THEN
            DO 120 I = 1, NCELLS
               IF (TABLE(I).GT.0.0D0 .AND. AUXT(I).LE.0.0D0) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99987)
                  GO TO 220
               END IF
  120       CONTINUE
         END IF
         IF (IND.EQ.5) THEN
C
C           Compute sum of squares from variances
C
            DO 140 I = 1, NCELLS
               IF (NZV .AND. TABLE(I).GT.0.0D0 .AND. ICOUNT(I).LE.1)
     *             THEN
                  IERROR = 4
               ELSE IF (NZV .AND. TABLE(I).LE.0.0D0 .AND. ICOUNT(I)
     *                  .GT.1) THEN
                  IERROR = 4
               ELSE IF ( .NOT. NZV .AND. TABLE(I).GT.0.0D0 .AND. AUXT(I)
     *                  .LE.1.0D0) THEN
                  IERROR = 4
               ELSE IF ( .NOT. NZV .AND. TABLE(I).LE.0.0D0 .AND. AUXT(I)
     *                  .GT.1.0D0) THEN
                  IERROR = 4
               END IF
               IF (IERROR.EQ.4) THEN
                  WRITE (P01REC,FMT=99986)
                  GO TO 220
               ELSE IF (TABLE(I).LT.0.0D0) THEN
                  TABLE(I) = -TABLE(I)
               ELSE IF (TABLE(I).GT.0.0D0) THEN
                  IF (NZV) THEN
                     TABLE(I) = TABLE(I)*(DBLE(ICOUNT(I))-1.0D0)
                  ELSE
                     TABLE(I) = TABLE(I)*(AUXT(I)-1.0D0)
                  END IF
               END IF
  140       CONTINUE
         END IF
      END IF
C
C        Loop through data
C
      DO 180 I = 1, N
         IF (WTL) THEN
            IF (WT(I).LT.0.0D0) THEN
               IERROR = 2
               WRITE (P01REC,FMT=99990) I
               GO TO 220
            ELSE IF (WT(I).EQ.0.0D0) THEN
               GO TO 180
            END IF
         END IF
         IF (SINGLE) THEN
            K = 1
         ELSE
            K = 0
            DO 160 J = NFAC, 1, -1
               IF (ISF(J).GT.0) THEN
                  IF (IFAC(I,J).LE.0) THEN
                     IERROR = 2
                     WRITE (P01REC,FMT=99993) I, J
                     GO TO 220
                  ELSE IF (IFAC(I,J).GT.LFAC(J)) THEN
                     IERROR = 2
                     NREC = 2
                     WRITE (P01REC,FMT=99992) I, J, LFAC(J)
                     GO TO 220
                  ELSE
                     K = K + 1
                     IWK(K) = IFAC(I,J)
                  END IF
               END IF
  160       CONTINUE
            CALL G04CAU(.FALSE.,.TRUE.,IWK(ICS),NDIM,K,IWK,IFAULT)
         END IF
         ICOUNT(K) = ICOUNT(K) + 1
         IF (IND.EQ.1) THEN
            TABLE(K) = TABLE(K) + 1
         ELSE IF (IND.EQ.2) THEN
            IF (WTL) THEN
               TABLE(K) = TABLE(K) + Y(I)*WT(I)
            ELSE
               TABLE(K) = TABLE(K) + Y(I)
            END IF
         ELSE IF (IND.EQ.3) THEN
            IF (WTL) THEN
               AUXT(K) = AUXT(K) + WT(I)
               TABLE(K) = TABLE(K) + WT(I)*(Y(I)-TABLE(K))/AUXT(K)
            ELSE
               AUXT(K) = AUXT(K) + 1.0D0
               TABLE(K) = TABLE(K) + (Y(I)-TABLE(K))/AUXT(K)
            END IF
         ELSE IF (IND.EQ.5) THEN
            IF (WTL) THEN
               SW = AUXT(K)
               DEV = Y(I) - AUXT(NCELLS+K)
               AUXT(K) = AUXT(K) + WT(I)
               AUXT(NCELLS+K) = AUXT(NCELLS+K) + WT(I)*DEV/AUXT(K)
               TABLE(K) = TABLE(K) + WT(I)*SW*DEV*DEV/AUXT(K)
            ELSE
               SW = AUXT(K)
               DEV = Y(I) - AUXT(NCELLS+K)
               AUXT(K) = AUXT(K) + 1.0D0
               AUXT(NCELLS+K) = AUXT(NCELLS+K) + DEV/AUXT(K)
               TABLE(K) = TABLE(K) + SW*DEV*DEV/AUXT(K)
            END IF
         ELSE IF (IND.EQ.6) THEN
            IF (ICOUNT(K).EQ.1) THEN
               TABLE(K) = Y(I)
            ELSE
               IF (TABLE(K).LT.Y(I)) TABLE(K) = Y(I)
            END IF
         ELSE IF (IND.EQ.7) THEN
            IF (ICOUNT(K).EQ.1) THEN
               TABLE(K) = Y(I)
            ELSE
               IF (TABLE(K).GT.Y(I)) TABLE(K) = Y(I)
            END IF
         END IF
  180 CONTINUE
C
C       Compute variance from sum of squares
C
      IF (IND.EQ.5) THEN
         DO 200 I = 1, NCELLS
            IF (NZV) THEN
               IF (ICOUNT(I).LE.1) THEN
                  IERROR = 3
                  WRITE (P01REC,FMT=99989)
                  TABLE(I) = -TABLE(I)
               ELSE
                  TABLE(I) = TABLE(I)/(DBLE(ICOUNT(I)-1))
               END IF
            ELSE
               IF (AUXT(I).LE.1.0D0) THEN
                  IERROR = 3
                  WRITE (P01REC,FMT=99989)
                  TABLE(I) = -TABLE(I)
               ELSE
                  TABLE(I) = TABLE(I)/(AUXT(I)-1.0D0)
               END IF
            END IF
  200    CONTINUE
      END IF
  220 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, N .lt. 2: N = ',I16)
99998 FORMAT (' ** On entry, NFAC .lt. 1: NFAC = ',I16)
99997 FORMAT (' ** On entry, LDF .lt. N: LDF = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, STAT is not valid: STAT = ',A1)
99995 FORMAT (' ** On entry, UPDATE is not valid: UPDATE = ',A1)
99994 FORMAT (' ** On entry, the ',I16,'th element of LFAC .le. 1')
99993 FORMAT (' ** On entry, the ',I16,',',I16,
     *       'th value of IFAC .le. 0')
99992 FORMAT (' ** On entry, the ',I16,',',I16,'th element of',/'     ',
     *       '  ','        IFAC .gt. LFAC(',I16,') = ',I16)
99991 FORMAT (' ** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99990 FORMAT (' ** On entry, the ',I16,'th element of WT .lt. 0')
99989 FORMAT (' ** The variance divisor .le. 0.0')
99988 FORMAT (' ** One entry with UPDATE = U,',/'    the value of NCEL',
     *       'LS is incompatible with the values of IWK')
99987 FORMAT (' ** One entry with UPDATE = U, AUXT has been changed')
99986 FORMAT (' ** One entry with UPDATE = U, AUXT or TABLE have been ',
     *       'changed')
99985 FORMAT (' ** On entry, MAXT is too small: MAXT = ',I16,' min  = ',
     *       I16)
      END
