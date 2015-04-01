      SUBROUTINE G11BCF(STAT,TABLE,NCELLS,NDIM,IDIM,ISDIM,STABLE,MAXST,
     *                  MCELLS,MDIM,MLEVEL,AUXT,IWK,WK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes table statistics
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G11BCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, MAXST, MCELLS, MDIM, NCELLS, NDIM
      CHARACTER         STAT
C     .. Array Arguments ..
      DOUBLE PRECISION  AUXT(*), STABLE(MAXST), TABLE(NCELLS),
     *                  WK(NCELLS)
      INTEGER           IDIM(NDIM), ISDIM(NDIM), IWK(3*NDIM),
     *                  MLEVEL(NDIM)
C     .. Local Scalars ..
      DOUBLE PRECISION  DEV, RN1, SW
      INTEGER           I, ICELLS, IERROR, IFAULT, II, IND, J, K, MCS,
     *                  NCS, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G04CAT, G04CAU, M01CAF
C     .. Intrinsic Functions ..
      INTRINSIC         NINT
C     .. Executable Statements ..
C
C     Carry out initial error checks
C
      NREC = 1
      IERROR = 1
      IF (NDIM.LT.2) THEN
         WRITE (P01REC,FMT=99999) NDIM
      ELSE
         IERROR = 0
      END IF
      IF (IERROR.EQ.1) GO TO 220
C
C     Check size of input table
C
      K = NDIM
      IWK(K) = 1
      MDIM = 0
      DO 20 I = NDIM, 1, -1
         IF (IDIM(I).LE.1) THEN
            WRITE (P01REC,FMT=99997) I
            IERROR = 2
            GO TO 220
         ELSE
            K = K + 1
            IWK(K) = IWK(K-1)*IDIM(I)
            IF (ISDIM(I).GT.0) MDIM = MDIM + 1
         END IF
   20 CONTINUE
      IF (IWK(K).NE.NCELLS) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99996)
         GO TO 220
      ELSE IF (MDIM.LT.1) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99995)
         GO TO 220
      ELSE IF (MDIM.EQ.NDIM) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99994)
         GO TO 220
      ELSE IF (MDIM.GT.MAXST) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99993) MAXST, MDIM
         GO TO 220
      END IF
C
C     Compute number of cells in reduced table.
C
      J = MDIM
      K = 2*NDIM + 1
      IWK(K) = 1
      DO 40 I = NDIM, 1, -1
         IF (ISDIM(I).GT.0) THEN
            MLEVEL(J) = IDIM(I)
            J = J - 1
            K = K + 1
            IWK(K) = IWK(K-1)*IDIM(I)
         END IF
   40 CONTINUE
      MCELLS = IWK(K)
      NCS = NDIM + 1
      MCS = 2*NDIM + 2
      IF (STAT.EQ.'T' .OR. STAT.EQ.'t') THEN
         IND = 2
      ELSE IF (STAT.EQ.'A' .OR. STAT.EQ.'a') THEN
         IND = 3
      ELSE IF (STAT.EQ.'M' .OR. STAT.EQ.'m') THEN
         IND = 4
         ICELLS = NCELLS/MCELLS
      ELSE IF (STAT.EQ.'V' .OR. STAT.EQ.'v') THEN
         IND = 5
         RN1 = NCELLS/MCELLS - 1
         RN1 = 1.0D0/RN1
      ELSE IF (STAT.EQ.'L' .OR. STAT.EQ.'l') THEN
         IND = 6
      ELSE IF (STAT.EQ.'S' .OR. STAT.EQ.'s') THEN
         IND = 7
      ELSE
         IERROR = 1
         WRITE (P01REC,FMT=99998) STAT
         GO TO 220
      END IF
      DO 60 I = 1, MCELLS
         STABLE(I) = 0.0D0
   60 CONTINUE
      IF (IND.NE.2) THEN
         DO 80 I = 1, MCELLS
            WK(I) = 0.0D0
   80    CONTINUE
      END IF
      IF (IND.EQ.5) THEN
         DO 100 I = 1, MCELLS
            AUXT(I) = 0.0D0
  100    CONTINUE
      END IF
      DO 140 I = 1, NCELLS
         II = I
         CALL G04CAU(.TRUE.,.FALSE.,IWK(NCS),NDIM,II,IWK,IFAULT)
         K = 1
         DO 120 J = 1, NDIM
            IF (ISDIM(J).GT.0) THEN
               IWK(K) = IWK(J)
               K = K + 1
            END IF
  120    CONTINUE
         CALL G04CAT(IWK,MDIM)
         CALL G04CAU(.FALSE.,.TRUE.,IWK(MCS),MDIM,K,IWK,IFAULT)
         IF (IND.EQ.2) THEN
            STABLE(K) = STABLE(K) + TABLE(I)
         ELSE IF (IND.EQ.3) THEN
            WK(K) = WK(K) + 1
            STABLE(K) = STABLE(K) + (TABLE(I)-STABLE(K))/WK(K)
         ELSE IF (IND.EQ.4) THEN
            STABLE(K) = STABLE(K) + 1.0D0
            WK((K-1)*ICELLS+NINT(STABLE(K))) = TABLE(I)
         ELSE IF (IND.EQ.5) THEN
            SW = WK(K)
            DEV = TABLE(I) - AUXT(K)
            WK(K) = WK(K) + 1.0D0
            AUXT(K) = AUXT(K) + DEV/WK(K)
            STABLE(K) = STABLE(K) + SW*DEV*DEV/WK(K)
         ELSE IF (IND.EQ.6) THEN
            IF (WK(K).EQ.0.0D0) THEN
               STABLE(K) = TABLE(I)
               WK(K) = 1.0D0
            ELSE
               IF (STABLE(K).LT.TABLE(I)) STABLE(K) = TABLE(I)
            END IF
         ELSE IF (IND.EQ.7) THEN
            IF (WK(K).EQ.0.0D0) THEN
               STABLE(K) = TABLE(I)
               WK(K) = 1.0D0
            ELSE
               IF (STABLE(K).GT.TABLE(I)) STABLE(K) = TABLE(I)
            END IF
         END IF
  140 CONTINUE
      IF (IND.EQ.4) THEN
         K = (ICELLS+1)/2
         IF (K*2.NE.ICELLS) THEN
            DO 160 I = 1, MCELLS
               IFAULT = 0
               CALL M01CAF(WK,(I-1)*ICELLS+1,I*ICELLS,'A',IFAULT)
               STABLE(I) = WK((I-1)*ICELLS+K)
  160       CONTINUE
         ELSE
            DO 180 I = 1, MCELLS
               IFAULT = 0
               CALL M01CAF(WK,(I-1)*ICELLS+1,I*ICELLS,'A',IFAULT)
               STABLE(I) = (WK((I-1)*ICELLS+K)+WK((I-1)*ICELLS+K+1))
     *                     /2.0D0
  180       CONTINUE
         END IF
      END IF
      IF (IND.EQ.5) THEN
         DO 200 I = 1, MCELLS
            STABLE(I) = STABLE(I)*RN1
  200    CONTINUE
      END IF
  220 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, NDIM .lt. 2: NDIM = ',I16)
99998 FORMAT (' ** On entry, STAT is not valid: STAT = ',A1)
99997 FORMAT (' ** On entry, the ',I16,'th element of IDIM .le. 1')
99996 FORMAT (' ** On entry, NCELLS is incompatible with IDIM')
99995 FORMAT (' ** On entry, no elements of ISDIM .gt. 0')
99994 FORMAT (' ** On entry, all elements of ISDIM .gt. 0')
99993 FORMAT (' ** On entry, MAXST (=',I16,') is too small, ','min val',
     *       'ue=',I16)
      END
