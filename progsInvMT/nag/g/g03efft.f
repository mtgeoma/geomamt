      SUBROUTINE G03EFF(WEIGHT,N,M,X,LDX,ISX,NVAR,K,CMEANS,LDC,WT,INC,
     *                  NIC,CSS,CSW,MAXIT,IWK,WK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     K-Means clustering
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      CHARACTER*6       SRNAME
      PARAMETER         (ZERO=0.0D0,SRNAME='G03EFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, K, LDC, LDX, M, MAXIT, N, NVAR
      CHARACTER         WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  CMEANS(LDC,NVAR), CSS(K), CSW(K), WK(2*K+N),
     *                  WT(*), X(LDX,M)
      INTEGER           INC(N), ISX(M), IWK(N+3*K), NIC(K)
C     .. Local Scalars ..
      INTEGER           COUNT, I, IERR, IFAIL1, NREC, WTSUM
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G03EFW, G03EFZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      COUNT = 0
      IERR = 1
      NREC = 1
      IF (WEIGHT.NE.'W' .AND. WEIGHT.NE.'w' .AND. WEIGHT.NE.'U' .AND.
     *    WEIGHT.NE.'u') THEN
         WRITE (REC,FMT=99999) WEIGHT
      ELSE IF (N.LT.2) THEN
         WRITE (REC,FMT=99998) N
      ELSE IF (NVAR.LT.1) THEN
         WRITE (REC,FMT=99997) NVAR
      ELSE IF (M.LT.NVAR) THEN
         WRITE (REC,FMT=99996) M, NVAR
      ELSE IF (K.LE.1) THEN
         WRITE (REC,FMT=99995) K
      ELSE IF (LDX.LT.N) THEN
         WRITE (REC,FMT=99994) LDX, N
      ELSE IF (LDC.LT.K) THEN
         WRITE (REC,FMT=99993) LDC, K
      ELSE IF (MAXIT.LE.0) THEN
         WRITE (REC,FMT=99992) MAXIT
      ELSE
         IERR = 0
         IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            WTSUM = 0
            DO 20 I = 1, N
               IF (WT(I).GT.ZERO) THEN
                  WTSUM = WTSUM + 1
               ELSE IF (WT(I).LT.ZERO) THEN
                  IERR = 2
                  WRITE (REC,FMT=99991) I
               END IF
   20       CONTINUE
            IF (WTSUM.LT.2) THEN
               IERR = 2
               WRITE (REC,FMT=99990)
            END IF
         END IF
         IF (IERR.EQ.0) THEN
            DO 40 I = 1, M
               IF (ISX(I).GT.0) THEN
                  COUNT = COUNT + 1
               END IF
   40       CONTINUE
            IF (COUNT.NE.NVAR) THEN
               IERR = 3
               WRITE (REC,FMT=99989)
            END IF
         END IF
         IF (IERR.EQ.0) THEN
            IFAIL1 = 0
            IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
               CALL G03EFW(N,M,X,LDX,ISX,NVAR,K,CMEANS,LDC,WT,INC,CSS,
     *                     CSW,MAXIT,IWK,IWK(N+1),WK,IWK(N+K+1),
     *                     IWK(2*K+N+1),IFAIL1)
               DO 60 I = 1, K
                  NIC(I) = 0
   60          CONTINUE
               DO 80 I = 1, N
                  IF (WT(I).GT.ZERO) NIC(INC(I)) = NIC(INC(I)) + 1
   80          CONTINUE
            ELSE
               CALL G03EFZ(N,M,X,LDX,ISX,NVAR,K,CMEANS,LDC,INC,NIC,CSS,
     *                     MAXIT,IWK,WK,WK(K+1),IWK(N+1),WK(2*K+1),
     *                     IWK(K+N+1),IWK(2*K+N+1),IFAIL1)
               DO 100 I = 1, K
                  CSW(I) = DBLE(NIC(I))
  100          CONTINUE
            END IF
            IF (IFAIL1.EQ.1) THEN
               IERR = 4
               WRITE (REC,FMT=99988)
            ELSE IF (IFAIL1.EQ.2) THEN
               IERR = 5
               NREC = 2
               WRITE (REC,FMT=99987) MAXIT
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, WEIGHT is not valid: WEIGHT = ',A1)
99998 FORMAT (' ** On entry, N.lt.2:  N =',I16)
99997 FORMAT (' ** On entry, NVAR.lt.1:  NVAR =',I16)
99996 FORMAT (' ** On entry, M.lt.NVAR:  M =',I16,'  NVAR =',I16)
99995 FORMAT (' ** On entry, K.le.1:  K =',I16)
99994 FORMAT (' ** On entry, LDX.lt.N:  LDX =',I16,'  N =',I16)
99993 FORMAT (' ** On entry, LDC.lt.K:  LDC =',I16,'  K =',I16)
99992 FORMAT (' ** On entry, MAXIT.le.0:  MAXIT =',I16)
99991 FORMAT (' ** On entry, the value of WT(',I16,').lt.0.0')
99990 FORMAT (' ** On entry, WT has less than two positive values')
99989 FORMAT (' ** On entry, number of positive values in ISX does not',
     *       ' equal NVAR')
99988 FORMAT (' ** At least one cluster is empty after the initial ass',
     *       'ignment')
99987 FORMAT (' ** Convergence has not been achieved within the maximu',
     *       'm number of ',/'    iterations    MAXIT =',I16)
      END
