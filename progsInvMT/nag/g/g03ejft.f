      SUBROUTINE G03EJF(N,CD,IORD,DORD,K,DLEVEL,IC,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes an indicator of group allocation from results of G03ECF
C
C     Either the number of clusters required, K, or the distance at
C     which clusters are to taken, DLEVEL, can be entered.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03EJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DLEVEL
      INTEGER           IFAIL, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CD(N-1), DORD(N)
      INTEGER           IC(N), IORD(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DI, DMAX
      INTEGER           I, IERROR, IPT, M
      LOGICAL           TIE
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IERROR = 0
      IF (N.LT.2) THEN
         IERROR = 1
         WRITE (REC,FMT=99990) N
      ELSE IF (K.GT.N) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) K, N
      ELSE IF (K.LE.0 .AND. DLEVEL.LE.0.0D0) THEN
         IERROR = 1
         WRITE (REC,FMT=99998)
      END IF
      IF (IERROR.EQ.0) THEN
         DI = CD(1)
         DMAX = DORD(1)
         DO 20 I = 2, N - 1
            IF (CD(I).LT.DI) THEN
               IERROR = 2
               WRITE (REC,FMT=99992)
               GO TO 220
            ELSE
               DI = CD(I)
            END IF
            IF (DORD(I).GT.DMAX) DMAX = DORD(I)
   20    CONTINUE
         IF (CD(N-1).NE.DORD(N) .OR. DMAX.NE.DORD(N)) THEN
            IERROR = 2
            WRITE (REC,FMT=99991)
            GO TO 220
         END IF
         IF (K.GT.0) THEN
            IF (K.EQ.1) THEN
               IERROR = 3
               WRITE (REC,FMT=99994)
               DO 40 I = 1, N
                  IC(I) = 1
   40          CONTINUE
            ELSE IF (K.EQ.N) THEN
               IERROR = 3
               WRITE (REC,FMT=99993)
               DO 60 I = 1, N
                  IC(I) = I
   60          CONTINUE
            ELSE
C
C              Find distance from number of clusters
C
               IPT = N - K
               DLEVEL = CD(IPT)
C
C              Check for ties
C
               TIE = .FALSE.
               M = IPT + 1
               DO 80 I = M, N - 1
                  IF (CD(I).EQ.DLEVEL) THEN
                     K = K - 1
                     IPT = IPT + 1
                     TIE = .TRUE.
                  ELSE
                     GO TO 100
                  END IF
   80          CONTINUE
  100          CONTINUE
               IF (TIE) THEN
                  IERROR = 4
                  WRITE (REC,FMT=99997)
               END IF
            END IF
         ELSE
C
C           Find number of clusters from distance
C
            DO 120 I = 1, N - 1
               IF (CD(I).GT.DLEVEL) THEN
                  IPT = I - 1
                  K = N - IPT
                  GO TO 140
               END IF
  120       CONTINUE
            IPT = N - 1
  140       CONTINUE
            IF (IPT.EQ.0) THEN
               IERROR = 3
               WRITE (REC,FMT=99996) DLEVEL
               DO 160 I = 1, N
                  IC(I) = I
  160          CONTINUE
            ELSE IF (IPT.EQ.N-1) THEN
               IERROR = 3
               WRITE (REC,FMT=99995) DLEVEL
               DO 180 I = 1, N
                  IC(I) = 1
  180          CONTINUE
            END IF
         END IF
         IF (IERROR.NE.3) THEN
            M = 1
            IC(1) = 1
            DO 200 I = 2, N
               IF (DORD(I-1).LE.DLEVEL) THEN
                  IC(IORD(I)) = M
               ELSE
                  M = M + 1
                  IC(IORD(I)) = M
               END IF
  200       CONTINUE
         END IF
      END IF
C
  220 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (' ** On entry, K.gt.N: K = ',I16,' N = ',I16)
99998 FORMAT (' ** On entry, K.lt.0 and DLEVEL.le.0.0.')
99997 FORMAT (' ** More than one cluster is formed at one distance.')
99996 FORMAT (' ** No clustering takes place below DLEVEL, DLEVEL = ',
     *       D13.5)
99995 FORMAT (' ** All data merged into one cluster at DLEVEL, DLEVEL ',
     *       '= ',D13.5)
99994 FORMAT (' ** All data is merged when K = 1.')
99993 FORMAT (' ** No clustering is performed when K = N.')
99992 FORMAT (' ** On entry the values of CD are not in increasing ord',
     *       'er')
99991 FORMAT (' ** On entry the values of DORD and CD are not compatib',
     *       'le')
99990 FORMAT (' ** On entry, N.lt.2: N = ',I16)
      END
