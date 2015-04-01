      SUBROUTINE G03ECF(METHOD,N,D,ILC,IUC,CD,IORD,DORD,IWK,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Performs Hierarchical Cluster Analysis using distance matrix D.
C
C     METHOD indicates which type is performed
C     METHOD = 1 - single link
C     METHOD = 2 - complete link
C     METHOD = 3 - group average
C     METHOD = 4 - centroid
C     METHOD = 5 - median
C     METHOD = 6 - minimum variance
C
C     ILC, IUC and CD give information about the clustering process
C
C     IORD and DORD give information for printing the dendrogram
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G03ECF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, METHOD, N
C     .. Array Arguments ..
      DOUBLE PRECISION  CD(N-1), D(N*(N-1)/2), DORD(N)
      INTEGER           ILC(N-1), IORD(N), IUC(N-1), IWK(2*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DMIN
      INTEGER           I, IERROR, IFAULT, ITH, JTH, NN
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G03ECW
      INTEGER           P01ABF
      EXTERNAL          G03ECW, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G03ECX, G03ECY, G03ECZ
C     .. Executable Statements ..
      IERROR = 0
      IF (METHOD.LT.1 .OR. METHOD.GT.6) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) METHOD
      ELSE IF (N.LT.2) THEN
         IERROR = 1
         WRITE (REC,FMT=99998) N
      ELSE
         NN = N*(N-1)/2
         DO 20 I = 1, NN
            IF (D(I).LT.0.0D0) THEN
               IERROR = 2
            END IF
   20    CONTINUE
         IF (IERROR.EQ.2) THEN
            WRITE (REC,FMT=99997)
         END IF
      END IF
      IF (IERROR.EQ.0) THEN
         DO 40 I = 1, N
            IWK(I) = I
            IWK(N+I) = 1
   40    CONTINUE
         DO 60 I = 1, N - 1
            CALL G03ECX(N,D,IWK,ITH,JTH,DMIN)
            IF (I.NE.N-1) THEN
               CALL G03ECY(METHOD,G03ECW,N,D,IWK,IWK(N+1),ITH,JTH,DMIN)
            END IF
            ILC(I) = JTH
            IUC(I) = ITH
            CD(I) = DMIN
            IF (I.GT.1) THEN
               IF (DMIN.LT.CD(I-1)) THEN
                  IERROR = 3
                  WRITE (REC,FMT=99996)
               END IF
            END IF
   60    CONTINUE
         CALL G03ECZ(N,ILC,IUC,CD,IORD,DORD,IWK,IFAULT)
         IF (IFAULT.NE.0) THEN
            IERROR = 4
         END IF
      END IF
      CONTINUE
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, METHOD.lt.1 .or. METHOD.gt.6: METHOD = ',
     *       I16)
99998 FORMAT (1X,'** On entry, N.lt.2: N = ',I16)
99997 FORMAT (1X,'** On entry, at least one element of D is negative.')
99996 FORMAT (1X,'** Minimum cluster distance not increasing, dendrogr',
     *       'am invalid.')
      END
