      DOUBLE PRECISION FUNCTION G03ECW(ITYPE,DKI,DKJ,DIJ,NI,NJ,NK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes the distance (DKI) from cluster k to cluster
C     formed by merging clusters i and j.
C     ITYPE indicates method:
C           1 - slingle link
C           2 - complete link
C           3 - group average
C           4 - centroid
C           5 - median
C           6 - minimum variance
C
C     Note the values of NI, NJ and NK are asummed gt 0
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DIJ, DKI, DKJ
      INTEGER                          ITYPE, NI, NJ, NK
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE, MAX, MIN
C     .. Executable Statements ..
C
      IF (ITYPE.EQ.1) THEN
         G03ECW = MIN(DKI,DKJ)
      ELSE IF (ITYPE.EQ.2) THEN
         G03ECW = MAX(DKI,DKJ)
      ELSE IF (ITYPE.EQ.3) THEN
         G03ECW = (DBLE(NI)*DKI+DBLE(NJ)*DKJ)/DBLE(NI+NJ)
      ELSE IF (ITYPE.EQ.4) THEN
         G03ECW = (DBLE(NI)*DKI+DBLE(NJ)*DKJ-DBLE(NI)*DBLE(NJ)
     *            *DIJ/DBLE(NI+NJ))/DBLE(NI+NJ)
      ELSE IF (ITYPE.EQ.5) THEN
         G03ECW = 0.5D0*DKI + 0.5D0*DKJ - 0.25D0*DIJ
      ELSE IF (ITYPE.EQ.6) THEN
         G03ECW = (DBLE(NI+NK)*DKI+DBLE(NJ+NK)*DKJ-DBLE(NK)*DIJ)
     *            /DBLE(NI+NJ+NK)
      END IF
      RETURN
      END
