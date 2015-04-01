      DOUBLE PRECISION FUNCTION D02ZAF(N,V,W,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 17 REVISED. IER-1547 (JUN 1995).
C
C     VP CHANGED AUG 15th 1994
C
C     OLD NAME VNORM
C
C-----------------------------------------------------------------------
C VP
C INORM = 1         AVERAGE L2 NORM
C       = 2         MAXIMUM NORM
C       = 3         AVERAGE L1 NORM
C  THE CHOICE OF NORM IS DEFINED BY THE USER ON ENTRY TO SPRINT AND
C  PASSED TO THIS FUNCTION BY THE COMMON BLOCK /HD02NM/.
C  NB ONLY D03PFF/PLF/PSF HAVE THE L1 NORM OPTION (INSTEAD OF THE MAX NO
C
C  ALL NORMS ARE WEIGHTED VECTOR NORMS OF
C  THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS
C  CONTAINED IN THE ARRAY W OF LENGTH N.
C  E.G. FOR THE AVERAGED L2 NORM
C       D02ZAF = SQRT( (1/N) * SUM( V(I)*W(I) )**2 )
C      THE L2 NORM CAN BE OBTAINED BY APPROPRIATE SCALING OF THE WEIGHTS
C  FOR THE MAX NORM
C       D02ZAF = MAX '' V(I)*W(I) ''
C                I
C      THE MODIFIED MAX NORM CAN BE OBTAINED , AS FOR THE L2 NORM, BY
C      APPROPRIATE SCALING OF THE WEIGHTS.
C  FOR THE AVERAGED L1 NORM
C       D02ZAF = (1/N) * SUM( V(I)*W(I) )
C      THE L1 NORM CAN BE OBTAINED BY APPROPRIATE SCALING OF THE WEIGHTS
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='D02ZAF')
C     .. Scalar Arguments ..
      INTEGER                          IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION                 V(N), W(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION                 BIG, QTYMIN, ROOTN
      INTEGER                          INORM
C     .. Local Scalars ..
      DOUBLE PRECISION                 BIGTST, SUM, TERM, TERMAX, VMAX
      INTEGER                          I, IERR, LSMALL
C     .. Local Arrays ..
      CHARACTER                        P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX, SQRT
C     .. Common blocks ..
      COMMON                           /HD02NM/INORM
      COMMON                           /JD02NM/BIG, QTYMIN, ROOTN
C     .. Save statement ..
      SAVE                             /HD02NM/, /JD02NM/
C     .. Executable Statements ..
      IERR = 0
      BIGTST = SQRT(0.5D0*BIG)
      VMAX = 0.0D0
      DO 20 I = 1, N
         VMAX = MAX(VMAX,ABS(V(I)*W(I)))
   20 CONTINUE
      IF (VMAX.GE.BIGTST) THEN
         D02ZAF = BIGTST
         IERR = 1
         GO TO 60
C VP
      ELSE IF (VMAX.EQ.0.0D0 .OR. INORM.EQ.2) THEN
         D02ZAF = VMAX
         GO TO 60
      END IF
      LSMALL = 0
      SUM = 0.0D0
      TERMAX = 0.0D0
C205
C205  CAN'T SEE HOW TO VECTORISE THE NEXT LOOP
C205
      DO 40 I = 1, N
         TERM = ABS(V(I)/VMAX*W(I))
         IF (TERM.LE.QTYMIN) THEN
            LSMALL = 1
            TERMAX = MAX(TERMAX,TERM)
         ELSE
            IF (INORM.EQ.1) THEN
               SUM = SUM + TERM**2
            ELSE IF (INORM.EQ.3) THEN
               SUM = SUM + TERM
            END IF
         END IF
   40 CONTINUE
      IF (LSMALL.EQ.1 .AND. SUM.EQ.0.0D0) THEN
         D02ZAF = TERMAX*VMAX/ROOTN
      ELSE
         IF (SUM.GT.0.0D0) THEN
            IF (INORM.EQ.1) THEN
               D02ZAF = SQRT(SUM)/ROOTN*VMAX
            ELSE IF (INORM.EQ.3) THEN
               D02ZAF = SUM/N*VMAX
            END IF
         ELSE
            D02ZAF = 0.0D0
         END IF
      END IF
   60 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
