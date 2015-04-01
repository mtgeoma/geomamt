      SUBROUTINE G10CAF(ITYPE,N,Y,SMOOTH,ROUGH,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     EDA SMOOTHERS
C
C     ITYPE  specifies the smoother to be used
C            ITYPE=0 specifies 4253H, twice
C            ITYPE=1 specifies 3RSSH, twice
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G10CAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, ITYPE, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ROUGH(N), SMOOTH(N), Y(N)
C     .. Local Scalars ..
      INTEGER           I, IERROR, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, G10CAY, G10CAZ
C     .. Executable Statements ..
C
      IERROR = 0
      IF (ITYPE.LT.0 .OR. ITYPE.GT.1) THEN
         IERROR = 1
         NREC = 1
         WRITE (REC,FMT=99999) ITYPE
      ELSE IF (N.LE.6) THEN
         IERROR = 2
         NREC = 1
         WRITE (REC,FMT=99998) N
      ELSE
         CALL DCOPY(N,Y,1,SMOOTH,1)
         IF (ITYPE.EQ.0) CALL G10CAZ(SMOOTH,N)
         IF (ITYPE.EQ.1) CALL G10CAY(SMOOTH,N)
C
C        compute rough from first smoothing
C
         DO 20 I = 1, N
            ROUGH(I) = Y(I) - SMOOTH(I)
   20    CONTINUE
C
C        Re-rough smoothers ("twicing")
C
         IF (ITYPE.EQ.0) CALL G10CAZ(ROUGH,N)
         IF (ITYPE.EQ.1) CALL G10CAY(ROUGH,N)
         DO 40 I = 1, N
            SMOOTH(I) = SMOOTH(I) + ROUGH(I)
            ROUGH(I) = Y(I) - SMOOTH(I)
   40    CONTINUE
      END IF
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** On entry, ITYPE.lt.0  or  ITYPE.gt.1: ITYPE =',I16)
99998 FORMAT (' ** On entry, N.le.6: N = ',I16)
      END
