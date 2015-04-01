      SUBROUTINE D02MVF(NEQMAX,NY2DIM,MAXORD,CONST,TCRIT,HMIN,HMAX,H0,
     *                  MAXSTP,MXHNIL,NORM,RWORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C***********************************************************************
C   THIS ROUTINE INITIALISES THE DASSL IMPLENTATION OF B.D.F
C   METHODS AS USED IN THE D02N SUB CHAPTER
C INPUT PARAMETERS
C
C  NEQMAX - INTEGER
C   THE MAXIMUM NUMBER OF EQUATIONS
C   NY2DIM     THE SECOND DIMENSION OF THE YSAVE ARRAY , IT SHOULD BE
C              SET TO MAXORD + 3, UNLESS MAXORD IS ZERO ,SEE BELOW.
C   MAXORD     THE MAXIMUM ORDER TO BE USED. THIS SHOULD BE POSITIVE AND
C              THIS SHOULD BE LESS THAN 6 . WHEN MAXORD IS SET
C              TO ZERO THE DEFAULT VALUES OF  5
C              WILL BE USED BY THE CODE (NY2DIM MUST BE LARGE ENOUGH
C              FOR THESE DEFAULT VALUES).
C
C
C          N.B. THE SIZE OF THE YSAVE ARRAY PASSED INTO NITE3  MUST BE
C              YSAVE(NYH, MAXORD + 3). WHERE NYH IS THE NUMBER OF
C              DIFFERENTIAL ALGEBRAIC EQUATIONS OR AN UPPER BOUND ON
C              THIS NUMBER .
C----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02MVF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H0, HMAX, HMIN, TCRIT
      INTEGER           IFAIL, MAXORD, MAXSTP, MXHNIL, NEQMAX, NY2DIM
      CHARACTER*1       NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  CONST(3), RWORK(50+4*NEQMAX)
C     .. Local Scalars ..
      INTEGER           I, IDEV, IERR, MXORD1
      LOGICAL           REPORT
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02MVZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IERR = 0
      REPORT = IFAIL .LE. 0
      CALL X04AAF(0,IDEV)
      IF (NEQMAX.LT.1) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99996) NEQMAX
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (MAXORD.LT.0) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) MAXORD
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      IF (MAXORD.EQ.0) THEN
         MXORD1 = 5
      ELSE
         IF (MAXORD.GT.5) THEN
            IERR = 1
            IF (REPORT) THEN
               WRITE (REC(1),FMT=99998)
               CALL X04BAF(IDEV,REC(1))
            END IF
         END IF
         MXORD1 = MIN(MAXORD,5)
      END IF
      IF (NY2DIM.LT.(MXORD1+3)) THEN
         IERR = 1
         I = MXORD1 + 3
         IF (REPORT) THEN
            WRITE (REC,FMT=99997) NY2DIM, I
            CALL X04BAF(IDEV,REC(1))
            CALL X04BAF(IDEV,REC(2))
         END IF
         GO TO 20
      END IF
C
C...CALL D02MVZ - ROUTINE TO INITIALISE RWORK...
      CALL D02MVZ(NEQMAX,NY2DIM,MXORD1,CONST,TCRIT,HMIN,HMAX,H0,MAXSTP,
     *            MXHNIL,NORM,RWORK,IERR,REPORT)
C
C...ERROR RETURN...
C
   20 CONTINUE
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
C------END-OF-ROUTINE-D02MVF-------------------------------------------
C
99999 FORMAT (' ** On entry, MAXORD.lt.0 : MAXORD =',I16)
99998 FORMAT (' ** On entry, MAXORD.gt.5, MAXORD = 5 assumed')
99997 FORMAT (' ** On entry, NY2DIM.lt.MAXORD+3 : NY2DIM =',I16,/'    ',
     *       'MAXORD+3 =',I16)
99996 FORMAT (' ** On entry, NEQMAX.lt.1 : NEQMAX =',I16)
      END
