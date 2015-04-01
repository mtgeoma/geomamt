      SUBROUTINE D02MWY(NEQMAX,NY2DIM,METOD1,METOD2,THETA,CONST,TCRIT,
     *                  HMIN,HMAX,H0,MAXSTP,MXHNIL,NORM,RWORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     THETA B INTEGRATOR FOR NAG LIBRARY.
C     THIS IS THE VERSION FOR USE WITH THE REVISED SPRINT ROUTINES.
C***********************************************************************
C**   D02THF - SETUP ROUTINE FOR THE THETA METHOD INTEGRATOR D02MWX   **
C**   VERSION 1  JULY 1986        T.D.                UNIV OF LEEDS   **
C***********************************************************************
C***********************************************************************
C INPUT PARAMETERS
C
C  NEQMAX - INTEGER
C   THE MAXIMUM NUMBER OF EQUATIONS
C
C  NY2DIM - INTEGER
C   THE SECOND DIMENSION OF THE WORK ARRAY YSAVE USED BY D02MWX
C   YSAVE(NEQN, NY2DIM); NY2DIM MUST BE >= 4
C
C  METOD1 - CHARACTER*1
C          = 'N' IF A NEWTON METHOD IS TO BE USED AT THE START
C            OF THE INTEGRATION.
C          = 'F' IF FUNCTIONAL ITERATION IS TO BE USED AT THE
C            START OF THE INTEGRATION.
C
C  METOD2 - CHARACTER*1
C          = 'N' IF NO ATTEMPT TO SWITCH WILL BE MADE - THE CODE
C            WILL STAY WITH THE METHOD DEFINED BY METOD1.
C          = 'S' FOR SWITCHING BETWEEN NEWTON METHOD AND FUNCTIONAL
C            ITERATION TO BE USED TO TRY IMPROVE EFFIENCY.
C
C  THETA -  DOUBLE PRECISION
C           A DOUBLE PRECISION VARIABLE THAT CONTAINS
C           THE VALUE OF THETA TO BE USED BY THE CODE
C           THE RECOMMENDED VALUE OF THETA IS
C    0.55   BUT ANY VALUE BETWEEN 0.51 AND 0.90 MAY
C           BE USED.
C      N.B. IN ORDER TO USE INTEGRATION FORMULAE WITH THETA = 1.0
C           USE THE BACKWARD DIFFERENCE FORMULA OF ORDER 1 AND TO
C           USE THETA = 0.5 USE THE ADAMS FORMULA OF ORDER 2 -SEE
C           THE DOCUMENTATION OF MODULE SPGEAR FOR FURTHER DETAILS.
C
C***********************************************************************
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02MWY')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H0, HMAX, HMIN, TCRIT, THETA
      INTEGER           IFAIL, MAXSTP, MXHNIL, NEQMAX, NY2DIM
      CHARACTER*1       METOD1, METOD2, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  CONST(3), RWORK(50+4*NEQMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION  THETB
      INTEGER           IDEV, IERR, NMETH1, NMETH2
      LOGICAL           REPORT
      CHARACTER*1       METHD1, METHD2
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02MWZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IERR = 0
      CALL X04AAF(0,IDEV)
      IF (NY2DIM.LT.4) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99999) NY2DIM
            CALL X04BAF(IDEV,REC(1))
         END IF
      ELSE IF (THETA.LT.0.51D0 .OR. THETA.GT.1.0D0) THEN
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99998)
            CALL X04BAF(IDEV,REC(1))
         END IF
         THETB = 0.55D0
      ELSE
         THETB = THETA
      END IF
      RWORK(44) = THETB
C
C sfnofc ??
      RWORK(28) = 4.0D0
C
      IF (METOD1.EQ.'N' .OR. METOD1.EQ.'n') THEN
         NMETH1 = 30
      ELSE IF (METOD1.EQ.'F' .OR. METOD1.EQ.'f') THEN
         NMETH1 = 3
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99997) METOD1
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      RWORK(26) = DBLE(NMETH1)
C
      IF (METOD2.EQ.'N' .OR. METOD2.EQ.'n') THEN
         NMETH2 = 40
      ELSE IF (METOD2.EQ.'S' .OR. METOD2.EQ.'s') THEN
         NMETH2 = 4
      ELSE
         IERR = 1
         IF (REPORT) THEN
            WRITE (REC(1),FMT=99996) METOD2
            CALL X04BAF(IDEV,REC(1))
         END IF
      END IF
      RWORK(27) = DBLE(NMETH2)
C...CALL D02MWZ - ROUTINE TO INITIALISE RWORK...
      CALL D02MWZ(NEQMAX,NY2DIM,THETA,CONST,TCRIT,HMIN,HMAX,H0,MAXSTP,
     *            MXHNIL,NORM,RWORK,IERR,REPORT)
C
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,REC)
      RETURN
C
99999 FORMAT (' ** On entry, NY2DIM.lt.4 : NY2DIM =',I16)
99998 FORMAT (' ** On entry, THETA is outside allowed range, THETA =',
     *       ' 0.55 assumed')
99997 FORMAT (' ** On entry, METOD1 is not valid : METOD1 = ',A1)
99996 FORMAT (' ** On entry, METOD2 is not valid : METOD2 = ',A1)
      END
