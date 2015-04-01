      SUBROUTINE E01DAU(NXDATA,XDATA,X,J)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE BNSCH     BINARY SEARCH FOR INTERVAL INDEX
C     ================
C
C     CREATED 08 05 78.  UPDATED 21 06 82.  RELEASE 00/04
C
C     AUTHORS ... MAURICE G. COX, PAULINE E. M. CURTIS AND
C                 J. GEOFFREY HAYES
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX  TW11 0LW, ENGLAND
C
C     (C)  CROWN COPYRIGHT 1978-1982
C
C     **********************************************************
C
C     BNSCH.  BINARY SEARCH ROUTINE TO FIND THE INTERVAL
C     INDEX  J  SUCH THAT  XDATA(J) .LE. X .LT. XDATA(J + 1)
C
C     INPUT PARAMETERS
C        NXDATA   NUMBER OF ABSCISSA VALUES
C        XDATA    ABSCISSA VALUES
C        X        ABSCISSA VALUE
C
C     OUTPUT  PARAMETER
C        J        INTERVAL INDEX ASSOCIATED WITH  X
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           J, NXDATA
C     .. Array Arguments ..
      DOUBLE PRECISION  XDATA(NXDATA)
C     .. Local Scalars ..
      INTEGER           JHI, JLO, MIDDLE
      LOGICAL           BELOW
C     .. Executable Statements ..
C
C     INITIAL BOUNDS FOR  J
C
      JLO = 1
      JHI = NXDATA
C
C     BISECT CURRENT INTERVAL
C
   20 MIDDLE = (JLO+JHI)/2
C
C     TERMINATE SEARCH IF INTERVAL BETWEEN LOWER AND
C     UPPER BOUNDS IS SUFFICIENTLY SMALL
C
      IF (JHI-JLO.LE.1) GO TO 40
C
C     SET NEW BOUNDS
C
      BELOW = X .LT. XDATA(MIDDLE)
      IF (BELOW) JHI = MIDDLE
      IF ( .NOT. BELOW) JLO = MIDDLE
      GO TO 20
   40 J = JLO
      RETURN
C
C     END E01DAU
C
      END
