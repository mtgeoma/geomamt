      SUBROUTINE D02PVX(OUTCH,MCHEPS,DWARF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ SUBROUTINE ENVIRN $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  The RK suite requires some environmental parameters that are provided
C  by this subroutine. The values provided with the distribution codes
C  are those appropriate to the IEEE standard. They must be altered, if
C  necessary, to those appropriate to the computing system you are using
C  before calling the codes of the suite.
C
C        ===============================================================
C        ===============================================================
C        TO MAKE SURE THAT THESE MACHINE AND INSTALLATION DEPENDENT
C        QUANTITIES ARE SPECIFIED PROPERLY, THE DISTRIBUTION VERSION
C        WRITES A MESSAGE ABOUT THE MATTER TO THE STANDARD OUTPUT
C        CHANNEL AND TERMINATES THE RUN. THE VALUES PROVIDED IN THE
C        DISTRIBUTION VERSION SHOULD BE ALTERED, IF NECESSARY, AND THE
C        "WRITE" AND "STOP" STATEMENTS COMMENTED OUT.
C        ===============================================================
C        ===============================================================
C
C  OUTPUT VARIABLES
C
C     OUTCH     - INTEGER
C                 Standard output channel
C     MCHEPS    - DOUBLE PRECISION
C                 MCHEPS is the largest positive number such that
C                 1.0D0 + MCHEPS = 1.0D0.
C     DWARF     - DOUBLE PRECISION
C                 DWARF is the smallest positive number.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
C  The following six statements are to be Commented out after
C  verification that the machine and installation dependent quantities
C  are specified correctly. If you pass copies of RKSUITE on to others,
C  please give them the whole distribution version of RKSUITE, and in
C  particular, give them a version of D02PVX that does not have the
C  following six statements Commented out.
C      WRITE(*,*) ' Before using RKSUITE, you must verify that the  '
C      WRITE(*,*) ' machine- and installation-dependent quantities  '
C      WRITE(*,*) ' specified in the subroutine D02PVX are correct, '
C      WRITE(*,*) ' and then Comment these WRITE statements and the '
C      WRITE(*,*) ' STOP statement out of D02PVX.                   '
C      STOP
C
C  The following values are appropriate to IEEE arithmetic with the
C  typical standard output channel.
C
C      OUTCH = 6
C      MCHEPS = 1.11D-16
C      DWARF = 2.23D-308
C
C-----------------------------------------------------------------------
C  If you have the routines D1MACH and I1MACH on your system, you could
C  replace the preceding statements by the following ones to obtain the
C  appropriate machine dependent numbers. The routines D1MACH and I1MACH
C  are public domain software. They are available from NETLIB.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DWARF, MCHEPS
      INTEGER           OUTCH
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          X04AAF
C     .. Executable Statements ..
C
      CALL X04AAF(0,OUTCH)
      MCHEPS = X02AJF()
      DWARF = X02AMF()
C
      RETURN
      END
