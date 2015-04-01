      SUBROUTINE D03PHG(NRMESH,DXMESH,TRMESH,IPMINF,XRATIO,CONST,NPTS,
     *                  IND,ISET)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C     ******* LOTS OF CHANGES BY VP *******
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine sets up COMMON blocks and checks input parameters.
C     Called by D03PPF or D03PRF.
C     All parameters except ISET must be supplied to D03PPF or D03PRF by
C     the user.
C     N.B. The user is allowed to change the following parameters
C        between calls to D03PPF/PRF (i.e. when restarting the
C        integration in time):
C        NRMESH, DXMESH, TRMESH, XRATIO, CONST.
C
C
C     NRMESH  < 0  means that the new mesh will be determined solely by
C                using the parameter DXMESH. As testing for a new mesh
C                after every timestep can be expensive the mesh is
C                tested every [NRMESH] steps.
C
C           = 0  remeshing will take place just once at the end of
C                the first time step when t > TRMESH.
C
C           > 0  remeshing will take place every NRMESH timesteps
C                (with NO testing using DXMESH).
C
C     DXMESH  : parameter used to determine when remeshing should take
C             place when NRMESH is set < 0
C             a possible new mesh is calculated at the end of every
C             [NRMESH] time step but is only adopted if
C
C               (new)    (old)            (old)  (old)
C              X      > X     + DXMESH ( X    - X     )
C               i        i                i+1    i
C
C          OR
C
C               (new)    (old)             (old)   (old)
C              X      < X      - DXMESH ( X    -  X      )
C               i        i                 i       i-1
C
C     DXMESH thus controls how much the mesh can drift from one
C     remesh to the next.     (SEE NOTE BELOW)
C
C     N.B.   DXMESH imposes a LOWER bound on the difference between one
C         mesh and the next, i.e. new mesh is only adopted if the
C         change is large enough to satisfy above criteria.
C         Note that DXMESH imposes no upper bound. Large changes should
C         be avoided by regular updating.
C         Note also that DXMESH = 0 (i.e. no lower bound) is valid,
C         but it would be more efficient to choose NRMESH positive.
C                                                  --- VP ---
C
C    TRMESH ; parameter used to determine when remeshing will take place
C           when NRMESH = 0. Remeshing will occur when at the first
C           time level reached when  t , the integration time , is
C           greater than TRMESH  (with no remeshing thereafter).
C           (Useful if want to repeat for another time interval with
C            an improved mesh, or, as in Sprint manual, if know speed
C            of wave front... ?).
C
C    IPMINF ; controls print out of mesh information. IPMINF=0 no print,
C         = 1 summary print of mesh characterisitics e.g. mesh ratios &
C           max. integral of monitor per interval (approx.=input const)
C         =2 print of monitors, old & new meshes and mesh sizes
C
C     VP (SEPT 1991) NEXT PARAMETER NO LONGER PRESENT...
C   NPUSER ; this allows the number of mesh points to be varied when the
C          remesh routines are called. If it is set to zero then the
C          remesh routines will pick the number of points required,
C          subject to the MAXNPT constraint in the setup routine D03PHR.
C          This parameter is set non-zero if NRMESH < 0.
C
C  XRATIO ; input bound on adjacent mesh ratio ( > 1 typically 1.5 TO 3)
C          the remsh routines will try and ensure that
C
C             X    - X    <  XRATIO (X - X   )
C              i+1    i               i   i-1
C     and
C             X    - X    <  XRATIO (X   - X   )
C              i      i-1             i+1   i
C
C
C     CONST ;  remeshing parameter described in user manual part 2
C          if in doubt set const to   2/(NPTS -1 ) . This routine
C          will check to ensure that
C
C            0.1/ (NPTS-1)  <   CONST < 10 / (NPTS-1)
C
C    ISET   -  set to -1 if there is an error in one of input parameters
C           set to  0 otherwise
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONST, DXMESH, TRMESH, XRATIO
      INTEGER           IND, IPMINF, ISET, NPTS, NRMESH
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, DYMESH, KONST, STPRAT, TOUT, UROUND
      INTEGER           ICOUNT, IOVFLO, JPMINF, MRMESH, NOPUSR
      LOGICAL           REMESH
      CHARACTER*6       PDCODE, RMTYPE
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
C     .. External Subroutines ..
      EXTERNAL          D02NNQ
C     .. Common blocks ..
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /GD03PC/KONST, STPRAT, JPMINF, NOPUSR
      COMMON            /HD03PC/TOUT, REMESH, MRMESH, ICOUNT
      COMMON            /MD03PC/DYMESH
C     .. Save statement ..
      SAVE              /CD03PC/, /GD03PC/, /HD03PC/, /MD03PC/, /FD03PC/
C     .. Executable Statements ..
      ISET = 0
C
      IF (NRMESH.LT.0) THEN
         RMTYPE = 'REMSET'
      ELSE
         RMTYPE = 'SETOFF'
      END IF
C
      REMESH = .TRUE.
C
C     first call ..
      IF (IND.EQ.0) THEN
         ICOUNT = 0
         GO TO 20
      END IF
C
C     later calls ..
C
      IF (NRMESH.EQ.0) THEN
         IF (TRMESH.GT.TOUT) THEN
            ICOUNT = 0
         END IF
      END IF
C
      IF (NRMESH.NE.MRMESH) THEN
         ICOUNT = 0
      END IF
C
   20 MRMESH = NRMESH
      JPMINF = IPMINF
C
      IF (IPMINF.LT.0 .OR. IPMINF.GT.2) THEN
         CALL D02NNQ(
     *' Routine entered with illegal value for IPMINF (=I1). IPMINF
     *  should be in the range 0 to 2.',1,1,IPMINF,0,0,0.0D0,0.0D0)
         ISET = -1
         JPMINF = 0
      END IF
C
      IF (XRATIO.LE.1.0D0) THEN
         CALL D02NNQ(
     *' Routine entered with an illegal value for XRATIO (=R1). XRATIO
     *  should be greater than 1.0D0 ',1,0,0,0,1,XRATIO,0.0D0)
         ISET = -1
C         XRATIO = 1.5D0
      END IF
C
      STPRAT = XRATIO
      NOPUSR = NPTS
C
      TEMP = CONST*(NPTS-1)
      IF (TEMP.LT.0.1D0 .OR. TEMP.GT.10.D0) THEN
         CALL D02NNQ(
     *' Routine entered with an illegal value of CONST, which
     *  should lie between 0.1/(NPTS-1) and 10/(NPTS-1)',1,0,0,0,0,
     *               0.0D0,0.0D0)
         ISET = -1
CRWB next line commented out to preserve input only nature of CONST
CRWB 22/4/93
C         CONST = 2.D0/(NPTS-1)
      END IF
      KONST = CONST
C
      IF (MRMESH.LT.0.0D0) THEN
C
C        (VP)   CONST = 2.D0/(NPTS-1)      (User can choose this)
C        NOPUSR = 1                 (no longer used)
         DYMESH = DXMESH
         IF (DYMESH.LT.0.0D0) THEN
            CALL D02NNQ(
     *' Routine entered with an illegal value for DXMESH (=R1). DXMESH
     *  should be greater than 0.0D0 ',1,0,0,0,1,DXMESH,0.0D0)
            ISET = -1
            DYMESH = 1.5D0
         END IF
      END IF
C
      IF (MRMESH.EQ.0) THEN
         TOUT = TRMESH
      END IF
C
      RETURN
      END
