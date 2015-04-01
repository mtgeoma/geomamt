      SUBROUTINE D03PRM(X,NPTS,XFIX,NXFIX,IXFIX,NF)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C  ---------------------------------------------------------------------
C     Routine to determine which of the mesh points are fixed.
C     FIXSET routine from SPRINT.
C
C     Parameter list
C     **************
C
C     NPTS    ; The number of mesh points
C     X       ; Array of length NPTS containing the mesh points
C     XFIX    ; Array of length NFIX containing the mesh points
C             which must remain fixed.
C     IXFIX   ; The array that is set by this routine as follows
C             IXFIX(I) = K when X(K) = XFIX(I)
C     NF      ; = 1 on normal exit
C           ; = -1 on exit if NXFIX > 0 and no K can be found such that
C               X(K) = XFIX(I)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      INTEGER           NF, NPTS, NXFIX
C     .. Array Arguments ..
      DOUBLE PRECISION  X(NPTS), XFIX(*)
      INTEGER           IXFIX(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           IOVFLO
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP, URTEMP
      INTEGER           I, J, K
      CHARACTER*200     ERRMSG
C     .. External Subroutines ..
      EXTERNAL          D02NNQ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
C     .. Save statement ..
      SAVE              /CD03PC/
C     .. Executable Statements ..
      URTEMP = UROUND*1000.0D0
      NF = 1
      IF (NXFIX.EQ.0) RETURN
      K = 1
      DO 60 I = 1, NXFIX
         IXFIX(I) = 0.0D0
         DO 20 J = 2, NPTS - 1
            TEMP = ABS(XFIX(I)-X(J))
            IF (TEMP.LE.URTEMP) THEN
               IXFIX(I) = J
               K = J + 1
               GO TO 40
            END IF
   20    CONTINUE
   40    CONTINUE
         IF (IXFIX(I).EQ.0.0D0) THEN
            NF = -1
            ERRMSG =
     *' Routine cannot match fixed mesh point
     *  (=I1), with value (=R1), with any of the supplied
     *  spatial mesh points. '
            CALL D02NNQ(ERRMSG,1,1,I,0,1,XFIX(I),0.0D0)
         END IF
   60 CONTINUE
      RETURN
      END
