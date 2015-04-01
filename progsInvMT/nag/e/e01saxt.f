      SUBROUTINE E01SAX(KK,X,Y,IADJ,IEND,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: ADNODE.
C     This routine adds node KK to a triangulation of a set
C     of points in the plane producing a new triangulation.  A
C     sequence of edge swaps is then applied to the mesh,
C     resulting in an optimal triangulation.  E01SAX is part
C     of an interpolation package which also provides routines
C     to initialize the data structure, plot the mesh, and
C     delete arcs.
C
C     Input Parameters -   KK - index of the node to be added
C                           to the mesh.  KK .ge. 4.
C
C                     X,Y - vectors of coordinates of the
C                           nodes in the mesh.  (X(I),Y(I))
C                           defines node I for I = 1,..,KK.
C
C                    IADJ - set of adjacency lists of nodes
C                           1,..,KK-1.
C
C                    IEND - pointers to the ends of
C                           adjacency lists in IADJ for
C                           each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY.
C
C     KK, X, and Y are not altered by this routine.
C
C     Output Parameters - IADJ,IEND - updated with the addition
C                                 of node KK as the last
C                                 entry.
C
C                           IER - error indicator
C                                 IER = 0 if no errors
C                                         were encountered.
C                                 IER = 1 if all nodes
C                                         (including KK) are
C                                         collinear.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     K =        local copy of KK
C     KM1 =      K - 1
C     I1,I2,I3 = vertices of a triangle containing K
C     INDKF =    IADJ index of the first neighbor of K
C     INDKL =    IADJ index of the last neighbor of K
C     NABOR1 =   first neighbor of K before any swaps occur
C     IO1,IO2 =  adjacent neighbors of K defining an arc to
C              be tested for a swap
C     IN1 =      vertex opposite K -- first neighbor of IO2
C              which precedes IO1.  IN1,IO1,IO2 are in
C              counterclockwise order.
C     INDK1 =    index of IO1 in the adjacency list for K
C     IND2F =    index of the first neighbor of IO2
C     IND21 =    index of IO1 in the adjacency list for IO2
C     XK,YK =    X(K), Y(K)
C
C     .. Scalar Arguments ..
      INTEGER           IER, KK
C     .. Array Arguments ..
      DOUBLE PRECISION  X(KK), Y(KK)
      INTEGER           IADJ(*), IEND(KK)
C     .. Local Scalars ..
      DOUBLE PRECISION  XK, YK
      INTEGER           I1, I2, I3, IN1, IND21, IND2F, INDK1, INDKF,
     *                  INDKL, IO1, IO2, K, KM1, NABOR1
C     .. External Functions ..
      INTEGER           E01SAR
      LOGICAL           E01SAT
      EXTERNAL          E01SAR, E01SAT
C     .. External Subroutines ..
      EXTERNAL          E01SAS, E01SAU, E01SAV, E01SAW
C     .. Executable Statements ..
      IER = 0
      K = KK
C
C     Initialization
C
      KM1 = K - 1
      XK = X(K)
      YK = Y(K)
C
C     Add node K to the mesh
C
      CALL E01SAW(KM1,XK,YK,X,Y,IADJ,IEND,I1,I2,I3)
      IF (I1.EQ.0) THEN
C
C        All nodes are collinear
C
         IER = 1
      ELSE
         IF (I3.LE.0) THEN
            CALL E01SAU(K,I1,I2,IADJ,IEND)
         ELSE
            CALL E01SAV(K,I1,I2,I3,IADJ,IEND)
         END IF
C
C        Initialize variables for optimization of the mesh
C
         INDKF = IEND(KM1) + 1
         INDKL = IEND(K)
         NABOR1 = IADJ(INDKF)
         IO2 = NABOR1
         INDK1 = INDKF + 1
         IO1 = IADJ(INDK1)
   20    CONTINUE
C
C        Begin loop -- find the vertex opposite K
C
         IND2F = 1
         IF (IO2.NE.1) IND2F = IEND(IO2-1) + 1
         IND21 = E01SAR(IO2,IO1,IADJ,IEND)
         IF (IND2F.EQ.IND21) THEN
C
C           IN1 is the last neighbor of IO2
C
            IND21 = IEND(IO2)
            IN1 = IADJ(IND21)
            IF (IN1.EQ.0) GO TO 40
         ELSE
            IN1 = IADJ(IND21-1)
         END IF
C
C        Swap test -- if a swap occurs, two new arcs are opposite K
C              and must be tested.  INDK1 and INDKF must be
C              decremented.
C
         IF (E01SAT(IN1,K,IO1,IO2,X,Y)) THEN
            CALL E01SAS(IN1,K,IO1,IO2,IADJ,IEND)
            IO1 = IN1
            INDK1 = INDK1 - 1
            INDKF = INDKF - 1
            GO TO 20
         END IF
C
C        No swap occurred.  Reset IO2 and IO1, and test for
C        termination.
C
   40    IF (IO1.NE.NABOR1) THEN
            IO2 = IO1
            INDK1 = INDK1 + 1
            IF (INDK1.GT.INDKL) INDK1 = INDKF
            IO1 = IADJ(INDK1)
            IF (IO1.NE.0) GO TO 20
         END IF
      END IF
      RETURN
      END
