      SUBROUTINE E01SAY(N,X,Y,IADJ,IEND,IER)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C     Original name: TRMESH.
C     This routine creates a Thiessen triangulation of N
C     arbitrarily spaced points in the plane referred to as
C     nodes.  The triangulation is optimal in the sense that it
C     is as nearly equiangular as possible.  E01SAY is part of
C     an interpolation package which also provides subroutines
C     to reorder the nodes, add a new node, delete an arc, plot
C     the mesh, and print the data structure.
C     Unless the nodes are already ordered in some reasonable
C     fashion, they should be sorted into ascending order for
C     increased efficiency before calling E01SAY.
C
C     Input Parameters -     N - number of nodes in the mesh.
C                            N .ge. 3.
C
C                      X,Y - N-vectors of coordinates.
C                            (X(I),Y(I)) defines node I.
C
C                     IADJ - vector of length .ge. 6*N-9.
C
C                     IEND - vector of length .ge. N.
C
C     N, X, and Y are not altered by this routine.
C
C     Output Parameters - IADJ - adjacency lists of neighbors in
C                            counterclockwise order.  The
C                            list for node I+1 follows that
C                            for node I where X and Y define
C                            the order.  The value 0 denotes
C                            the boundary (or a pseudo-node
C                            at infinity) and is always the
C                            last neighbor of a boundary
C                            node.  IADJ is unchanged if IER
C                            .ne. 0.
C
C                     IEND - pointers to the ends of
C                            adjacency lists (sets of
C                            neighbors) in IADJ.  The
C                            neighbors of node 1 begin in
C                            iadj(1).  For K .gt. 1, the
C                            neighbors of node K begin in
C                            IADJ(IEND(K-1)+1) and K has
C                            IEND(K) - IEND(K-1) neighbors
C                            including (possibly) the
C                            boundary.  IADJ(IEND(K)) .eq. 0
C                            iff node K is on the boundary.
C                            IEND is unchanged if IER = 1.
C                            If IER = 2 IEND contains the
C                            indices of a sequence of N
C                            nodes ordered from left to
C                            right where left and right are
C                            defined by assuming node 1 is
C                            to the left of node 2.
C
C                      IER - error indicator
C                            IER = 0 if no errors were
C                                    encountered.
C                            IER = 1 if N .lt. 3.
C                            IER = 2 if N .ge. 3 and all
C                                    nodes are collinear.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NN =          local copy of N
C     K =           node (index) to be inserted into IEND
C     KM1 =         K-1 - (variable) length of IEND
C     NL,NR =       IEND(1), IEND(KM1) -- leftmost and rightmost
C                 nodes in IEND as viewed from the right of
C                 1-2 when IEND contains the initial ordered
C                 set of nodal indices
C     XL,YL,XR,YR = X and Y coordinates of NL and NR
C     DXR,DYR =     XR-XL, YR-YL
C     XK,YK =       X and Y coordinates of node K
C     DXK,DYK =     XK-XL, YK-YL
C     CPROD =       vector cross product of NL-NR and NL-K --
C                 used to determine the position of node K
C                 with respect to the line defined by the
C                 nodes in IEND
C     SPROD =       scalar product used to determine the
C                 interval containing node K when K is on
C                 the line defined by the nodes in IEND
C     IND,INDX =    indices for IEND and IADJ, respectively
C     N0,ITEMP =    temporary nodes (indices)
C     IERR =        dummy parameter for call to E01SAX
C     KM1D2,KMI,I = KM1/2, K-I, do-loop index -- used in IEND
C                 reordering loop
C     KMIN =        first node index sent to E01SAX
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IER, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
      INTEGER           IADJ(6*N-9), IEND(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CPROD, DXK, DXR, DYK, DYR, SPROD, XK, XL, XR,
     *                  YK, YL, YR
      INTEGER           I, IERR, IND, INDX, ITEMP, K, KM1, KM1D2, KMI,
     *                  KMIN, N0, NL, NN, NR
C     .. External Subroutines ..
      EXTERNAL          E01SAQ, E01SAX
C     .. Executable Statements ..
      NN = N
      IER = 1
      IF (NN.GE.3) THEN
         IER = 0
C
C        Initialize IEND, NL, NR, and K
C
         IEND(1) = 1
         IEND(2) = 2
         XL = X(1)
         YL = Y(1)
         XR = X(2)
         YR = Y(2)
         K = 2
   20    CONTINUE
C
C        Begin loop on nodes 3,4,...
C
         DXR = XR - XL
         DYR = YR - YL
   40    CONTINUE
C
C        Next loop begins here if NL and NR are unchanged
C
         IF (K.EQ.NN) THEN
            GO TO 260
         ELSE
            KM1 = K
            K = KM1 + 1
            XK = X(K)
            YK = Y(K)
            DXK = XK - XL
            DYK = YK - YL
            CPROD = DXR*DYK - DXK*DYR
            IF (CPROD.GT.ZERO) THEN
               GO TO 120
            ELSE IF (CPROD.LT.ZERO) THEN
               GO TO 160
            ELSE
C
C              Node K lies on the line containing nodes 1,2,...,K-1.
C              Set SPROD to (NL-NR,NL-K).
C
               SPROD = DXR*DXK + DYR*DYK
               IF (SPROD.GT.ZERO) THEN
C
C                 Node K is to the right of NL.  Find the leftmost node
C                 N0 which lies to the right of K.
C                 Set SPROD to (N0-NL,N0-K).
C
                  DO 60 IND = 2, KM1
                     N0 = IEND(IND)
                     SPROD = (XL-X(N0))*(XK-X(N0)) + (YL-Y(N0))
     *                       *(YK-Y(N0))
                     IF (SPROD.GE.ZERO) GO TO 80
   60             CONTINUE
                  GO TO 100
C
C                 Node K lies between IEND(IND-1) and IEND(IND).
C                 Insert K in IEND.
C
   80             CALL E01SAQ(IND,KM1,1,IEND)
                  IEND(IND) = K
                  GO TO 40
               END IF
            END IF
         END IF
C
C        Node K is to the left of NL.  Insert K as the first
C        (leftmost) node in IEND and set NL to K.
C
         CALL E01SAQ(1,KM1,1,IEND)
         IEND(1) = K
         XL = XK
         YL = YK
         GO TO 20
C
C        Node K is to the right of NR.  Insert K as the last
C        (rightmost) node in IEND and set NR to K.
C
  100    IEND(K) = K
         XR = XK
         YR = YK
         GO TO 20
C
C        Node K is to the left of NL-NR.  Reorder IEND so that NL
C        is the leftmost node as viewed from K.
C
  120    KM1D2 = KM1/2
         DO 140 I = 1, KM1D2
            KMI = K - I
            ITEMP = IEND(I)
            IEND(I) = IEND(KMI)
            IEND(KMI) = ITEMP
  140    CONTINUE
C
C        Node K is to the right of NL-NR.  Create a triangulation
C        consisting of nodes 1,2,...,K.
C
  160    NL = IEND(1)
         NR = IEND(KM1)
C
C        Create the adjacency lists for the first K-1 nodes.
C        Insert neighbors in reverse order.  Each node has four
C        neighbors except NL and NR which have three.
C
         DO 180 IND = 1, KM1
            N0 = IEND(IND)
            INDX = 4*N0
            IF (N0.GE.NL) INDX = INDX - 1
            IF (N0.GE.NR) INDX = INDX - 1
            IADJ(INDX) = 0
            INDX = INDX - 1
            IF (IND.LT.KM1) IADJ(INDX) = IEND(IND+1)
            IF (IND.LT.KM1) INDX = INDX - 1
            IADJ(INDX) = K
            IF (IND.NE.1) IADJ(INDX-1) = IEND(IND-1)
  180    CONTINUE
C
C        Create the adjacency list for node K
C
         INDX = 5*KM1 - 1
         IADJ(INDX) = 0
         DO 200 IND = 1, KM1
            INDX = INDX - 1
            IADJ(INDX) = IEND(IND)
  200    CONTINUE
C
C        Replace IEND elements with pointers to IADJ
C
         INDX = 0
         DO 220 IND = 1, KM1
            INDX = INDX + 4
            IF (IND.EQ.NL .OR. IND.EQ.NR) INDX = INDX - 1
            IEND(IND) = INDX
  220    CONTINUE
         INDX = INDX + K
         IEND(K) = INDX
C
C        Add the remaining nodes to the triangulation
C
         IF (K.NE.NN) THEN
            KMIN = K + 1
            DO 240 K = KMIN, NN
               CALL E01SAX(K,X,Y,IADJ,IEND,IERR)
  240       CONTINUE
         END IF
         RETURN
C
C        All nodes are collinear
C
  260    IER = 2
      END IF
      RETURN
      END
