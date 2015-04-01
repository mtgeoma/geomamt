      SUBROUTINE E01SAU(KK,I1,I2,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: BDYADD.
C     This routine adds a boundary node to a triangulation
C     of a set of KK-1 points in the plane.  IADJ and IEND are
C     updated with the insertion of node KK.
C
C     Input Parameters -   KK - index of an exterior node to be
C                           added.  KK .ge. 4.
C
C                      I1 - first (rightmost as viewed from
C                           KK) boundary node in the mesh
C                           which is visible from KK - the
C                           line segment KK-I1 intersects
C                           no arcs.
C
C                      I2 - last (leftmost) boundary node
C                           which is visible from KK.
C
C                    IADJ - set of adjacency lists of nodes
C                           in the mesh.
C
C                    IEND - pointers to the ends of
C                           adjacency lists in IADJ for
C                           each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY and must contain
C     the vertices I1 and I2.  I1 and I2 may be determined by
C     E01SAW.
C
C     KK, I1, and I2 are not altered by this routine.
C
C     Output Parameters - IADJ,IEND - updated with the addition
C                                 of node KK as the last
C                                 entry.  Node KK will be
C                                 connected to I1, I2, and
C                                 all boundary nodes between
C                                 them.  No optimization of
C                                 the mesh is performed.
C
C     Module referenced by E01SAU - E01SAQ
C
C     ***********************************************************
C
C     Local Parameters -
C
C     K =            local copy of KK
C     KM1 =          K - 1
C     NRIGHT,NLEFT = local copies of I1, I2
C     NF,NL =        indices of IADJ bounding the portion of the
C                  array to be shifted
C     N1 =           IADJ index of the first neighbor of NLEFT
C     N2 =           IADJ index of the last neighbor of NRIGHT
C     I =            do-loop index
C     IMIN,IMAX =    bounds on do-loop index -- first and last
C                  elements of IEND to be incremented
C     KEND =         pointer to the last neighbor of K in IADJ
C     NEXT =         next boundary node to be connected to KK
C     INDX =         index for IADJ
C
C     .. Scalar Arguments ..
      INTEGER           I1, I2, KK
C     .. Array Arguments ..
      INTEGER           IADJ(*), IEND(KK)
C     .. Local Scalars ..
      INTEGER           I, IMAX, IMIN, INDX, K, KEND, KM1, N1, N2, NEXT,
     *                  NF, NL, NLEFT, NRIGHT
C     .. External Subroutines ..
      EXTERNAL          E01SAQ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
      K = KK
      KM1 = K - 1
      NRIGHT = I1
      NLEFT = I2
C
C     Initialize variables
C
      NL = IEND(KM1)
      N1 = 1
      IF (NLEFT.NE.1) N1 = IEND(NLEFT-1) + 1
      N2 = IEND(NRIGHT)
      NF = MAX(N1,N2)
C
C     Insert K as a neighbor of max(NRIGHT,NLEFT)
C
      CALL E01SAQ(NF,NL,2,IADJ)
      IADJ(NF+1) = K
      IMIN = MAX(NRIGHT,NLEFT)
      DO 20 I = IMIN, KM1
         IEND(I) = IEND(I) + 2
   20 CONTINUE
C
C     Initialize KEND and insert K as a neighbor of
C     min(NRIGHT,NLEFT)
C
      KEND = NL + 3
      NL = NF - 1
      NF = MIN(N1,N2)
      CALL E01SAQ(NF,NL,1,IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = MIN(NRIGHT,NLEFT)
      DO 40 I = IMIN, IMAX
         IEND(I) = IEND(I) + 1
   40 CONTINUE
C
C     INSERT NRIGHT AS THE FIRST NEIGHBOR OF K
C
      IADJ(KEND) = NRIGHT
C
C     Initialize INDX for loop on boundary nodes between nright
C     and NLEFT
C
      INDX = IEND(NRIGHT) - 2
   60 CONTINUE
      NEXT = IADJ(INDX)
      IF (NEXT.NE.NLEFT) THEN
C
C        Connect NEXT and K
C
         KEND = KEND + 1
         IADJ(KEND) = NEXT
         INDX = IEND(NEXT)
         IADJ(INDX) = K
         INDX = INDX - 1
         GO TO 60
      END IF
C
C     Insert NLEFT and 0 as the last neighbors of K
C
      IADJ(KEND+1) = NLEFT
      KEND = KEND + 2
      IADJ(KEND) = 0
      IEND(K) = KEND
      RETURN
      END
