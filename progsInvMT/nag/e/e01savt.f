      SUBROUTINE E01SAV(KK,I1,I2,I3,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: INTADD.
C     This routine adds an interior node to a triangulation
C     of a set of KK-1 points in the plane.  IADJ and IEND are
C     updated with the insertion of node KK in the triangle
C     whose vertices are I1, I2, and I3.
C
C     Input Parameters -        KK - index of node to be
C                                inserted.  KK .ge. 4.
C
C                     I1,I2,I3 - indices of the vertices of
C                                a triangle containing node
C                                KK -- in counterclockwise
C                                order.
C
C                         IADJ - set of adjacency lists
C                                of nodes in the mesh.
C
C                         IEND - pointers to the ends of
C                                adjacency lists in IADJ for
C                                each node in the mesh.
C
C     IADJ and IEND may be created by E01SAY and must contain
C     the vertices I1, I2, and I3.  I1,I2,I3 may be determined
C     by E01SAW.
C
C     KK, I1, I2, and I3 are not altered by this routine.
C
C     Output Parameters - IADJ,IEND - updated with the addition
C                                 of node KK as the last
C                                 entry.  Node KK will be
C                                 connected to nodes I1, I2,
C                                 and I3.  No optimization
C                                 of the mesh is performed.
C
C     Module referenced by E01SAV - E01SAQ
C
C     ***********************************************************
C
C     Local Parameters -
C
C     K =           local copy of KK
C     KM1 =         K - 1
C     N =           vector containing I1, I2, I3
C     NFT =         pointers to the tops of the 3 sets of IADJ
C                 elements to be shifted downward
C     IP1,IP2,IP3 = permutation indices for N and NFT
C     INDX =        index for IADJ and N
C     NF,NL =       indices of first and last entries in IADJ
C                 to be shifted down
C     N1,N2 =       first 2 vertices of a new triangle --
C                 (N1,N2,KK)
C     IMIN,IMAX =   bounds on do-loop index -- first and last
C                 elements of IEND to be incremented
C     I =           do-loop index
C     ITEMP =       temporary storage location
C
C     .. Scalar Arguments ..
      INTEGER           I1, I2, I3, KK
C     .. Array Arguments ..
      INTEGER           IADJ(*), IEND(KK)
C     .. Local Scalars ..
      INTEGER           I, IMAX, IMIN, INDX, IP1, IP2, IP3, ITEMP, K,
     *                  KM1, N1, N2, NF, NL
C     .. Local Arrays ..
      INTEGER           N(3), NFT(3)
C     .. External Subroutines ..
      EXTERNAL          E01SAQ
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      K = KK
C
C     Initialization
C
      N(1) = I1
      N(2) = I2
      N(3) = I3
C
C     set up NFT
C
      DO 40 I = 1, 3
         N1 = N(I)
         INDX = MOD(I,3) + 1
         N2 = N(INDX)
         INDX = IEND(N1) + 1
   20    CONTINUE
C
C        Find the index of N2 as a neighbor of N1
C
         INDX = INDX - 1
         IF (IADJ(INDX).NE.N2) GO TO 20
         NFT(I) = INDX + 1
   40 CONTINUE
C
C     Order the vertices by decreasing magnitude.
C     N(IP(I+1)) precedes N(IP(I)) in IEND for
C     I = 1,2.
C
      IP1 = 1
      IP2 = 2
      IP3 = 3
      IF (N(2).GT.N(1)) THEN
         IP1 = 2
         IP2 = 1
      END IF
      IF (N(3).GT.N(IP1)) THEN
         IP3 = IP1
         IP1 = 3
      END IF
      IF (N(IP3).GT.N(IP2)) THEN
         ITEMP = IP2
         IP2 = IP3
         IP3 = ITEMP
      END IF
C
C     Add node K to the adjacency lists of each vertex and
C     update IEND.  For each vertex, a set of IADJ elements
C     is shifted downward and K is inserted.  Shifting starts
C     at the end of the array.
C
      KM1 = K - 1
      NL = IEND(KM1)
      NF = NFT(IP1)
      IF (NF.LE.NL) CALL E01SAQ(NF,NL,3,IADJ)
      IADJ(NF+2) = K
      IMIN = N(IP1)
      IMAX = KM1
      DO 60 I = IMIN, IMAX
         IEND(I) = IEND(I) + 3
   60 CONTINUE
C
      NL = NF - 1
      NF = NFT(IP2)
      CALL E01SAQ(NF,NL,2,IADJ)
      IADJ(NF+1) = K
      IMAX = IMIN - 1
      IMIN = N(IP2)
      DO 80 I = IMIN, IMAX
         IEND(I) = IEND(I) + 2
   80 CONTINUE
C
      NL = NF - 1
      NF = NFT(IP3)
      CALL E01SAQ(NF,NL,1,IADJ)
      IADJ(NF) = K
      IMAX = IMIN - 1
      IMIN = N(IP3)
      DO 100 I = IMIN, IMAX
         IEND(I) = IEND(I) + 1
  100 CONTINUE
C
C     Add node K to IEND and its neighbors to IADJ
C
      INDX = IEND(KM1)
      IEND(K) = INDX + 3
      DO 120 I = 1, 3
         INDX = INDX + 1
         IADJ(INDX) = N(I)
  120 CONTINUE
      RETURN
      END
