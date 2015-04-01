      SUBROUTINE E01SAS(NIN1,NIN2,NOUT1,NOUT2,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: SWAP.
C     This subroutine swaps the diagonals in a convex quadri-
C     lateral.
C
C     Input Parameters -  NIN1,NIN2,NOUT1,NOUT2 - nodal indices
C                            of a pair of adjacent triangles
C                            which form a convex quadrilat-
C                            eral.  NOUT1 and NOUT2 are con-
C                            nected by an arc which is to be
C                            replaced by the arc NIN1-NIN2.
C                            (NIN1,NOUT1,NOUT2) must be tri-
C                            angle vertices in counterclock-
C                            wise order.
C
C     The above parameters are not altered by this routine.
C
C                iadj,iend - triangulation data structure
C                            (see subroutine E01SAY).
C
C     output parameters - IADJ,IEND - updated with the arc
C                                 replacement.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     IN =        NIN1 and NIN2 ordered by increasing magnitude
C               (the neighbors of IN(1) precede those of
C               IN(2) in IADJ)
C     IO =        NOUT1 and NOUT2 in increasing order
C     IP1,IP2 =   permutation of (1,2) such that IO(IP1)
C               precedes IO(IP2) as a neighbor of IN(1)
C     J,K =       permutation of (1,2) used as indices of in
C               and IO
C     NF,NL =     IADJ indices boundary a portion of the array
C               to be shifted
C     I =         IEND index
C     IMIN,IMAX = bounds on the portion of IEND to be incre-
C               mented or decremented
C
C     .. Scalar Arguments ..
      INTEGER           NIN1, NIN2, NOUT1, NOUT2
C     .. Array Arguments ..
      INTEGER           IADJ(*), IEND(*)
C     .. Local Scalars ..
      INTEGER           I, IMAX, IMIN, IP1, IP2, J, K, NF, NL
C     .. Local Arrays ..
      INTEGER           IN(2), IO(2)
C     .. External Functions ..
      INTEGER           E01SAR
      EXTERNAL          E01SAR
C     .. External Subroutines ..
      EXTERNAL          E01SAQ
C     .. Executable Statements ..
      IN(1) = NIN1
      IN(2) = NIN2
      IO(1) = NOUT1
      IO(2) = NOUT2
      IP1 = 1
C
C     Order the indices so that IN(1) .lt. IN(2) and IO(1) .lt.
C     IO(2), and choose IP1 and IP2 such that (IN(1),IO(IP1),
C     IO(IP2)) forms a triangle.
C
      IF (IN(1).GE.IN(2)) THEN
         IN(1) = IN(2)
         IN(2) = NIN1
         IP1 = 2
      END IF
      IF (IO(1).GE.IO(2)) THEN
         IO(1) = IO(2)
         IO(2) = NOUT1
         IP1 = 3 - IP1
      END IF
      IP2 = 3 - IP1
      IF (IO(2).LT.IN(1)) THEN
C
C        The vertices are ordered (IO(1),IO(2),IN(1),IN(2)).
C        Delete IO(2) by shifting up by 1
C
         NF = 1 + E01SAR(IO(1),IO(2),IADJ,IEND)
         NL = -1 + E01SAR(IO(2),IO(1),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,-1,IADJ)
         IMIN = IO(1)
         IMAX = IO(2) - 1
         DO 20 I = IMIN, IMAX
            IEND(I) = IEND(I) - 1
   20    CONTINUE
C
C        Delete IO(1) by shifting up by 2 and insert IN(2)
C
         NF = NL + 2
         NL = -1 + E01SAR(IN(1),IO(IP2),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,-2,IADJ)
         IADJ(NL-1) = IN(2)
         IMIN = IO(2)
         IMAX = IN(1) - 1
         DO 40 I = IMIN, IMAX
            IEND(I) = IEND(I) - 2
   40    CONTINUE
C
C        Shift up by 1 and insert IN(1)
C
         NF = NL + 1
         NL = -1 + E01SAR(IN(2),IO(IP1),IADJ,IEND)
         CALL E01SAQ(NF,NL,-1,IADJ)
         IADJ(NL) = IN(1)
         IMIN = IN(1)
         IMAX = IN(2) - 1
         DO 60 I = IMIN, IMAX
            IEND(I) = IEND(I) - 1
   60    CONTINUE
      ELSE IF (IN(2).LT.IO(1)) THEN
C
C        The vertices are ordered (IN(1),IN(2),IO(1),IO(2)).
C        Delete IO(1) by shifting down by 1
C
         NF = 1 + E01SAR(IO(1),IO(2),IADJ,IEND)
         NL = -1 + E01SAR(IO(2),IO(1),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,1,IADJ)
         IMIN = IO(1)
         IMAX = IO(2) - 1
         DO 80 I = IMIN, IMAX
            IEND(I) = IEND(I) + 1
   80    CONTINUE
C
C        Delete IO(2) by shifting down by 2 and insert IN(1)
C
         NL = NF - 2
         NF = 1 + E01SAR(IN(2),IO(IP2),IADJ,IEND)
         IF (NF.LE.NL) CALL E01SAQ(NF,NL,2,IADJ)
         IADJ(NF+1) = IN(1)
         IMIN = IN(2)
         IMAX = IO(1) - 1
         DO 100 I = IMIN, IMAX
            IEND(I) = IEND(I) + 2
  100    CONTINUE
C
C        Shift down by 1 and insert IN(2)
C
         NL = NF - 1
         NF = 1 + E01SAR(IN(1),IO(IP1),IADJ,IEND)
         CALL E01SAQ(NF,NL,1,IADJ)
         IADJ(NF) = IN(2)
         IMIN = IN(1)
         IMAX = IN(2) - 1
         DO 120 I = IMIN, IMAX
            IEND(I) = IEND(I) + 1
  120    CONTINUE
      ELSE
C
C        IN(1) and IO(1) precede IN(2) and IO(2).  For (J,K) =
C        (1,2) and (2,1), delete IO(K) as a neighbor of IO(J)
C        by shifting a portion of IADJ either up or down and
C        and insert IN(K) as a neighbor of IN(J).
C
         DO 180 J = 1, 2
            K = 3 - J
            IF (IN(J).GT.IO(J)) THEN
C
C              The neighbors of IO(J) precede those of IN(J) -- shift
C              up by 1
C
               NF = 1 + E01SAR(IO(J),IO(K),IADJ,IEND)
               NL = -1 + E01SAR(IN(J),IO(IP2),IADJ,IEND)
               IF (NF.LE.NL) CALL E01SAQ(NF,NL,-1,IADJ)
               IADJ(NL) = IN(K)
               IMIN = IO(J)
               IMAX = IN(J) - 1
               DO 140 I = IMIN, IMAX
                  IEND(I) = IEND(I) - 1
  140          CONTINUE
            ELSE
C
C              The neighbors of IN(J) precede those of IO(J) -- shift
C              down by 1
C
               NF = 1 + E01SAR(IN(J),IO(IP1),IADJ,IEND)
               NL = -1 + E01SAR(IO(J),IO(K),IADJ,IEND)
               IF (NF.LE.NL) CALL E01SAQ(NF,NL,1,IADJ)
               IADJ(NF) = IN(K)
               IMIN = IN(J)
               IMAX = IO(J) - 1
               DO 160 I = IMIN, IMAX
                  IEND(I) = IEND(I) + 1
  160          CONTINUE
            END IF
C
C           Reverse (IP1,IP2) for (J,K) = (2,1)
C
            IP1 = IP2
            IP2 = 3 - IP1
  180    CONTINUE
      END IF
      RETURN
      END
