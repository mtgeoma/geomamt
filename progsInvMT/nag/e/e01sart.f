      INTEGER FUNCTION E01SAR(NVERTX,NABOR,IADJ,IEND)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: INDEX.
C     This function returns the index of NABOR in the
C     adjacency list for NVERTX.
C
C     Input Parameters - NVERTX - node whose adjacency list is
C                             to be searched.
C
C                     NABOR - node whose index is to be
C                             returned.  NABOR must be
C                             connected to NVERTX.
C
C                      IADJ - set of adjacency lists.
C
C                      IEND - pointers to the ends of
C                             adjacency lists in IADJ.
C
C     Input parameters are not altered by this function.
C
C     Output Parameter -  E01SAR - IADJ(INDEX) = NABOR.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     NB =   local copy of NABOR
C     INDX = index for IADJ
C
C     .. Scalar Arguments ..
      INTEGER                 NABOR, NVERTX
C     .. Array Arguments ..
      INTEGER                 IADJ(*), IEND(*)
C     .. Local Scalars ..
      INTEGER                 INDX, NB
C     .. Executable Statements ..
      NB = NABOR
C
C     Initialization
C
      INDX = IEND(NVERTX) + 1
   20 CONTINUE
C
C     Search the list of NVERTX neighbors for NB
C
      INDX = INDX - 1
      IF (IADJ(INDX).NE.NB) GO TO 20
C
      E01SAR = INDX
      RETURN
      END
