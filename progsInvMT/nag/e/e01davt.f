      SUBROUTINE E01DAV(NORDER,KNOT,LKNOT,X,BASIS)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE BBV       EVALUATION OF B-SPLINE BASIS
C     ==============       FROM BASIC 3-TERM RECURRENCE
C                          RELATION
C
C     CREATED 09 01 81.  UPDATED 21 06 82.  RELEASE 00/07
C
C     AUTHORS ... MAURICE G. COX AND PAULINE E. M. CURTIS.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1981-1982
C
C     **********************************************************
C
C     BBV.  DETERMINES THE  NORDER  VALUES AT  X  OF THE
C     NORMALIZED B-SPLINES OF ORDER  NORDER  WHICH ARE
C     NONZERO IN THE KNOT INTERVAL CONTAINING  X.
C
C     INPUT PARAMETERS
C        NORDER   ORDER (DEGREE + 1)
C        KNOT     KNOTS RELATING TO B-SPLINES OF ORDER
C                    NORDER  THAT ARE NONZERO IN THE
C                    INTERVAL CONTAINING  X,  I.E. THE
C                    NORDER  KNOTS TO THE IMMEDIATE LEFT AND
C                    THE  NORDER  KNOTS TO THE IMMEDIATE
C                    RIGHT OF  X
C        LKNOT    DIMENSION OF  KNOT.  .GE.  2*NORDER.
C        X        ABSCISSA VALUE
C
C     OUTPUT PARAMETERS
C        BASIS    NONZERO B-SPLINE BASIS VALUES
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           LKNOT, NORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  BASIS(NORDER), KNOT(LKNOT)
C     .. Local Scalars ..
      INTEGER           IKNOT, JORDER
C     .. External Subroutines ..
      EXTERNAL          E01DAT
C     .. Executable Statements ..
C
      IKNOT = NORDER + 1
      DO 20 JORDER = 1, NORDER
         IKNOT = IKNOT - 1
         CALL E01DAT(JORDER,KNOT(IKNOT),2*JORDER,X,BASIS)
   20 CONTINUE
      RETURN
C
C     END E01DAV
C
      END
