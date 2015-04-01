      SUBROUTINE E01DAM(NORDER,NKNOTS,XMIN,XMAX,LAMBDA,LLMBDA,JINTVL,
     *                  KNOT,LKNOT)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE KNTAS     ASSEMBLE KNOTS RELATING TO
C     ================     B-SPLINES THAT ARE NONZERO
C                          WITHIN SPECIFIED INTERVAL
C
C     CREATED 08 01 81.  UPDATED 23 06 82.  RELEASE 00/04
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1982
C
C     **********************************************************
C
C     KNTAS.  ASSEMBLES, FROM THE KNOT SET
C
C        ( L(1 - N), L(2 - N), ..., L(Q) ),
C
C     WHERE  N = NORDER,  Q = NKNOTS + NORDER,  AND
C
C              ( XMIN        (              J  .LT.  1     )
C     L(J)  =  ( LAMBDA(J)   (     1  .LE.  J  .LE.  NKNOTS)
C              ( XMAX        (NKNOTS  .LT.  J              )
C
C     THE CONTIGUOUS SUBSET
C
C        (KNOT(1), KNOT(2), ..., KNOT(2*NORDER))
C
C     WHICH RELATES TO THE B-SPLINES OF ORDER  NORDER  THAT
C     ARE NONZERO WITHIN INTERVAL NUMBER  JINTVL.
C
C     INPUT PARAMETERS
C        NORDER   ORDER (DEGREE + 1)
C        NKNOTS   NUMBER OF INTERIOR KNOTS
C        XMIN,
C        XMAX     LOWER AND UPPER ENDPOINTS OF INTERVAL
C        LAMBDA   INTERIOR KNOTS
C        LLMBDA   DIMENSION OF  LAMBDA.  .GE. MAX(NKNOTS, 1)
C        JINTVL   INTERVAL NUMBER
C
C     OUTPUT (AND ASSOCIATED DIMENSION) PARAMETERS
C        KNOT     ASSEMBLED KNOT SET
C        LKNOT    DIMENSION OF  KNOT.  .GE.  2*NORDER.
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           JINTVL, LKNOT, LLMBDA, NKNOTS, NORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  KNOT(LKNOT), LAMBDA(LLMBDA)
C     .. Local Scalars ..
      INTEGER           J, K, N2
C     .. Executable Statements ..
C
      N2 = 2*NORDER
      J = JINTVL - NORDER
      DO 20 K = 1, N2
         J = J + 1
         IF (J.LT.1) KNOT(K) = XMIN
         IF (J.GE.1 .AND. J.LE.NKNOTS) KNOT(K) = LAMBDA(J)
         IF (J.GT.NKNOTS) KNOT(K) = XMAX
   20 CONTINUE
      RETURN
C
C     END E01DAM
C
      END
