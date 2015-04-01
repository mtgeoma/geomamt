      SUBROUTINE E01DAT(NORDER,KNOT,LKNOT,X,BASIS)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE BRV       ONE MAJOR STEP OF BASIC 3-TERM
C     ==============       B-SPLINE RECURRENCE RELATION
C
C     CREATED 08 01 81.  UPDATED 21 06 82.  RELEASE 00/08
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1981-1982
C
C     **********************************************************
C
C     BRV.  GIVEN THE  NORDER - 1  VALUES AT  X  OF THE
C     NORMALIZED B-SPLINES OF ORDER  NORDER - 1  WHICH ARE
C     NONZERO IN THE KNOT INTERVAL CONTAINING  X,  BRV
C     DETERMINES THE  NORDER  VALUES AT  X  OF THE
C     NORMALIZED B-SPLINES OF ORDER  NORDER
C     WHICH ARE NONZERO IN THIS INTERVAL.
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
C     INPUT/OUTPUT PARAMETERS
C        BASIS    NONZERO B-SPLINE BASIS VALUES.
C                    NORDER - 1  VALUES ON ENTRY,
C                    NORDER      VALUES ON EXIT.
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           LKNOT, NORDER
C     .. Array Arguments ..
      DOUBLE PRECISION  BASIS(NORDER), KNOT(LKNOT)
C     .. Local Scalars ..
      DOUBLE PRECISION  BNOW, BPREV, ONE, ZERO
      INTEGER           I, KL, KR
C     .. Data statements ..
C
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
C     .. Executable Statements ..
C
C     BASIS OF ORDER UNITY
C
      IF (NORDER.EQ.1) BASIS(1) = ONE
      IF (NORDER.EQ.1) GO TO 60
C
C        BASIS OF ORDER GREATER THAN ONE
C
   20 KL = NORDER + 1
      KR = 2*NORDER
      BNOW = ZERO
      DO 40 I = 2, NORDER
         KL = KL - 1
         KR = KR - 1
         BPREV = BNOW
         BNOW = BASIS(KL-1)/(KNOT(KR)-KNOT(KL))
         BASIS(KL) = (X-KNOT(KL))*BNOW + (KNOT(KR+1)-X)*BPREV
   40 CONTINUE
      BASIS(1) = (KNOT(NORDER+1)-X)*BNOW
   60 RETURN
C
C     END E01DAT
C
      END
