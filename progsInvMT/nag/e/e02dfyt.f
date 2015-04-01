      SUBROUTINE E02DFY(NXORDR,NXKNTS,XMIN,XMAX,XLAM,LXLAM,NYORDR,
     *                  NYKNTS,YMIN,YMAX,YLAM,LYLAM,C,IC1,IC2,LC,MX,X,
     *                  MY,Y,S,IS1,IS2,LS,YBASIS,LYBSIS,XBASIS,YDEP,
     *                  LYDEP,IYINT,IYBSI,LIYBSI)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A modified DASL routine.
C     The original routine handled interior knots only, and assumed
C     the exterior knots all to be equal to XMIN or XMAX depending
C     on the end of the knot range. The modified routine works with
C     the complete set of knots, interior plus exterior, so that
C     the user may supply exterior knots that are not equal.
C     This is not required if the spline is computed by any
C     current Nag routine, which all produce equal exterior
C     knots, but may be useful if a user has his own knot
C     set and spline. This is for compatibility with E02DBF,
C     which may be withdrawn.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE B2VR0     BIVARIATE POLYNOMIAL SPLINE
C     ================     EVALUATION ON RECTANGULAR MESH
C                          - BASE (WORKING) VERSION
C
C     CREATED 23 08 80.  UPDATED 25 12 82.  RELEASE 00/19
C
C     AUTHORS ... MAURICE G. COX AND PAULINE E. M. CURTIS.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1980-1982
C
C     **********************************************************
C
C     E02DFY.  AN ALGORITHM FOR EVALUATING A BIVARIATE
C     POLYNOMIAL SPLINE  S(X, Y)  FROM ITS B-SPLINE
C     REPRESENTATION AT ALL VERTICES OF A RECTANGULAR MESH.
C
C     THE ROUTINE IS DESIGNED TO BE FAST, BUT THIS IS AT THE
C     EXPENSE OF SOME WORKSPACE.  IF SPEED IS LESS IMPORTANT
C     THAN SPACE,  E02DFY  MAY BE CALLED SEVERAL TIMES WITH
C     ONE MESHLINE OF VALUES, FOR EXAMPLE, AT EACH CALL.
C
C     INPUT PARAMETERS
C        NXORDR   ORDER (DEGREE + 1) OF SPLINE  S  IN  X
C        NXKNTS   NUMBER OF INTERIOR AND EXTERIOR X-KNOTS
C        XMIN,
C        XMAX     LOWER AND UPPER ENDPOINTS OF  X-INTERVAL
C        XLAM     INTERIOR AND EXTERIOR X-KNOTS
C        LXLAM    DIMENSION OF  XLAM
C        NYORDR   ORDER (DEGREE + 1) OF SPLINE  S  IN  Y
C        NYKNTS   NUMBER OF INTERIOR AND EXTERIOR Y-KNOTS
C        YMIN,
C        YMAX     LOWER AND UPPER ENDPOINTS OF  Y-INTERVAL
C        YLAM     INTERIOR AND EXTERIOR Y-KNOTS
C        LYLAM    DIMENSION OF YLAM
C        C        B-SPLINE COEFFICIENTS OF  S
C        IC1,
C        IC2      INDEX INCREMENTS OF  C
C        LC       DIMENSION OF  C
C        MX       NUMBER OF  X-MESHLINES
C        X        X-MESHLINE VALUES
C        MY       NUMBER OF  Y-MESHLINES
C        Y        Y-MESHLINE VALUES
C
C     OUTPUT (AND ASSOCIATED DIMENSION) PARAMETERS
C        S        VALUES OF SPLINE
C        IS1,
C        IS2      INDEX INCREMENTS OF  S
C        LS       DIMENSION OF  S
C
C     WORKSPACE (AND ASSOCIATED DIMENSION) PARAMETERS
C        YBASIS   NON-ZERO B-SPLINE BASIS VALUES FOR ALL
C                    Y-VALUES
C        LYBSIS   DIMENSION OF  YBASIS.
C                    .GE. MY*NYORDR
C        XBASIS   NON-ZERO B-SPLINE BASIS VALUES FOR CURRENT
C                    X-VALUES
C        YDEP     Y-DEPENDENT PARTS OF BIVARIATE SPLINE
C                     S(X, Y)  NEEDED IN FORMING  S(X, Y)
C                     FOR EACH X-VALUE
C        LYDEP    DIMENSION OF  YDEP.
C                    .GE. NYKNTS + NYORDR.
C        IYINT    INTERVAL INDICES CORRESPONDING TO Y-VALUES
C        IYBSI    THE SET OF INDICES OF B-SPLINES WITH
C                    RESPECT TO  Y  WHICH CONTAIN WITHIN
C                    THEIR SUPPORT THE INTERVALS SPECIFIED
C                    BY  IYINT
C        LIYBSI   DIMENSION OF  IYBSI.
C                    .GE. NYKNTS + NYORDR.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN, YMAX, YMIN
      INTEGER           IC1, IC2, IS1, IS2, LC, LIYBSI, LS, LXLAM,
     *                  LYBSIS, LYDEP, LYLAM, MX, MY, NXKNTS, NXORDR,
     *                  NYKNTS, NYORDR
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LC), S(LS), X(MX), XBASIS(NXORDR),
     *                  XLAM(LXLAM), Y(MY), YBASIS(LYBSIS), YDEP(LYDEP),
     *                  YLAM(LYLAM)
      INTEGER           IYBSI(LIYBSI), IYINT(MY)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           IC, IND, IS, ISREF, IX, IXBAS, IY, IYBAS, J,
     *                  JSTOP, JSTRT, JXINT, JYINT, LKNOT, NINDEX
C     .. External Subroutines ..
      EXTERNAL          E01DAV, E02DFV, E02DFX
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     DETERMINE (1) THE SET OF INTERVAL INDICES
C     CORRESPONDING TO THE SPECIFIED  Y-VALUES, AND (2) THE
C     INDICES OF B-SPLINES WHICH CONTAIN WITHIN THEIR
C     SUPPORT THE INTERVALS SPECIFIED BY THE INTERVAL
C     INDICES IN (1), AND THE NUMBER OF SUCH INDICES
C
      CALL E02DFX(NYORDR,NYKNTS+NYORDR,YLAM(5),LYLAM-NYORDR,YMAX,MY,Y,
     *            IYINT,NINDEX,IYBSI)
C
C     FORM THE NON-ZERO B-SPLINE BASIS VALUES (IN  Y) FOR
C     EACH SPECIFIED  Y-VALUE
C
      LKNOT = 2*MAX(NXORDR,NYORDR)
      IYBAS = 1 - NYORDR
      DO 20 IY = 1, MY
         IYBAS = IYBAS + NYORDR
         JYINT = IYINT(IY)
C        Call to KNTAS removed from here.
         CALL E01DAV(NYORDR,YLAM(JYINT+1),LKNOT,Y(IY),YBASIS(IYBAS))
   20 CONTINUE
C
C     COUNT OVER THE SPECIFIED  X-VALUES
C
      ISREF = 1 - IS1 - IS2
      JXINT = 0
      DO 120 IX = 1, MX
         ISREF = ISREF + IS1
C
C        DETERMINE THE INTERVAL INDEX CORRESPONDING TO
C        X(IX)
C
         CALL E02DFV(NXKNTS,XLAM(5),LXLAM-NXORDR,XMAX,X(IX),JXINT)
C
C        FORM THE NON-ZERO B-SPLINE BASIS VALUES (IN  X)  AT
C        X(IX)
C
C        Call to KNTAS removed from here.
         CALL E01DAV(NXORDR,XLAM(JXINT+1),LKNOT,X(IX),XBASIS)
C
C        FORM IN  YDEP  ALL  Y-DEPENDENT PARTS OF THE
C        BIVARIATE SPLINE NEEDED IN COMPUTING THE REQUIRED
C        VALUES OF  S(X, Y)  AT  X = X(IX)
C
         DO 60 IND = 1, NINDEX
            J = IYBSI(IND)
            IC = 1 + (JXINT-1)*IC1 + (J-1)*IC2
            T = ZERO
            DO 40 IXBAS = 1, NXORDR
               IC = IC + IC1
               T = T + XBASIS(IXBAS)*C(IC)
   40       CONTINUE
            YDEP(J) = T
   60    CONTINUE
C
C        FINALLY FORM THE VALUES OF
C        S(X(IX), Y(IY)), FOR  IY = 1, 2, ... , MY
C
         IS = ISREF
         IYBAS = 0
         DO 100 IY = 1, MY
            IS = IS + IS2
            JSTRT = IYINT(IY) + 1
            JSTOP = IYINT(IY) + NYORDR
            T = ZERO
            DO 80 J = JSTRT, JSTOP
               IYBAS = IYBAS + 1
               T = T + YBASIS(IYBAS)*YDEP(J)
   80       CONTINUE
            S(IS) = T
  100    CONTINUE
  120 CONTINUE
      RETURN
C
C     END E02DFY
C
      END
