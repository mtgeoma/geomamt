      SUBROUTINE E01DAR(N,IBANDW,IXNZST,XROW,NSETS,YROW,IY1,LYROW,
     *                  LASTCL,UFCTR,LUFCTR,THETA,IT1,IT2,LTHETA,WRK,
     *                  LWRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     A DASL routine, unmodified except for name, and calls to
C     DASL routines VSCADD and VDCOPY replaced by functionally
C     equivalent NAG.
C
C     **********************************************************
C
C     D A S L  -  DATA APPROXIMATION SUBROUTINE LIBRARY
C
C     SUBROUTINE ERB       ROW ELIMINATION FOR BAND UPPER
C     ==============       TRIANGULAR MATRIX STORED BY ROWS
C
C     CREATED 08 09 81.  UPDATED 24 06 82.  RELEASE 00/10
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     (C)  CROWN COPYRIGHT 1981-1982
C
C     **********************************************************
C
C     ERB.  UPDATES, USING ELIMINATION WITHOUT PIVOTING, A
C     BAND UPPER TRIANGULAR SYSTEM  (U, THETA)  BY ONE ROW
C     (XROW, YROW)
C
C     INPUT PARAMETERS
C        N        ORDER OF BAND UPPER TRIANGULAR MATRIX  U
C        IBANDW   BANDWIDTH OF  U
C        IXNZST   COLUMN POSITION OF FIRST NONZERO IN  XROW
C        XROW     OBSERVATIONAL ROW (BASIS FUNCTION VALUES)
C                    DESTROYED ON EXIT
C        NSETS    NUMBER OF RIGHT HAND SIDES (DEPENDENT
C                    VARIABLES)
C        YROW     OBSERVATIONAL ROW (DEPENDENT VARIABLE
C                    VALUES)
C        IY1      INDEX INCREMENT OF  YROW
C        LYROW    DIMENSION OF  YROW
C
C     INPUT/OUTPUT (AND ASSOCIATED) PARAMETERS
C        LASTCL   LAST COLUMN OF  UFCTR  CONTAINING
C                    NONZERO ELEMENTS
C        UFCTR    BAND UPPER TRIANGULAR MATRIX, STORED BY
C                    ROWS
C        LUFCTR   DIMENSION OF  UFCTR.
C                   .GE. IBANDW*(2*N - IBANDW + 1)/2
C        THETA    TRANSFORMED RIGHT HAND SIDE VECTORS
C        IT1,
C        IT2      INDEX INCREMENTS OF  THETA
C        LTHETA   DIMENSION OF  THETA
C
C     WORKSPACE (AND ASSOCIATED) PARAMETERS
C        WRK      WORKSPACE
C        LWRK     DIMENSION OF  WRK.  .GE. NSETS.
C
C     ----------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           IBANDW, IT1, IT2, IXNZST, IY1, LASTCL, LTHETA,
     *                  LUFCTR, LWRK, LYROW, N, NSETS
C     .. Array Arguments ..
      DOUBLE PRECISION  THETA(LTHETA), UFCTR(LUFCTR), WRK(LWRK),
     *                  XROW(N), YROW(LYROW)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, SCALE, ZERO
      INTEGER           IT, IU, IU0, IU1, J, NMBW
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Data statements ..
C
C     KEY LOCAL VARIABLES
C        IU0,
C        IU1      CONSTANTS IN FORMULA FOR  J-TH  DIAGONAL
C                    ELEMENT OF  UFCTR  (J .GT. N - IBANDW)
C        SCALE    NEGATED MULTIPLIER IN ELIMINATION STEP
C
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
C     .. Executable Statements ..
C
      NMBW = N - IBANDW
      IU0 = -NMBW**2 - 3*N + IBANDW
      IU1 = 2*N + 3
      LASTCL = MAX(LASTCL,MIN(IXNZST+IBANDW-1,N))
C
C     COPY  YROW  TO  WRK
C
C      CALL VDCOPY ( NSETS, YROW, IY1, LYROW, ONE,
C     *              WRK, 1, LWRK )
C     Replace call to DASL routine VDCOPY above by BLAS DCOPY.
      CALL DCOPY(NSETS,YROW,IY1,WRK,1)
C
C     ANNIHILATE SUCCESSIVELY ELEMENTS  IXNZST  TO  LASTCL
C
      DO 40 J = IXNZST, LASTCL
         IF (XROW(J).EQ.ZERO) GO TO 40
         IT = 1 + (J-1)*IT1
C
C        INDEX OF DIAGONAL ELEMENT ON MAIN PART OF BAND ...
C
         IF (J.LE.NMBW) IU = IBANDW*(J-1) + 1
C
C        ... OR IN BOTTOM RIGHT TRIANGLE OF ORDER  IBANDW
C
         IF (J.GT.NMBW) IU = (IU0+(IU1-J)*J)/2
         IF (UFCTR(IU).NE.ZERO) GO TO 20
C
C           TRANSFER STEP
C
C            CALL VDCOPY ( LASTCL - J + 1, XROW(J),
C        *                    1, N - J + 1, ONE, UFCTR(IU),
C        *                    1, LUFCTR - IU + 1 )
C        Replace call to DASL routine VDCOPY above by BLAS DCOPY.
         CALL DCOPY(LASTCL-J+1,XROW(J),1,UFCTR(IU),1)
C            CALL VDCOPY ( NSETS, WRK, 1, LWRK, ONE,
C        *                    THETA(IT), IT2, LTHETA - IT + 1 )
C        Replace call to DASL routine VDCOPY above by BLAS DCOPY.
         CALL DCOPY(NSETS,WRK,1,THETA(IT),IT2)
         GO TO 60
C
C           ELIMINATION STEP
C
   20    SCALE = -XROW(J)/UFCTR(IU)
C            CALL VSCADD ( SCALE,
C        *                    MIN0 ( IBANDW, LASTCL - J + 1 ),
C        *                    UFCTR(IU), 1, LUFCTR - IU + 1,
C        *                    XROW(J), 1, N - J + 1 )
C           Replace DASL routine VSCADD call above by BLAS DAXPY.
         CALL DAXPY(MIN(IBANDW,LASTCL-J+1),SCALE,UFCTR(IU),1,XROW(J),1)
C            CALL VSCADD ( SCALE, NSETS, THETA(IT), IT2,
C        *                    LTHETA - IT2 + 1, WRK, 1, LWRK )
C        Replace DASL routine VSCADD call above by BLAS DAXPY.
         CALL DAXPY(NSETS,SCALE,THETA(IT),IT2,WRK,1)
   40 CONTINUE
   60 RETURN
C
C     END E01DAR
C
      END
