      SUBROUTINE E01SAQ(NFRST,NLAST,KK,IARR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     ***********************************************************
C
C                                               Robert Renka
C                                       Oak Ridge Natl. Lab.
C
C     Original name: SHIFTD.
C     This routine shifts a set of contiguous elements of an
C     integer array KK positions downward (upward if KK .lt. 0).
C     The loops are unrolled in order to increase efficiency.
C
C     Input parameters - NFRST,NLAST - bounds on the portion of
C                                  IARR to be shifted.  All
C                                  elements between and
C                                  including the bounds are
C                                  shifted unless NFRST .gt.
C                                  NLAST, in which case no
C                                  shift occurs.
C
C                             KK - number of positions each
C                                  element is to be shifted.
C                                  if KK .lt. 0 shift up.
C                                  if KK .gt. 0 shift down.
C
C                           IARR - integer array of length
C                                  .ge. NLAST + max(KK,0).
C
C     NFRST, NLAST, and KK are not altered by this routine.
C
C     Output Parameter -        IARR - shifted array.
C
C     ***********************************************************
C
C     Local Parameters -
C
C     INC =  do-loop increment (unrolling factor) -- if INC is
C          changed, statements must be added to or deleted
C          from the do-loops
C     K =    local copy of KK
C     NF =   local copy of NFRST
C     NL =   local copy of NLAST
C     NLP1 = NL + 1
C     NS =   number of shifts
C     NSL =  number of shifts done in unrolled do-loop (multiple
C          of INC)
C     I =    do-loop index and index for IARR
C     IBAK = index for downward shift of IARR
C     INDX = index for IARR
C     IMAX = bound on do-loop index
C
C     .. Scalar Arguments ..
      INTEGER           KK, NFRST, NLAST
C     .. Array Arguments ..
      INTEGER           IARR(*)
C     .. Local Scalars ..
      INTEGER           I, IBAK, IMAX, INC, INDX, K, NF, NL, NLP1, NS,
     *                  NSL
C     .. Data statements ..
      DATA              INC/5/
C     .. Executable Statements ..
      K = KK
      NF = NFRST
      NL = NLAST
      IF (NF.LE.NL .AND. K.NE.0) THEN
         NLP1 = NL + 1
         NS = NLP1 - NF
         NSL = INC*(NS/INC)
         IF (K.LT.0) THEN
C
C           Shift upward starting from the top
C
            IF (NSL.GT.0) THEN
               IMAX = NLP1 - INC
               DO 20 I = NF, IMAX, INC
                  INDX = I + K
                  IARR(INDX) = IARR(I)
                  IARR(INDX+1) = IARR(I+1)
                  IARR(INDX+2) = IARR(I+2)
                  IARR(INDX+3) = IARR(I+3)
                  IARR(INDX+4) = IARR(I+4)
   20          CONTINUE
            END IF
C
C           Perform the remaining NS-NSL shifts one at a time
C
            I = NSL + NF
   40       CONTINUE
            IF (I.LE.NL) THEN
               INDX = I + K
               IARR(INDX) = IARR(I)
               I = I + 1
               GO TO 40
            END IF
         ELSE
C
C           Shift downward starting from the bottom
C
            IF (NSL.GT.0) THEN
               DO 60 I = 1, NSL, INC
                  IBAK = NLP1 - I
                  INDX = IBAK + K
                  IARR(INDX) = IARR(IBAK)
                  IARR(INDX-1) = IARR(IBAK-1)
                  IARR(INDX-2) = IARR(IBAK-2)
                  IARR(INDX-3) = IARR(IBAK-3)
                  IARR(INDX-4) = IARR(IBAK-4)
   60          CONTINUE
            END IF
C
C           Perform the remaining NS-NSL shifts one at a time
C
            IBAK = NLP1 - NSL
   80       CONTINUE
            IF (IBAK.GT.NF) THEN
               IBAK = IBAK - 1
               INDX = IBAK + K
               IARR(INDX) = IARR(IBAK)
               GO TO 80
            END IF
         END IF
      END IF
      RETURN
      END
