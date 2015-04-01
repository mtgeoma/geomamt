      SUBROUTINE Y90RUW(N,P,NNZ,A,IROW,ICOL,IDIMA,ISTR,ISTC,AP)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
*-----------------------------------------------------------------------
*
*         ====================================
*         *  Y90RUW :  Auxiliary for Y90RUF  *
*         ====================================
*
*     Purpose
*     =======
*     Maps row P of a real symmetric matrix stored in dense format in AP
*     into linked-list format.
*
*-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           IDIMA, N, NNZ, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IDIMA), AP(N)
      INTEGER           ICOL(IDIMA), IROW(IDIMA), ISTC(N), ISTR(N)
C     .. Local Scalars ..
      INTEGER           COL, I, J, ROW
C     .. Executable Statements ..
*
*     Copy non-fill values back into linked list form of row P.
*     Set values in dense format to zero to indicate that they
*     have been copied.
*
      I = ISTR(P)
   20 CONTINUE
      J = I
   40 CONTINUE
      J = ICOL(J)
      IF (J.GT.0) GO TO 40
      J = -J
      A(I) = AP(J)
      IF (J.NE.P) AP(J) = ZERO
      I = IROW(I)
      IF (I.GT.0) GO TO 20
*
      I = ISTC(P)
   60 CONTINUE
      J = I
   80 CONTINUE
      J = IROW(J)
      IF (J.GT.0) GO TO 80
      J = -J
      A(I) = AP(J)
      AP(J) = ZERO
      I = ICOL(I)
      IF (I.GT.0) GO TO 60
*
*     Add fill elements to linked-list form of row P.
*
      DO 100 I = 1, N
         IF (AP(I).NE.ZERO) THEN
            NNZ = NNZ + 1
            A(NNZ) = AP(I)
            ROW = P
            COL = I
            IF (I.GT.P) THEN
               ROW = I
               COL = P
            END IF
            IROW(NNZ) = ISTR(ROW)
            ISTR(ROW) = NNZ
            ICOL(NNZ) = ISTC(COL)
            ISTC(COL) = NNZ
         END IF
  100 CONTINUE
*
*     End of subroutine Y90RUW
*
      RETURN
      END
