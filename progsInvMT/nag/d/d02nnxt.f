      SUBROUTINE D02NNX(N,ITOL,RTOL,ATOL,YCUR,EWT,IEWSET)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     OLD NAME EWSET
C
C-----------------------------------------------------------------------
C THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR EWT ACCORDING TO
C     EWT(I) = 1.0 /MAX ( RTOL(I)*ABS(YCUR(I)) , ATOL(I) ),  I = 1,...,N
C WITH THE SUBSCRIPT ON RTOL AND/OR ATOL POSSIBLY REPLACED BY 1 ABOVE,
C DEPENDING ON THE VALUE OF ITOL. THE ERROR WEIGHTS ARE FIRST FORMED
C AND THEN INVERTED . IN THE CASE WHEN A ZERO WEIGHT IS FORMED THE
C ROUTINE RETURNS WITH IEWSET = -1 , OTHERWISE IEWSET =1.
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IEWSET, ITOL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), EWT(N), RTOL(*), YCUR(N)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DRELPR, DUNFLO
      INTEGER           IOVFLO
C     .. Local Scalars ..
      DOUBLE PRECISION  ATOLI, RTOLI
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /FD02NM/DUNFLO, DRELPR, IOVFLO
C     .. Save statement ..
      SAVE              /FD02NM/
C     .. Executable Statements ..
      IEWSET = 0
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
C205
C205  THE FOLLOWING LOOP DOES NOT APPEAR TO BE VECTORISABLE
C205  SINCE WE REQUIRE TO SAVE THE LOOP VARIABLE WHEN
C205  EWT(I) .LT. DUNFLO
C205
      DO 20 I = 1, N
         IF (ITOL.GE.3) RTOLI = RTOL(I)
         IF (ITOL.EQ.2 .OR. ITOL.EQ.4) ATOLI = ATOL(I)
         EWT(I) = RTOLI*ABS(YCUR(I)) + ATOLI
         IF (EWT(I).LT.DUNFLO) IEWSET = -I
   20 CONTINUE
      IF (IEWSET.LT.0) GO TO 60
C
      IEWSET = 1
C
C   INVERT THE NON-ZERO ERROR WEIGHTS
C
      DO 40 I = 1, N
         EWT(I) = 1.0D0/EWT(I)
   40 CONTINUE
   60 CONTINUE
C
      RETURN
      END
