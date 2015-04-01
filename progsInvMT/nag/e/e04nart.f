      SUBROUTINE E04NAR(ORTHOG,N,X,LENX,INCX,Y,LENY,INCY,CS,SN)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C *********************************************************************
C     IF  ORTHOG  IS TRUE, E04NAR  APPLIES A PLANE ROTATION.  OTHERWISE,
C     E04NAR COMPUTES THE TRANSFORMATION (X Y)*E  AND RETURNS THE RESULT
C     IN  (X Y),  WHERE THE 2 BY 2 MATRIX  E  IS DEFINED BY  CS  AND  SN
C     AS FOLLOWS...
C
C     E  =  ( 1  SN )  IF  CS .GT. ZERO,    E  =  (     1 )  OTHERWISE.
C           (     1 )                             ( 1  SN )
C
C     VERSION 1, APRIL 5 1983.
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C
C     THIS VERSION 17-JANUARY-1984.  SVEN.
C *********************************************************************
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CS, SN
      INTEGER           INCX, INCY, LENX, LENY, N
      LOGICAL           ORTHOG
C     .. Array Arguments ..
      DOUBLE PRECISION  X(LENX), Y(LENY)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
C     .. External Subroutines ..
      EXTERNAL          DSWAP, F06FPF, DAXPY
C     .. Executable Statements ..
      IF (ORTHOG) GO TO 20
      ZERO = 0.0D0
      IF (CS.LE.ZERO) CALL DSWAP(N,X,INCX,Y,INCY)
      IF (SN.NE.ZERO) CALL DAXPY(N,SN,X,INCX,Y,INCY)
      RETURN
   20 CONTINUE
      CALL F06FPF(N,X,INCX,Y,INCY,CS,SN)
      RETURN
C
C     END OF E04NAR  ( ELM )
      END
