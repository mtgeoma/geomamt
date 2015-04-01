      DOUBLE PRECISION FUNCTION F06ERF(NZ,X,INDX,Y)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PURPOSE
C
C         DDOTI COMPUTES THE VECTOR INNER PRODUCT OF
C             A REAL SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX)
C         WITH
C             A REAL VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX
C         ARE REFERENCED.
C
C     ARGUMENTS
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       REAL        ARRAY CONTAINING THE VALUES OF THE
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         DDOTI   REAL        REAL FUNCTION VALUE EQUAL TO THE
C                             VECTOR INNER PRODUCT.
C                             IF  NZ .LE. 0  DDOTI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     .. Entry Points ..
      DOUBLE PRECISION                 DDOTI
      ENTRY                            DDOTI(NZ,X,INDX,Y)
C     .. Scalar Arguments ..
      INTEGER                          NZ
C     .. Array Arguments ..
      DOUBLE PRECISION                 X(*), Y(*)
      INTEGER                          INDX(*)
C     .. Local Scalars ..
      INTEGER                          I
C     .. Executable Statements ..
C
      F06ERF = 0.0D0
      IF (NZ.LE.0) RETURN
C
      DO 20 I = 1, NZ
         F06ERF = F06ERF + X(I)*Y(INDX(I))
   20 CONTINUE
C
      RETURN
      END
