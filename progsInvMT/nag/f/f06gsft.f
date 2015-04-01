      COMPLEX*16  FUNCTION F06GSF(NZ,X,INDX,Y)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PURPOSE
C
C         ZDOTCI COMPUTES THE CONJUGATED VECTOR INNER PRODUCT OF
C             A COMPLEX SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX)
C         WITH
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX
C         ARE REFERENCED.
C
C     ARGUMENTS
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       COMPLEX     ARRAY CONTAINING THE VALUES OF THE
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         ZDOTCI   COMPLEX    COMPLEX FUNCTION VALUE EQUAL TO THE
C                             CONJUGATED VECTOR INNER PRODUCT.
C                             IF  NZ .LE. 0  ZDOTCI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     .. Entry Points ..
      COMPLEX*16                  ZDOTCI
      ENTRY                       ZDOTCI(NZ,X,INDX,Y)
C     .. Scalar Arguments ..
      INTEGER                     NZ
C     .. Array Arguments ..
      COMPLEX*16                  X(*), Y(*)
      INTEGER                     INDX(*)
C     .. Local Scalars ..
      INTEGER                     I
C     .. Intrinsic Functions ..
      INTRINSIC                   DCONJG
C     .. Executable Statements ..
C
      F06GSF = (0.0D0,0.0D0)
      IF (NZ.LE.0) RETURN
C
      DO 20 I = 1, NZ
         F06GSF = F06GSF + DCONJG(X(I))*Y(INDX(I))
   20 CONTINUE
C
      RETURN
      END
