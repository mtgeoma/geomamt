      SUBROUTINE F06GUF(NZ,Y,X,INDX)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PURPOSE
C
C         ZGTHR GATHERS THE SPECIFIED ELEMENTS FROM
C             A COMPLEX VECTOR  Y  IN FULL STORAGE FORM
C         INTO
C             A COMPLEX VECTOR  X  IN COMPRESSED FORM (X,INDX).
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN INDX
C         ARE REFERENCED.
C
C     ARGUMENTS
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO
C                             COMPRESSED FORM.
C         Y       COMPLEX     ARRAY, ON INPUT, WHICH CONTAINS THE
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS CORRESPONDING TO THE INDICES
C                             IN  INDX  WILL BE ACCESSED.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED
C                             FORM.
C
C     OUTPUT ...
C
C         X       COMPLEX     ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     .. Entry Points ..
      ENTRY             ZGTHR(NZ,Y,X,INDX)
C     .. Scalar Arguments ..
      INTEGER           NZ
C     .. Array Arguments ..
      COMPLEX*16        X(*), Y(*)
      INTEGER           INDX(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
C
      IF (NZ.LE.0) RETURN
C
      DO 20 I = 1, NZ
         X(I) = Y(INDX(I))
   20 CONTINUE
C
      RETURN
      END
