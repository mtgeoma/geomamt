      SUBROUTINE F06EVF(NZ,Y,X,INDX)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     PURPOSE
C
C         DGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM
C             A REAL VECTOR  Y  IN FULL STORAGE FORM
C         INTO
C             A REAL VECTOR  X  IN COMPRESSED FORM  (X,INDX).
C         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX
C         ARE REFERENCED OR MODIFIED.
C
C     ARGUMENTS
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED
C                             FORM.
C
C     UPDATED ...
C
C         Y       REAL        ARRAY, ON INPUT, WHICH CONTAINS THE
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE
C                             GATHERED COMPONENTS IN  Y  ARE SET TO
C                             ZERO. ONLY THE ELEMENTS CORRESPONDING TO
C                             THE INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       REAL        ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     .. Entry Points ..
      ENTRY             DGTHRZ(NZ,Y,X,INDX)
C     .. Scalar Arguments ..
      INTEGER           NZ
C     .. Array Arguments ..
      DOUBLE PRECISION  X(*), Y(*)
      INTEGER           INDX(*)
C     .. Local Scalars ..
      INTEGER           I
C     .. Executable Statements ..
C
      IF (NZ.LE.0) RETURN
C
      DO 20 I = 1, NZ
         X(I) = Y(INDX(I))
         Y(INDX(I)) = 0.0D0
   20 CONTINUE
C
      RETURN
      END
