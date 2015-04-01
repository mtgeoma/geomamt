      SUBROUTINE F04MCV(N1,N2,U,IU1,IU2,LU,SCALE,V,IV1,IV2,LV)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ******************************************************
C
C     NPL DATA FITTING LIBRARY ROUTINE RSCOPY
C
C     CREATED 20 09 78.  UPDATED 06 09 79.  RELEASE 00/09
C
C     AUTHOR ... MAURICE G. COX.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     ******************************************************
C
C     SCALE AND COPY  N1  BY  N2  RECTANGULAR ARRAY  U  TO  V.
C     U  AND  V  ARE DESCRIBED BY STANDARD RECTANGULAR DATA
C     STRUCTURES  (U, IU1, IU2, LU)  AND
C     (V, IV1, IV2, LV).
C
C     INPUT PARAMETERS
C        N1       NUMBER OF ELEMENTS IN EACH SLICE OF  U
C        N1       NUMBER OF SLICES OF  U
C        U        ARRAY TO BE SCALED AND COPIED
C        IU1      FIRST INDEX INCREMENT FOR  U
C        IU2      SECOND INDEX INCREMENT FOR  U
C        LU       DIMENSION OF  U
C        SCALE    SCALING FACTOR
C
C     OUTPUT (AND ASSOCIATED) PARAMETERS
C        V        SCALED AND COPIED ARRAY
C        IV1      FIRST INDEX INCREMENT FOR  V
C        IV2      SECOND INDEX INCREMENT FOR  V
C        LV       DIMENSION OF  V
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SCALE
      INTEGER           IU1, IU2, IV1, IV2, LU, LV, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  U(LU), V(LV)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE
      INTEGER           I, IU, IUREF, IV, IVREF, J
C     .. Data statements ..
      DATA              ONE/1.0D+0/
C     .. Executable Statements ..
      IUREF = 1 - IU1 - IU2
      IVREF = 1 - IV1 - IV2
      IF (SCALE.NE.ONE) GO TO 60
      DO 40 J = 1, N2
         IUREF = IUREF + IU2
         IU = IUREF
         IVREF = IVREF + IV2
         IV = IVREF
         DO 20 I = 1, N1
            IU = IU + IU1
            IV = IV + IV1
            V(IV) = U(IU)
   20    CONTINUE
   40 CONTINUE
      GO TO 120
   60 DO 100 J = 1, N2
         IUREF = IUREF + IU2
         IU = IUREF
         IVREF = IVREF + IV2
         IV = IVREF
         DO 80 I = 1, N1
            IU = IU + IU1
            IV = IV + IV1
            V(IV) = SCALE*U(IU)
   80    CONTINUE
  100 CONTINUE
  120 RETURN
      END
