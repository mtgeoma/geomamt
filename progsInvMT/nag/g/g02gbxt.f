      SUBROUTINE G02GBX(N,LINK,ETA,T,DER,A,WT,NO)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C
C     CALCULATES THE DERIVATIVE OF THE LINK FUNCTION
C     FOR BINOMIAL LINKS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A
      INTEGER           N, NO
      CHARACTER         LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  DER(N), ETA(N), T(N), WT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CONST, E, Z
      INTEGER           I
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, SQRT
C     .. Executable Statements ..
      IF (N.EQ.NO) THEN
         IF (LINK.EQ.'G' .OR. LINK.EQ.'g') THEN
            DO 20 I = 1, N
               E = EXP(ETA(I))
               DER(I) = T(I)*E/((1.0D0+E)*(1.0D0+E))
   20       CONTINUE
         ELSE IF (LINK.EQ.'P' .OR. LINK.EQ.'p') THEN
            CONST = 1.0D0/SQRT(2.0D0*X01AAF(Z))
            DO 40 I = 1, N
               Z = -0.5D0*ETA(I)*ETA(I)
               DER(I) = T(I)*CONST*EXP(Z)
   40       CONTINUE
         ELSE IF (LINK.EQ.'C' .OR. LINK.EQ.'c') THEN
            DO 60 I = 1, N
               E = -EXP(ETA(I))
               DER(I) = -T(I)*E*EXP(E)
   60       CONTINUE
         END IF
      ELSE
         IF (LINK.EQ.'G' .OR. LINK.EQ.'g') THEN
            DO 80 I = 1, N
               IF (WT(I).NE.0.0D0 .AND. T(I).NE.0.0D0) THEN
                  E = EXP(ETA(I))
                  DER(I) = T(I)*E/((1.0D0+E)*(1.0D0+E))
               ELSE
                  DER(I) = 0.0D0
               END IF
   80       CONTINUE
         ELSE IF (LINK.EQ.'P' .OR. LINK.EQ.'p') THEN
            CONST = 1.0D0/SQRT(2.0D0*X01AAF(Z))
            DO 100 I = 1, N
               IF (WT(I).NE.0.0D0 .AND. T(I).NE.0.0D0) THEN
                  Z = -0.5D0*ETA(I)*ETA(I)
                  DER(I) = T(I)*CONST*EXP(Z)
               ELSE
                  DER(I) = 0.0D0
               END IF
  100       CONTINUE
         ELSE IF (LINK.EQ.'C' .OR. LINK.EQ.'c') THEN
            DO 120 I = 1, N
               IF (WT(I).NE.0.0D0 .AND. T(I).NE.0.0D0) THEN
                  E = -EXP(ETA(I))
                  DER(I) = -T(I)*E*EXP(E)
               ELSE
                  DER(I) = 0.0D0
               END IF
  120       CONTINUE
         END IF
      END IF
      RETURN
      END
