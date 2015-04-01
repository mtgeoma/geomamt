      SUBROUTINE G02GBY(N,LINK,ETA,FV,T,A,WT,NO,IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-937 (APR 1991).
C
C     CALCULATES FITTED VALUE FROM LINEAR PREDICTOR FOR
C     BINOMIAL LINKS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A
      INTEGER           IND, N, NO
      CHARACTER         LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), FV(N), T(N), WT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, BIG, C, E, P
      INTEGER           I, IFA
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF, X02AJF, X02AMF
      EXTERNAL          S15ABF, X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, SQRT
C     .. Executable Statements ..
      IND = 1
      A = X02AJF()
      IF (N.EQ.NO) THEN
         IF (LINK.EQ.'G' .OR. LINK.EQ.'g') THEN
            B = -LOG(A)
            DO 20 I = 1, N
               IF (ABS(ETA(I)).GE.B) THEN
                  IND = 0
                  RETURN
               ELSE
                  E = EXP(ETA(I))
                  FV(I) = T(I)*E/(1.0D0+E)
               END IF
   20       CONTINUE
         ELSE IF (LINK.EQ.'P' .OR. LINK.EQ.'p') THEN
            B = 1.0D0 - A
            BIG = SQRT(-2.0D0*LOG(X02AMF()))
            DO 40 I = 1, N
               IF (ETA(I).LT.-BIG .OR. ETA(I).GT.BIG) THEN
                  IND = 0
                  RETURN
               ELSE
                  IFA = 0
                  P = S15ABF(ETA(I),IFA)
                  IF (P.LT.A .OR. P.GT.B) THEN
                     IND = 0
                     RETURN
                  ELSE
                     FV(I) = T(I)*P
                  END IF
               END IF
   40       CONTINUE
         ELSE IF (LINK.EQ.'C' .OR. LINK.EQ.'c') THEN
            B = LOG(A)
            C = LOG(-B)
            DO 60 I = 1, N
               IF (ETA(I).LT.B .OR. ETA(I).GT.C) THEN
                  IND = 0
                  RETURN
               ELSE
                  FV(I) = T(I)*(1.0D0-EXP(-EXP(ETA(I))))
               END IF
   60       CONTINUE
         END IF
      ELSE
         IF (LINK.EQ.'G' .OR. LINK.EQ.'g') THEN
            B = -LOG(A)
            DO 80 I = 1, N
               IF (WT(I).NE.0.0D0 .AND. T(I).NE.0.0D0) THEN
                  IF (ABS(ETA(I)).GE.B) THEN
                     IND = 0
                     RETURN
                  ELSE
                     E = EXP(ETA(I))
                     FV(I) = T(I)*E/(1.0D0+E)
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
   80       CONTINUE
         ELSE IF (LINK.EQ.'P' .OR. LINK.EQ.'p') THEN
            B = 1.0D0 - A
            BIG = SQRT(-2.0D0*LOG(X02AMF()))
            DO 100 I = 1, N
               IF (WT(I).NE.0.0D0 .AND. T(I).NE.0.0D0) THEN
                  IF (ETA(I).LT.-BIG .OR. ETA(I).GT.BIG) THEN
                     IND = 0
                     RETURN
                  ELSE
                     IFA = 0
                     P = S15ABF(ETA(I),IFA)
                     IF (P.LT.A .OR. P.GT.B) THEN
                        IND = 0
                        RETURN
                     ELSE
                        FV(I) = T(I)*P
                     END IF
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  100       CONTINUE
         ELSE IF (LINK.EQ.'C' .OR. LINK.EQ.'c') THEN
            B = LOG(A)
            C = LOG(-B)
            DO 120 I = 1, N
               IF (WT(I).NE.0.0D0 .AND. T(I).NE.0.0D0) THEN
                  IF (ETA(I).LT.B .OR. ETA(I).GT.C) THEN
                     IND = 0
                     RETURN
                  ELSE
                     FV(I) = T(I)*(1.0D0-EXP(-EXP(ETA(I))))
                  END IF
               ELSE
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               END IF
  120       CONTINUE
         END IF
      END IF
      RETURN
      END
