      SUBROUTINE G02GBU(N,LINK,Y,T,FV,ETA,WT,NO,IND)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     FINDS INTIAL VALUES OF FV AND ETA FOR BINOMIAL GLM
C
C     .. Scalar Arguments ..
      INTEGER           IND, N, NO
      CHARACTER         LINK
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), FV(N), T(N), WT(*), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  P
      INTEGER           I, IFA
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF
      EXTERNAL          G01CEF
C     .. Intrinsic Functions ..
      INTRINSIC         LOG
C     .. Executable Statements ..
      IND = 1
      IF (N.EQ.NO) THEN
         IF (LINK.EQ.'G' .OR. LINK.EQ.'g') THEN
            DO 20 I = 1, N
               FV(I) = T(I)*(Y(I)+0.5D0)/(T(I)+1.0D0)
               ETA(I) = LOG(FV(I)/(T(I)-FV(I)))
   20       CONTINUE
         ELSE IF (LINK.EQ.'P' .OR. LINK.EQ.'p') THEN
            DO 40 I = 1, N
               P = (Y(I)+0.5D0)/(T(I)+1.0D0)
               IFA = 1
               ETA(I) = G01CEF(P,IFA)
               IF (IFA.EQ.1) THEN
                  IND = 0
                  RETURN
               END IF
               FV(I) = T(I)*P
   40       CONTINUE
         ELSE IF (LINK.EQ.'C' .OR. LINK.EQ.'c') THEN
            DO 60 I = 1, N
               FV(I) = T(I)*(Y(I)+0.5D0)/(T(I)+1.0D0)
               ETA(I) = LOG(-LOG(1.0D0-FV(I)/T(I)))
   60       CONTINUE
         END IF
      ELSE
         IF (LINK.EQ.'G' .OR. LINK.EQ.'g') THEN
            DO 80 I = 1, N
               IF (WT(I).EQ.0.0D0 .OR. T(I).EQ.0.0D0) THEN
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               ELSE
                  FV(I) = T(I)*(Y(I)+0.5D0)/(T(I)+1.0D0)
                  ETA(I) = LOG(FV(I)/(T(I)-FV(I)))
               END IF
   80       CONTINUE
         ELSE IF (LINK.EQ.'P' .OR. LINK.EQ.'p') THEN
            DO 100 I = 1, N
               IF (WT(I).EQ.0.0D0 .OR. T(I).EQ.0.0D0) THEN
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               ELSE
                  P = (Y(I)+0.5D0)/(T(I)+1.0D0)
                  IFA = 1
                  ETA(I) = G01CEF(P,IFA)
                  IF (IFA.EQ.1) THEN
                     IND = 0
                     RETURN
                  END IF
                  FV(I) = T(I)*P
               END IF
  100       CONTINUE
         ELSE IF (LINK.EQ.'C' .OR. LINK.EQ.'c') THEN
            DO 120 I = 1, N
               IF (WT(I).EQ.0.0D0 .OR. T(I).EQ.0.0D0) THEN
                  FV(I) = 0.0D0
                  ETA(I) = 0.0D0
               ELSE
                  FV(I) = T(I)*(Y(I)+0.5D0)/(T(I)+1.0D0)
                  ETA(I) = LOG(-LOG(1.0D0-FV(I)/T(I)))
               END IF
  120       CONTINUE
         END IF
      END IF
      END
