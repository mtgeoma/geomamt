      SUBROUTINE D03FAZ(LBDCND,L,C1,MBDCND,M,C2,NBDCND,N,A,B,C,LDIMF,
     *                  MDIMF,F,IERROR,W,JWORK1,JWORK2,JTRIGL,JTRIGM)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C1, C2
      INTEGER           IERROR, JTRIGL, JTRIGM, JWORK1, JWORK2, L,
     *                  LBDCND, LDIMF, M, MBDCND, MDIMF, N, NBDCND
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), B(N), C(N), F(LDIMF,MDIMF,*), W(*)
C     .. Local Scalars ..
      INTEGER           I, ITRIGL, ITRIGM, IWORK1, IWORK2, IWYRT, J, K,
     *                  NH, NODD
C     .. Local Arrays ..
      DOUBLE PRECISION  SAVE(6)
C     .. External Subroutines ..
      EXTERNAL          D03FAY
C     .. Executable Statements ..
C
C     Check for invalid input.
C
C     IERROR = 0
C     IF (NBDCND.EQ.0) THEN
C        DO 20 K = 1, N
C           IF (A(K).NE.C(1)) IERROR = 9
C           IF (C(K).NE.C(1)) IERROR = 9
C           IF (B(K).NE.B(1)) IERROR = 9
C     20    CONTINUE
C     END IF
C     IF (NBDCND.EQ.1 .AND. (A(1).NE.0.0D0 .OR. C(N).NE.0.0D0))
C     *    IERROR = 10
C
C     Partition workspace
C
      IWYRT = 1 + L
      IWORK1 = IWYRT + M
      IWORK2 = IWORK1 + JWORK1
      ITRIGL = IWORK2 + JWORK2
      ITRIGM = ITRIGL + JTRIGL
      IF (NBDCND.EQ.0) THEN
C
C        Reorder unknowns when NBDCND = 0
C
         NH = (N+1)/2
         NODD = 1
         IF (2*NH.EQ.N) NODD = 2
         DO 80 I = 1, L
            DO 60 J = 1, M
               DO 20 K = 1, NH - 1
                  W(K) = F(I,J,NH-K) - F(I,J,NH+K)
                  W(NH+K) = F(I,J,NH-K) + F(I,J,NH+K)
   20          CONTINUE
               W(NH) = 2.0D0*F(I,J,NH)
               IF (NODD.EQ.2) W(N) = 2.0D0*F(I,J,N)
               DO 40 K = 1, N
                  F(I,J,K) = W(K)
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
         SAVE(1) = C(NH-1)
         SAVE(2) = A(NH)
         SAVE(3) = C(NH)
         SAVE(4) = B(NH-1)
         SAVE(5) = B(N)
         SAVE(6) = A(N)
         C(NH-1) = 0.0D0
         A(NH) = 0.0D0
         C(NH) = 2.0D0*C(NH)
         IF (NODD.EQ.1) THEN
            B(NH-1) = B(NH-1) - A(NH-1)
            B(N) = B(N) + A(N)
         ELSE
            A(N) = C(NH)
         END IF
      END IF
      CALL D03FAY(LBDCND,L,MBDCND,M,N,A,B,C,LDIMF,MDIMF,F,W,W(IWYRT),C1,
     *            C2,W(IWORK1),JWORK1,W(IWORK2),JWORK2,W(ITRIGL),JTRIGL,
     *            W(ITRIGM),JTRIGM)
      IF (NBDCND.EQ.0) THEN
         DO 160 I = 1, L
            DO 140 J = 1, M
               DO 100 K = 1, NH - 1
                  W(NH-K) = 0.5D0*(F(I,J,NH+K)+F(I,J,K))
                  W(NH+K) = 0.5D0*(F(I,J,NH+K)-F(I,J,K))
  100          CONTINUE
               W(NH) = 0.5D0*F(I,J,NH)
               IF (NODD.EQ.2) W(N) = 0.5D0*F(I,J,N)
               DO 120 K = 1, N
                  F(I,J,K) = W(K)
  120          CONTINUE
  140       CONTINUE
  160    CONTINUE
         C(NH-1) = SAVE(1)
         A(NH) = SAVE(2)
         C(NH) = SAVE(3)
         B(NH-1) = SAVE(4)
         B(N) = SAVE(5)
         A(N) = SAVE(6)
      END IF
      RETURN
      END
