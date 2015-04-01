      SUBROUTINE F02AXF(AR,IAR,AI,IAI,N,WR,VR,IVR,VI,IVI,WK1,WK2,WK3,
     *                  IFAIL)
C     MARK 3 RELEASE. NAG COPYRIGHT 1972.
C     MARK 4.5 REVISED
C     MARK 9 REVISED. IER-327 (SEP 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 14 REVISED. IER-739 (DEC 1989).
C
C     EIGENVALUES AND EIGENVECTORS OF A COMPLEX HERMITIAN MATRIX
C     1ST APRIL 1972
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AXF')
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IFAIL, IVI, IVR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), VI(IVI,N), VR(IVR,N),
     *                  WK1(N), WK2(N), WK3(N), WR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, MAX, SQ, SUM, XXXX
      INTEGER           I, ISAVE, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01BCF, F02AYF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      DO 40 I = 1, N
         IF (AI(I,I).NE.0.0D0) GO TO 140
         DO 20 J = 1, I
            VR(I,J) = AR(I,J)
            VI(I,J) = AI(I,J)
   20    CONTINUE
   40 CONTINUE
      CALL F01BCF(N,XXXX,VR,IVR,VI,IVI,WR,WK1,WK2,WK3)
      IFAIL = 1
      CALL F02AYF(N,X02AJF(),WR,WK1,VR,IVR,VI,IVI,IFAIL)
      IF (IFAIL.EQ.0) GO TO 60
      IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
C     NORMALISE
   60 DO 120 I = 1, N
         SUM = 0.0D0
         MAX = 0.0D0
         DO 80 J = 1, N
            SQ = VR(J,I)*VR(J,I) + VI(J,I)*VI(J,I)
            SUM = SUM + SQ
            IF (SQ.LE.MAX) GO TO 80
            MAX = SQ
            A = VR(J,I)
            B = VI(J,I)
   80    CONTINUE
         IF (SUM.EQ.0.0D0) GO TO 120
         SUM = 1.0D0/SQRT(SUM*MAX)
         DO 100 J = 1, N
            SQ = SUM*(VR(J,I)*A+VI(J,I)*B)
            VI(J,I) = SUM*(VI(J,I)*A-VR(J,I)*B)
            VR(J,I) = SQ
  100    CONTINUE
  120 CONTINUE
      RETURN
  140 IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
      END
