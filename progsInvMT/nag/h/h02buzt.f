      SUBROUTINE H02BUZ(MAXM,A,M,N,IOBJ,BL,BU,CRNAME)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Scalar Arguments ..
      INTEGER           IOBJ, M, MAXM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(MAXM,N), BL(N+M), BU(N+M)
      CHARACTER*8       CRNAME(N+M)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      IF (IOBJ.NE.M) THEN
         DO 40 I = IOBJ, M - 1
            BL(N+I) = BL(N+I+1)
            BU(N+I) = BU(N+I+1)
            CRNAME(N+I) = CRNAME(N+I+1)
            DO 20 J = 1, N
               A(I,J) = A(I+1,J)
   20       CONTINUE
   40    CONTINUE
      END IF
      RETURN
      END
