      SUBROUTINE D03PRQ(A,B,N,NOLD,NV)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C    Auxiliary routine to do put the contents of vector B into vector A.
C     The last NV components of A are passed form the last NV
C     components of the new vector A , after allowing for N > NOLD .
C     This allows remesh to be used with mixed ODE/PDE problems.
C
C     SPSWOP routine from SPRINT.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Scalar Arguments ..
      INTEGER           N, NOLD, NV
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*)
C     .. Scalars in Common ..
      INTEGER           IDEV, ITRACE
C     .. Local Scalars ..
      INTEGER           I, NOPNV, NPNV
      CHARACTER*100     REC
C     .. External Subroutines ..
      EXTERNAL          X04BAF
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
C     .. Save statement ..
      SAVE              /AD02NM/
C     .. Executable Statements ..
      NPNV = N + NV
      NOPNV = NOLD + NV
      IF (ITRACE.GE.3) THEN
         WRITE (REC,FMT=99999) (A(I),I=1,NOPNV)
         CALL X04BAF(IDEV,REC)
         WRITE (REC,FMT=99998) (B(I),I=1,NPNV)
         CALL X04BAF(IDEV,REC)
      END IF
      IF (NV.GT.0) THEN
         DO 20 I = 1, NV
            B(N+I) = A(NOLD+I)
   20    CONTINUE
      END IF
      DO 40 I = 1, NPNV
         A(I) = B(I)
   40 CONTINUE
      RETURN
C
99999 FORMAT (' OLD IS',1X,6D10.3)
99998 FORMAT (' NEW IS',1X,6D10.3)
      END
