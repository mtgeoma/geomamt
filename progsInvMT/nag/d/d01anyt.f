      DOUBLE PRECISION FUNCTION D01ANY(X,OMEGA,P2,P3,P4,INTEGR)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QWGTO.
C
C        THIS FUNCTION SUBPROGRAM IS USED IN CONJUNCTION
C        WITH THE ROUTINE  D01ANF,  AND DEFINES THE WEIGHT
C        FUNCTIONS.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 OMEGA, P2, P3, P4, X
      INTEGER                          INTEGR
C     .. Local Scalars ..
      DOUBLE PRECISION                 OMX
C     .. Intrinsic Functions ..
      INTRINSIC                        COS, SIN
C     .. Executable Statements ..
      OMX = OMEGA*X
      GO TO (20,40) INTEGR
   20 D01ANY = COS(OMX)
      GO TO 60
   40 D01ANY = SIN(OMX)
   60 RETURN
      END
