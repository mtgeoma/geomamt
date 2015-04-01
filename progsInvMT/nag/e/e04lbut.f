      SUBROUTINE E04LBU(N,NUMNEG,ISTATE)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBU (RESET) RESTORES TO THEIR NORMAL VALUES THE ISTATE
C     ENTRIES OF - 4 AND - 5 CORRESPONDING TO FIXED VARIABLES WITH
C     NEGATIVE MULTIPLIERS.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           N, NUMNEG
C     .. Array Arguments ..
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      INTEGER           I, ISI
C     .. Executable Statements ..
      NUMNEG = 0
      DO 20 I = 1, N
         ISI = ISTATE(I)
         IF (ISI.GT.-4) GO TO 20
         ISTATE(I) = ISI + 3
   20 CONTINUE
      RETURN
C
C     END OF E04LBU (RESET)
C
      END
