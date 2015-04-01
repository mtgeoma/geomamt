      SUBROUTINE G01ARQ(Y,N,MED,HL,HH,ADJL,ADJH,IADJL,IADJH,STEP)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Get general information about Y().  Useful for plot scaling.
C     Sorts Y() and returns it sorted.  Also returns :-
C       HL  =  low hinge         MED =  median      HH = hi hinge
C      ADJL =  low adjacent value                 ADJH = hi adj value
C     IADJL =  its index                         IADJH = its index
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ADJH, ADJL, HH, HL, MED, STEP
      INTEGER           IADJH, IADJL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  HFENCE, LFENCE
      INTEGER           IFAIL, J, K, TEMP1, TEMP2
C     .. External Subroutines ..
      EXTERNAL          M01CAF
C     .. Executable Statements ..
C
      IFAIL = 0
      CALL M01CAF(Y,1,N,'A',IFAIL)
C
C     M01CAF sorts the data into ascending order.  It cannot fail.
C
      K = N
      J = (K/2) + 1
C
      TEMP1 = N + 1 - J
      MED = (Y(J)+Y(TEMP1))/2.0D0
C
      K = (K+1)/2
      J = (K/2) + 1
      TEMP1 = K + 1 - J
      HL = (Y(J)+Y(TEMP1))/2.0D0
      TEMP1 = N - K + J
      TEMP2 = N + 1 - J
      HH = (Y(TEMP1)+Y(TEMP2))/2.0D0
C
      STEP = (HH-HL)*1.5D0
      HFENCE = HH + STEP
      LFENCE = HL - STEP
C
C        Find adjacent values
C
      IADJL = 0
   20 CONTINUE
      IADJL = IADJL + 1
      IF (IADJL.LT.1) THEN
         IADJL = 1
      ELSE IF (IADJL.GT.N) THEN
         IADJL = N
      ELSE IF (Y(IADJL).LT.LFENCE) THEN
         GO TO 20
      END IF
      ADJL = Y(IADJL)
C
      IADJH = N + 1
   40 CONTINUE
      IADJH = IADJH - 1
      IF (IADJH.LT.1) THEN
         IADJH = 1
      ELSE IF (IADJH.GT.N) THEN
         IADJH = N
      ELSE IF (Y(IADJH).GT.HFENCE) THEN
         GO TO 40
      END IF
      ADJH = Y(IADJH)
      RETURN
      END
