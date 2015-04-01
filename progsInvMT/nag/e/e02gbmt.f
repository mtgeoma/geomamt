      SUBROUTINE E02GBM(MAKE,IR,INDX,ARAY,MPL1)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     AN ADAPTATION OF D. E. KNUTH,S HEAPING
C     ROUTINES (SEE VOLUME 3 OF
C     THE ART OF COMPUTER PROGRAMMING  ).
C     IF  MAKE  IS  .TRUE.,  THE FULL HEAP BUILDING
C     PROCESS IS CARRIED OUT ON
C     ARAY(1),...,ARAY(IR) ,
C     AND THE VALUE OF  IR  IS UNCHANGED.
C     IF  MAKE  IS  .FALSE.,  ONE STEP OF THE SORTING
C     PROCESS IS CARRIED OUT TO PROVIDE THE NEXT
C     ELEMENT OF  ARAY  IN ORDER,  AND THE VARIABLE
C     IR  IS DECREASED BY ONE.  THE INTERRUPTION OF THE
C     SORTING PHASE IS BUILT IN VIA THE FLAG  ONCE.
C     INDX  IS AN INDEX VECTOR ASSOCIATED WITH
C     ARAY  WHICH MUST BE REARRANGED IN PARRALEL
C     WITH IT.
C     ***************
C
C     .. Scalar Arguments ..
      INTEGER           IR, MPL1
      LOGICAL           MAKE
C     .. Array Arguments ..
      DOUBLE PRECISION  ARAY(MPL1)
      INTEGER           INDX(MPL1)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           I, IL, IT, J
      LOGICAL           ONCE
C     .. Executable Statements ..
      IF (IR.GT.1) GO TO 20
      IF ( .NOT. MAKE) IR = 0
      RETURN
   20 CONTINUE
C
C     ***************
C     TEST WHETHER OR NOT THE INITIAL
C     HEAP IS TO BE BUILT
C     ***************
C
      IL = 1
      IF (MAKE) IL = (IR/2) + 1
      ONCE = .FALSE.
C
C     ***************
C     THE LOOP BEGINS HERE
C     ***************
C
   40 CONTINUE
      IF (IL.GT.1) GO TO 60
C
C     ***************
C     THE SORTING PHASE USES THIS BRANCH
C     ***************
C
      IF (MAKE .OR. ONCE) RETURN
      ONCE = .TRUE.
      IT = INDX(IR)
      T = ARAY(IR)
      INDX(IR) = INDX(1)
      ARAY(IR) = ARAY(1)
      IR = IR - 1
      IF (IR.GT.1) GO TO 80
      INDX(1) = IT
      ARAY(1) = T
      RETURN
   60 CONTINUE
C
C     ***************
C     THE HEAP-BUILDING PHASE USES THIS BRANCH
C     ***************
C
      IL = IL - 1
      IT = INDX(IL)
      T = ARAY(IL)
   80 CONTINUE
C
C     ***************
C     THE REMAINING STATEMENTS ARE COMMON
C     TO BOTH PHASES AND EMBODY THE
C     HEAP-RECTIFYING (SIFTING) SECTION
C     ***************
C
      J = IL
  100 CONTINUE
      I = J
      J = 2*J
      IF (J-IR) 120, 140, 160
  120 CONTINUE
      IF (ARAY(J).LE.ARAY(J+1)) GO TO 140
      J = J + 1
  140 CONTINUE
      IF (T.LE.ARAY(J)) GO TO 160
      INDX(I) = INDX(J)
      ARAY(I) = ARAY(J)
      GO TO 100
  160 CONTINUE
      INDX(I) = IT
      ARAY(I) = T
      GO TO 40
      END
