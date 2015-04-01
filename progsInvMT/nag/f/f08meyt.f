      SUBROUTINE F08MEY(F,G,H,SSMIN,SSMAX,SNR,CSR,SNL,CSL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             DLASV2(F,G,H,SSMIN,SSMAX,SNR,CSR,SNL,CSL)
C
C  Purpose
C  =======
C
C  DLASV2 computes the singular value decomposition of a 2-by-2
C  triangular matrix
C     [  F   G  ]
C     [  0   H  ].
C  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
C  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
C  right singular vectors for abs(SSMAX), giving the decomposition
C
C     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
C     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
C
C  Arguments
C  =========
C
C  F       (input) DOUBLE PRECISION
C          The (1,1) entry of the 2-by-2 matrix.
C
C  G       (input) DOUBLE PRECISION
C          The (1,2) entry of the 2-by-2 matrix.
C
C  H       (input) DOUBLE PRECISION
C          The (2,2) entry of the 2-by-2 matrix.
C
C  SSMIN   (output) DOUBLE PRECISION
C          abs(SSMIN) is the smaller singular value.
C
C  SSMAX   (output) DOUBLE PRECISION
C          abs(SSMAX) is the larger singular value.
C
C  SNL     (output) DOUBLE PRECISION
C  CSL     (output) DOUBLE PRECISION
C          The vector (CSL, SNL) is a unit left singular vector for the
C          singular value abs(SSMAX).
C
C  SNR     (output) DOUBLE PRECISION
C  CSR     (output) DOUBLE PRECISION
C          The vector (CSR, SNR) is a unit right singular vector for the
C          singular value abs(SSMAX).
C
C  Further Details
C  ===============
C
C  Any input parameter may be aliased with any output parameter.
C
C  Barring over/underflow and assuming a guard digit in subtraction, all
C  output quantities are correct to within a few units in the last
C  place (ulps).
C
C  In IEEE arithmetic, the code works correctly if one matrix entry is
C  infinite.
C
C  Overflow will not occur unless the largest singular value itself
C  overflows or is within a few ulps of overflow. (On machines with
C  partial overflow, like the Cray, overflow may occur if the largest
C  singular value is within a factor of 2 of overflow.)
C
C  Underflow is harmless if underflow is gradual. Otherwise, results
C  may correspond to a matrix modified by perturbations of size near
C  the underflow threshold.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  TWO
      PARAMETER         (TWO=2.0D0)
      DOUBLE PRECISION  FOUR
      PARAMETER         (FOUR=4.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
C     .. Local Scalars ..
      DOUBLE PRECISION  A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M,
     *                  MM, R, S, SLT, SRT, T, TEMP, TSIGN, TT
      INTEGER           PMAX
      LOGICAL           GASMAL, SWAP
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
C
      FT = F
      FA = ABS(FT)
      HT = H
      HA = ABS(H)
C
C     PMAX points to the maximum absolute entry of matrix
C       PMAX = 1 if F largest in absolute values
C       PMAX = 2 if G largest in absolute values
C       PMAX = 3 if H largest in absolute values
C
      PMAX = 1
      SWAP = (HA.GT.FA)
      IF (SWAP) THEN
         PMAX = 3
         TEMP = FT
         FT = HT
         HT = TEMP
         TEMP = FA
         FA = HA
         HA = TEMP
C
C        Now FA .ge. HA
C
      END IF
      GT = G
      GA = ABS(GT)
      IF (GA.EQ.ZERO) THEN
C
C        Diagonal matrix
C
         SSMIN = HA
         SSMAX = FA
         CLT = ONE
         CRT = ONE
         SLT = ZERO
         SRT = ZERO
      ELSE
         GASMAL = .TRUE.
         IF (GA.GT.FA) THEN
            PMAX = 2
            IF ((FA/GA).LT.X02AJF()) THEN
C
C              Case of very large GA
C
               GASMAL = .FALSE.
               SSMAX = GA
               IF (HA.GT.ONE) THEN
                  SSMIN = FA/(GA/HA)
               ELSE
                  SSMIN = (FA/GA)*HA
               END IF
               CLT = ONE
               SLT = HT/GT
               SRT = ONE
               CRT = FT/GT
            END IF
         END IF
         IF (GASMAL) THEN
C
C           Normal case
C
            D = FA - HA
            IF (D.EQ.FA) THEN
C
C              Copes with infinite F or H
C
               L = ONE
            ELSE
               L = D/FA
            END IF
C
C           Note that 0 .le. L .le. 1
C
            M = GT/FT
C
C           Note that abs(M) .le. 1/macheps
C
            T = TWO - L
C
C           Note that T .ge. 1
C
            MM = M*M
            TT = T*T
            S = SQRT(TT+MM)
C
C           Note that 1 .le. S .le. 1 + 1/macheps
C
            IF (L.EQ.ZERO) THEN
               R = ABS(M)
            ELSE
               R = SQRT(L*L+MM)
            END IF
C
C           Note that 0 .le. R .le. 1 + 1/macheps
C
            A = HALF*(S+R)
C
C           Note that 1 .le. A .le. 1 + abs(M)
C
            SSMIN = HA/A
            SSMAX = FA*A
            IF (MM.EQ.ZERO) THEN
C
C              Note that M is very tiny
C
               IF (L.EQ.ZERO) THEN
                  T = SIGN(TWO,FT)*SIGN(ONE,GT)
               ELSE
                  T = GT/SIGN(D,FT) + M/T
               END IF
            ELSE
               T = (M/(S+T)+M/(R+L))*(ONE+A)
            END IF
            L = SQRT(T*T+FOUR)
            CRT = TWO/L
            SRT = T/L
            CLT = (CRT+SRT*M)/A
            SLT = (HT/FT)*SRT/A
         END IF
      END IF
      IF (SWAP) THEN
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      END IF
C
C     Correct signs of SSMAX and SSMIN
C
      IF (PMAX.EQ.1) TSIGN = SIGN(ONE,CSR)*SIGN(ONE,CSL)*SIGN(ONE,F)
      IF (PMAX.EQ.2) TSIGN = SIGN(ONE,SNR)*SIGN(ONE,CSL)*SIGN(ONE,G)
      IF (PMAX.EQ.3) TSIGN = SIGN(ONE,SNR)*SIGN(ONE,SNL)*SIGN(ONE,H)
      SSMAX = SIGN(SSMAX,TSIGN)
      SSMIN = SIGN(SSMIN,TSIGN*SIGN(ONE,F)*SIGN(ONE,H))
      RETURN
C
C     End of F08MEY (DLASV2)
C
      END
