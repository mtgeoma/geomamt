      SUBROUTINE G05CAZ(INIT)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     called by G05CAF, G05CBF, G05CCF, G05CFZ, G05CGZ, G05DGF, G05FAF,
C     G05FBF AND G05FDF to ensure that the contents of common blocks
C     /AG05CA/, /BG05CA/, /CG05CA/ and /DG05CA/ are initialized.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     This version of G05CAZ must be used in conjunction with the
C     new auxiliary routine G05CAY which has been introduced at Mark 14.
C
C     These notes are intended to guide implementors through the text
C     changes necessary to implement the basic random number generator
C     routines G05CAY, G05CAZ, G05CBF, G05CCF, G05CFZ, G05CGZ. Please
C     follow these guidelines, and consult NAG Central Office if in any
C     doubt or difficulty. Please send a listing of your final text for
C     these routines to Central Office.
C
C     1.  Prepare code for G05CAY following guidelines supplied there.
C
C     2.  Read "DETAILS-NOTE-1" below.
C
C     3.  Activate all lines beginning CAnn, where nn is the value of
C         ILIM used in G05CAY.
C
C     ******************************************************************
C
C     ************************ DETAILS-NOTE-1 **************************
C
C     G05CAZ must be implemented consistently with G05CAY.
C
C     If G05CAY has been implemented simply by selecting suitable
C     variant code according to the value of ILIM, then a consistent
C     implementation of G05CAY may be obtained by using the variant
C     code supplied in comments beginning CAnn where the digits nn
C     are the value of ILIM.
C
C     If G05CAY has been implemented in machine code, it will still
C     be possible on many machines to implement G05CAZ in Fortran
C     and this will be satisfactory since it is not important for
C     G05CAZ to be particularly efficient. Essentially the code for
C     G05CAZ depends only on how the internal variable N is stored in
C     the array B in the common block /AG05CA/ and the code given
C     below should be applicable provided that N is stored in
C     accordance with a particular value of ILIM as defined in the
C     text of G05CAY.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER           LV
      PARAMETER         (LV=63)
      INTEGER           ILIM
      PARAMETER         (ILIM=4)
CA03  PARAMETER         (ILIM=3)
CA02  PARAMETER         (ILIM=2)
C     .. Scalar Arguments ..
      LOGICAL           INIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  GAMMA, NORMAL, VNORML
      INTEGER           DEFOPT, OPTION, POSSOP
C     .. Arrays in Common ..
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      LOGICAL           INIT2
C     .. External Subroutines ..
      EXTERNAL          G05CAY
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /BG05CA/NORMAL, GAMMA
      COMMON            /DG05CA/VNORML
C     .. Save statement ..
      SAVE              INIT2, /AG05CA/, /BG05CA/, /DG05CA/
C     .. Data statements ..
      DATA              INIT2/.TRUE./
C     .. Executable Statements ..
C
C     If INIT2 is not already .FALSE. , initialize /AG05CA/, /BG05CA/
C     and /DG05CA/ and set INIT2 to .FALSE.
C
      IF (INIT2) THEN
C
         B(0,1) =  6698
         B(0,2) =  7535
         B(0,3) = 26792
         B(0,4) = 30140
CA03     B(0,1) = 498218
CA03     B(0,2) = 172267
CA03     B(0,3) = 964506
CA02     B(0,1) = 246913578
CA02     B(0,2) = 987654312
         OPTION = 0
         DEFOPT = 0
         POSSOP = 0
C
         NORMAL = 1.0D0
         GAMMA = -1.0D0
         VNORML = 256.0D0
C
         INIT2 = .FALSE.
C
C        Initialize the buffer
C
         CALL G05CAY(.TRUE.)
      END IF
C
C     Set INIT to .FALSE. in any case
C
      INIT = .FALSE.
C
      RETURN
      END
