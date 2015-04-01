      SUBROUTINE G05CBF(I)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     initializes the internal variables used by the generator
C     routines to produce a repeatable sequence
C
C     G05CBF initializes the notional internal variable N in G05CAY
C     to 2*ABS(I)+1 and then calls G05CAY to fill the buffer with the
C     next LV pseudo-random numbers.
C
C     G05CBF also re-initializes the varibles NORMAL, GAMMA and VNORML
C     which are used by G05DDF, G05DGF and G05FDF respectively.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     This version of G05CBF must be used in conjunction with the
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
C     G05CBF must be implemented consistently with G05CAY.
C
C     If G05CAY has been implemented simply by selecting suitable
C     variant code according to the value of ILIM, then a consistent
C     implementation of G05CBF may be obtained by using the variant
C     code supplied in comments beginning CAnn where the digits nn
C     are the value of ILIM.
C
C     If G05CAY has been implemented in machine code, it will still
C     be possible on many machines to implement G05CBF in Fortran
C     and this will be satisfactory since it is not important for
C     G05CBF to be particularly efficient. Essentially the code for
C     G05CBF depends only on how the internal variable N is stored in
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
      INTEGER           IPW2, IPW2A
      PARAMETER         (IPW2=2**(60/ILIM),IPW2A=IPW2/4)
C     .. Scalar Arguments ..
      INTEGER           I
C     .. Scalars in Common ..
      DOUBLE PRECISION  GAMMA, NORMAL, VNORML
      INTEGER           DEFOPT, OPTION, POSSOP, KV
C     .. Arrays in Common ..
      DOUBLE PRECISION  RV(LV)
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      INTEGER           IA, IB, J
      LOGICAL           INIT
C     .. External Subroutines ..
      EXTERNAL          G05CAY, G05CAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /BG05CA/NORMAL, GAMMA
      COMMON            /CG05CA/RV, KV
      COMMON            /DG05CA/VNORML
C     .. Save statement ..
      SAVE              /AG05CA/, /BG05CA/, /CG05CA/, /DG05CA/
C     .. Executable Statements ..
C
C     Call G05CAZ to ensure that the initializations about to be
C     performed by G05CBF will not be overwritten by a subsequent
C     call t0 G05CAZ from G05CAY
C
      CALL G05CAZ(INIT)
C
C     Store 2*IABS(I)+1 in B(0,1:ILIM) as required by G05CAY
C
      IA = ABS(I)
      IB = IA/IPW2A
      B(0,1) = 4*(IA-IPW2A*IB) + 2
      DO 20 J = 2, ILIM
         IA = IB
         IB = IA/IPW2
         B(0,J) = IA - IPW2*IB
   20 CONTINUE
C
C     Re-initialize NORMAL, GAMMA and VNORML
C
      NORMAL = 1.0D0
      GAMMA = -1.0D0
      VNORML = 256.0D0
C
C     Re-initialize RV and set KV to skip first number
C
      CALL G05CAY(.TRUE.)
      KV = KV + 1
C
      RETURN
      END
