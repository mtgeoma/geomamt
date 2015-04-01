      SUBROUTINE G05CFZ(IA,XB)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     called by G05CFF to save the current values of the internal
C     variables used by the random-number-generating routines
C
C     The contents of common blocks /AG05CA/, /BG05CA/ and /D05CA/ must
C     be saved in the arrays IA and XB as follows (not all the contents
C     of /AG05CA/ are saved):
C
C     The 59-bit integer N which is held in B(KV,1:ILIM) in /AG05CA/
C     must be stored in the array IA with the most significant 11
C     bits in IA(1) and with the next most significant 12 bits
C     in turn stored in each of IA(2), IA(3), IA(4) and IA(5).
C
C     The integer OPTION which may have up to 30 bits and is also
C     held in /AG05CA/ must be stored with the most significant
C     6 bits (of the notional 30 bits) in IA(6) and with the next
C     most significant 12 bits in turn stored in each of IA(7)
C     and IA(8).
C
C     (N.B. it is assumed that an integer can store at least 12
C     bits.)
C
C     The //real// variable VNORML which is held in /DG05CA/ must be
C     stored in XB(1).
C
C     The //real// variable NORMAL which is held in /BG05CA/ must be
C     stored in XB(2).
C
C     The //real// variable GAMMA which is held in /BG05CA/ must be
C     stored in XB(3).
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     This version of G05CFZ must be used in conjunction with the
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
C     G05CFZ must be implemented consistently with G05CAY.
C
C     If G05CAY has been implemented simply by selecting suitable
C     variant code according to the value of ILIM, then a consistent
C     implementation of G05CFZ may be obtained by using the variant
C     code supplied in comments beginning CAnn where the digits nn
C     are the value of ILIM.
C
C     If G05CAY has been implemented in machine code, it will still
C     be possible on many machines to implement G05CFZ in Fortran
C     and this will be satisfactory since it is not important for
C     G05CFZ to be particularly efficient. Essentially the code for
C     G05CFZ depends only on how the internal variable N is stored in
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
      INTEGER           IBITS, IPW2
      PARAMETER         (IBITS=60/ILIM,IPW2=2**IBITS)
C     .. Array Arguments ..
      DOUBLE PRECISION  XB(3)
      INTEGER           IA(8)
C     .. Scalars in Common ..
      DOUBLE PRECISION  GAMMA, NORMAL, VNORML
      INTEGER           DEFOPT, OPTION, POSSOP, KV
C     .. Arrays in Common ..
      DOUBLE PRECISION  RV(LV)
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      INTEGER           I, IB, IS, IT, J, JBITS, JPW2
      LOGICAL           INIT
C     .. External Subroutines ..
      EXTERNAL          G05CAZ
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /BG05CA/NORMAL, GAMMA
      COMMON            /CG05CA/RV, KV
      COMMON            /DG05CA/VNORML
C     .. Save statement ..
      SAVE              /AG05CA/, /BG05CA/, /CG05CA/, /DG05CA/
C     .. Executable Statements ..
C
C     Call G05CAZ to ensure that the contents of /AG05CA/ and
C     /BG05CA/ have been initialized
C
      CALL G05CAZ(INIT)
C
C     Store current value of N in IA(1), . . . , IA(5).
C
      I = 0
      J = ILIM
      IT = 0
      IB = B(KV,ILIM)
      JBITS = -11
      JPW2 = 2048
   20 JBITS = JBITS + IBITS
      IF (JBITS.GE.0) GO TO 40
      JPW2 = JPW2/IPW2
      GO TO 100
   40 JPW2 = IPW2/JPW2
   60 IS = IB/JPW2
      I = I + 1
      IA(I) = IT + IS
      IF (I.EQ.5) GO TO 120
      IT = 0
      IB = IB - IS*JPW2
      JBITS = JBITS - 12
      IF (JBITS.LT.0) GO TO 80
      JPW2 = JPW2/4096
      GO TO 60
   80 JPW2 = 4096/JPW2
  100 IT = IT + IB*JPW2
      J = J - 1
      IB = B(KV,J)
      GO TO 20
C
C     Store OPTION in IA(6), IA(7) and IA(8).
C
  120 IB = OPTION/4096
      IA(8) = OPTION - 4096*IB
      IS = IB/4096
      IA(7) = IB - 4096*IS
      IB = IS/4096
      IA(6) = IS - 4096*IB
C
C     Store VNORML, NORMAL and GAMMA in XB
C
      XB(1) = VNORML
      XB(2) = NORMAL
      XB(3) = GAMMA
C
      RETURN
      END
