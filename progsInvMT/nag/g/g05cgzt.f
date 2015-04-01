      SUBROUTINE G05CGZ(IA,XB,IERR)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     called by G05CGF to restore the current values of the internal
C     variables used by the random-number-generating routines
C
C     The contents of common blocks /AG05CA/, /BG05CA/ and /DG05CA/
C     must be restored from the arrays IA and XB as follows (not all
C     the contents of /AG05CA/ are restored):
C
C     The 59-bit integer N which is held with its most significant
C     11 bits in IA(1) and with its next most significant 12 bits
C     in turn in each of IA(2), IA(3), IA(4) and IA(5) must be
C     stored in B(0,1:ILIM) in /AG05CA/ (see text of G05CAY for
C     details).
C
C     The integer OPTION which may have up to 30 bits and is held
C     with its most significant 6 bits in IA(6) and with its next
C     most significant 12 bits in turn in IA(7) and IA(8) must be
C     stored in /AG05CA/.
C
C     (N.B. it is assumed that an integer can store at least 12
C     bits.)
C
C     The //real// variable VNORML which is held in XB(1) must be
C     stored in /DG05CA/.
C
C     The //real// variable NORMAL which is held in XB(2) must be
C     stored in /BG05CA/.
C
C     The //real// variable GAMMA which is held in XB(3) must be
C     stored in /BG05CA/.
C
C     ******************** ADVICE FOR IMPLEMENTORS *********************
C
C     This version of G05CGZ must be used in conjunction with the
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
C     G05CGZ must be implemented consistently with G05CAY.
C
C     If G05CAY has been implemented simply by selecting suitable
C     variant code according to the value of ILIM, then a consistent
C     implementation of G05CGZ may be obtained by using the variant
C     code supplied in comments beginning CAnn where the digits nn
C     are the value of ILIM.
C
C     If G05CAY has been implemented in machine code, it will still
C     be possible on many machines to implement G05CGZ in Fortran
C     and this will be satisfactory since it is not important for
C     G05CGZ to be particularly efficient. Essentially the code for
C     G05CGZ depends only on how the internal variable N is stored in
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
C     .. Scalar Arguments ..
      INTEGER           IERR
C     .. Array Arguments ..
      DOUBLE PRECISION  XB(3)
      INTEGER           IA(8)
C     .. Scalars in Common ..
      DOUBLE PRECISION  GAMMA, NORMAL, VNORML
      INTEGER           DEFOPT, OPTION, POSSOP
C     .. Arrays in Common ..
      INTEGER           B(0:LV,ILIM)
C     .. Local Scalars ..
      INTEGER           I, IB, IS, IT, J, JBITS, JPW2
      LOGICAL           INIT
C     .. External Subroutines ..
      EXTERNAL          G05CAY, G05CAZ
C     .. Common blocks ..
      COMMON            /AG05CA/B, OPTION, POSSOP, DEFOPT
      COMMON            /BG05CA/NORMAL, GAMMA
      COMMON            /DG05CA/VNORML
C     .. Save statement ..
      SAVE              /AG05CA/, /BG05CA/, /DG05CA/
C     .. Executable Statements ..
C
C     Call G05CAZ to ensure that the initializations about to be
C     performed by G05CGZ will not be overwritten by a subsequent
C     call to G05CAZ from G05CAY
C
      CALL G05CAZ(INIT)
C
C     Restore N from IA(1), . . . , IA(5).
C
      I = 1
      J = ILIM
      IT = 0
      IB = IA(1)
      JBITS = 11 - IBITS
      JPW2 = 2048
      GO TO 80
   20 JBITS = JBITS + 12
      IF (JBITS.GE.0) GO TO 40
      JPW2 = JPW2/4096
      GO TO 120
   40 JPW2 = 4096/JPW2
   60 IS = IB/JPW2
      B(0,J) = IT + IS
      J = J - 1
      IF (J.EQ.0) GO TO 140
      IT = 0
      IB = IB - IS*JPW2
      JBITS = JBITS - IBITS
   80 IF (JBITS.LT.0) GO TO 100
      JPW2 = JPW2/IPW2
      GO TO 60
  100 JPW2 = IPW2/JPW2
  120 IT = IT + IB*JPW2
      I = I + 1
      IB = IA(I)
      GO TO 20
C
C     Restore OPTION from IA(6), IA(7) and IA(8).
C
  140 OPTION = (IA(6)*64+IA(7))*4096 + IA(8)
C
C     Restore VNORML, NORMAL and GAMMA from XB
C
      VNORML = XB(1)
      NORMAL = XB(2)
      GAMMA = XB(3)
C
C     Re-initialize the buffer
C
      CALL G05CAY(.TRUE.)
C
      IERR = 0
      RETURN
      END
