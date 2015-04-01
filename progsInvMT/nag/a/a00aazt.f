      SUBROUTINE A00AAZ(MSG)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Returns information about the particular implementation of the
C     NAG Fortran Library in use.
C
C     **************************************************************
C
C     Implementors must insert the correct details for each
C     distinct implementation.
C
C     **************************************************************
C
C     .. Array Arguments ..
      CHARACTER*80      MSG(20)
C     .. Executable Statements ..
      MSG(1) = ' *** Start of NAG Library implementation details ***'
      MSG(2) = ' '
      MSG(3) = ' Implementation title: Sun Solaris / 64-bit'
      MSG(4) = '            Precision: Fortran Double Precision'
      MSG(5) = '         Product Code: FLSO617DA'
      MSG(6) = '                 Mark: 17A'
      MSG(7) = ' '
      MSG(8) = ' Created using:'
      MSG(9) = '     hardware -   Sun Ultra Enterprise 2'
      MSG(10) = '     op. sys. -   Solaris 2.5.1'
      MSG(11) = '     compiler -   Sun f77 4.0'
      MSG(12) = ' '
      MSG(13) = ' '
      MSG(14) = ' '
      MSG(15) = ' Applicable to:'
      MSG(16) = '     hardware -   all SPARC Ultra systems'
      MSG(17) = '     op. sys. -   Solaris 2.5 or later'
      MSG(18) = '     compiler -   Sun f77 4.0 or later'
      MSG(19) = ' '
      MSG(20) = ' *** End of NAG Library implementation details ***'
      RETURN
      END
