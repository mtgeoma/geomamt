C     ==================================================================
                             PROGRAM D2ANIZ
C                           **************** 
C     ==================================================================
C     Modelling of MT fields in 2-D structures with anisotropy.
C     Version runs all in RAM!
C
C	Corrected <CORRECTED> on November 6, 1999 by J.Pek
C     ------------------------------------------------------------------
C     For the theoretical and algorithmic background see the paper:
C       Pek, J. and Verner, T., 1997: Finite difference modelling
C       of magnetotelluric fields in 2-D anisotropic media, 
C       Geophys. J. Int., 128, 505-521.
C     ==================================================================
C_D   DATA INPUT
C_D   ^^^^^^^^^^
C_D   FDAT (format A40) = path and name of the file containing input
C_D                       data - given from the KEYBOARD after the
C_D                       prompt 'GIVE DATA FILE [.DAT] ===>' appears.
C_D                       The extension .DAT must be typed in explicit-
C_D                       ely, it is not concatenated automatically to
C_D                       the file name!
C_D   FRES (format A40) = path and name for the result file into which
C_D                       all the results are directed - given from the
C_D                       KEYBOARD after the prompt 'GIVE OUTPUT FILE
C_D                       [.RES] =>' appears. The extension .RES is
C_D                       only a proposal, it is not added automatically
C_D                       to the given file name!
C_D
C_D   All input data which specify the periods, mesh parameters and
C_D   the structure involved are given in the FDAT file which has the
C_D   following structure:
C_D   NPER (format I5) = number of periods used for modelling, maximum
C_D                      number of periods is set to NPMAX=30, but in
C_D                      principle can be changed to any integer without
C_D                      affecting the size of the executable
C_D   PER(I),I=1,NPER (format 5F10.5) = periods of the field in sec
C_D   N,M,MD (format 3I5) = number of FD mesh lines in horizontal
C_D                       direction (N), in vertical direction (M),
C_D                       and within the air layer (MD) including the
C_D                       earth-air interface line. Maximum mesh size
C_D                       in this version is NMAX=151 for N and MMAX=59
C_D                       for M. Increasing NMAX and MMAX causes the
C_D                       size of the exe to increase quite rapidly,
C_D                       particularly when MMAX (vertical mesh extent)
C_D                       is changed!
C_D   SY(J),J=1,N1 (format 10F5.2) = horizontal FD mesh steps in km,
C_D                       N1 is for N-1
C_D   SZ(K),K=1,M1 (format 10F5.2) = vertical FD mesh steps in km, M1
C_D                       is for M-1
C_D   ..................................................................
C_D   Next the conductivity map of the structure follows in this form:
C_D
C_D         ÚÄ> 1.column
C_D         <--------- (N-1) characters --------->
C_D         11111111111111.....1111111111111111111<ÄÄÄ¿
C_D         11111111111111.....1111111111111111111    ³
C_D          :    :    :    :    :    :    :    :   (M-1) lines
C_D         22222222233333.....3333344444444444444    ³
C_D         44444444444444.....5555555555555555555<ÄÄÄÙ
C_D
C_D   The format for each line of the conductivity map is A80, positions
C_D   exceeding (N-1) are automatically filled with blanks.
C_D   Each position within this map characterizes the conductivity
C_D   within the corresponding FD mesh cell using a character. The
C_D   symbols within the cells are mapped onto a succession of
C_D   integers in the following way (not too convenient, should be
C_D   changed so that the user can give any ASCII symbol without
C_D   observing the ASCII sequence of the symbols, but this change
C_D   have not been made yet):
C_D
C_D   1  2  3  4  5  6  7  8  9  :  ;  <  =  >  ?  @  A  B  C  D
C_D   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 =NCMAX
C_D   ..................................................................
C_D   NRO (format I5) = number of resistivities within the structure,
C_D                     maximum NCMAX=20 in this version
C_D   ÚÄÄÄÄÄÄÄÄÄ loop over I=1,NRO ÄÄÄÄÄÄÄ¿
C_D   ³                                   V
C_D   ³ LAN(I),RON(I,1),RON(I,2),RON(I,3),USTR(I),UDIP(I),U3(I)
C_D   ³   (format I5,6F10.2) = resistivity of the I-th domain of the
C_D   ³                 conductivity map. LAN(I) is 0 for an isotropic
C_D   ³                 domain and 1 for an anisotropic one. The conduc-
C_D   ³                 tivity tensor for the particular domain is com-
C_D   ³                 puted from its principal values åx=1./RON(I,1),
C_D   ³                 åy=1./RON(I,2), åz=1./RON(I,3) by rotating this
C_D   ³                 diagonal matrix by USTR(I) around z-axis (stri-
C_D   ³                 ke), next by UDIP(I) around the new x-axis (dip)
C_D   ³                 and, finally, by U3(I) around the latest z-axis.
C_D   ³                                   V
C_D   ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
C_D   END OF DATA INPUT
C     ^^^^^^^^^^^^^^^^^
C
C     This version of the program does not need the RAM disk to be de-
C     fined. The program must be, however, compiled with a Fortran   
C     compiler which allows the whole available RAM to be addressed.
C     This program has been successfully compiled with LAHEY-Fortran
C     compiler F77L-EM/32, for the current setting of parameters
C     NMAX=151, MMAX=59 the size of the resulting executable code   
C     was about 17.2 MB (single precision version 8.7 MB).
C     In general, the memory required for running the code is given
C     by the size of the finite difference matrix APOM, which really
C     is 
C                size(APOM)=(N-2)*(2*M-ND-2)*(2*M-ND)
C     of complex*16 numbers. So, the parameters NMAX,MMAX can be
C     altered to fit both the model and the computer capacities.
C
C     In comparison with earlier versions a double precision 
C     version of the subroutine EM1ANZ is delivered with this code,
C     and also the Gaussian elimination is performed in double
C     precision arithmetics (it increases the accuracy, but also
C     causes the code to be quite a fatty).
C
C     Some Fortran compilers do not like the C-convention formats
C     in I/O operations (like write(*,'(a\)') with the backslash).
C     Then, do remove them, it only causes the prompt and the input
C     cursor to appear on different lines.
C
C     We tried to make the structure of the output file as clear as
C     possible, but... In principle, the output file repeats the
C     input data at the beginning. The results are given in a form
C     of clusters for each period and each point on the surface.
C     The symbols for the MT/MV functions are:
C     EX, EY, HX, HY, HZ - components of the elmg field, normalized
C                          relatively to the boundary conditions
C     EN, HN - leftmost (normal) boundary value of the respective
C              field component (used to output relative field values
C              on the surface of the model)
C     ZXX, ZXY, ZYX, ZYY - impedances in practical units (mV/km)/nT
C     ROAXY, ROAYX - apparent resistivities XY (E-polarization in
C                    the degenerate isotropic case) and YX (H-pola-
C                    rization in the degenerate isotropic case)
C     PH(ZXY), PH(ZYX) - phase of the respective impedance
C     SWIFT - Swift's principal direction 
C     ANIZ - impedance anisotropy, always > 1
C     SKEW - classical MT skew parameter
C     WX, WY - geomagnetic transfer functions, defined by the
C              linear relation HZ = WX * HX + WY * HY
C     TR, TI - modules of the real and imaginary geomagnetic transfer
C              functions, defined as 
C              Real Transfer Function = [ REAL(WX), REAL(WY) ]
C              Imag Transfer Function = [ AIMAG(WX), AIMAG(WY) ]
C     PH(TR), PH(TI) - azimuth of the real and imaginary geomagnetic
C                      transfer functions from the X axis
C     X is the structural strike of the 2-D structure!
C     The output is formatted to fit the 132 characters printer, so
C     a piece of the results is hidden beyond the margins of the
C     80 characters screen!
C
C     DISABLED!!!
C     Running the program always produces an auxiliary file DARRT.RES
C     with MT functions along the profile with a simple formatting.
C     It was used to peck out the values for subsequent graphic
C     processing. If no use can be made of this file, ignore it!
C
C     Notice, that a lot of exceptional situations is not treated in
C     this program. E.g. if too short periods are considered and the
C     depth of the mesh is too large, so that the field disappears
C     long before it reaches the bottom of the model, the program will
C     crash with an overflow error caused by the boundary conditions
C     routine. It is only an example of a bug we know about, but there
C     will quite certainly be many more, which we would like to be
C     be informed about - please use either VCV@IG.CAS.CZ or
C     JPK@IG.CAS.CZ, i.e. Vaclav Cerv or Josef Pek, to contact us.
C
C-----------------------------------------------------------------------
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20,NPMAX=30,NBMAX=10)
      PARAMETER(NEHMAX=2*(NMAX-2)*(MMAX-2),MEHMAX=2*(MMAX-2)+3)
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(OMMI0=8.E-7*PI*PI)
      COMPLEX BOULE,BOURE,BOUUE,BOUDE,BOULH,BOURH,BOUUH,BOUDH
      COMPLEX BEXL,BEXR,BEYL,BEYR,BHXL,BHXR,BHYL,BHYR
      COMPLEX*16 APOM,BPOM
      COMPLEX SEX,SEY,SHX,SHY,SHZ,SEZ
      COMPLEX AHES
      DIMENSION PER(NPMAX)
      CHARACTER FDAT*40,FRES*40
C
      COMMON /NHS/ NEH
      COMMON /RHS/ BPOM(NEHMAX)
      COMMON /MAT/ APOM(NEHMAX,MEHMAX)
      COMMON /CAP/ OMMI,PERIOD
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /MED/ IRG(NMAX-1,MMAX-1),SG(NCMAX,3,3)
      COMMON /DIS/ SYY(NMAX),SZZ(MMAX)
      COMMON /BOM/ NLBL,HBL(NBMAX),SGBL(NBMAX,3,3),
     &             NLBR,HBR(NBMAX),SGBR(NBMAX,3,3)
      COMMON /BUE/ BOULE(MMAX),BOURE(MMAX),BOUUE(NMAX),BOUDE(NMAX)
      COMMON /BUH/ BOULH(MMAX),BOURH(MMAX),BOUUH(NMAX),BOUDH(NMAX)
      COMMON /BOS/ BEXL(2),BEYL(2),BHXL(2),BHYL(2),
     /             BEXR(2),BEYR(2),BHXR(2),BHYR(2)
      COMMON /SFD/ SEX(2,NMAX),SEY(2,NMAX,2),SEZ(2,NMAX,2),
     &             SHX(2,NMAX),SHY(2,NMAX),SHZ(2,NMAX)
      COMMON /SCO/ AHES(NMAX-2)
C-----------------------------------------------------------------------
C
    2 WRITE(*,'(A)')' GIVE DATA FILE [.DAT] ===> '
      READ(*,'(A)')FDAT
      OPEN(10,FILE=FDAT,STATUS='OLD',ERR=2)
C
    3 WRITE(*,'(A)')' GIVE OUTPUT FILE [.RES] => '
      READ(*,'(A)')FRES
      OPEN(11,FILE=FRES,STATUS='NEW',ERR=3)
C
      WRITE(11,2001)
c	open(11,file='prdel1.dat')
C
      READ(10,1000)NPER
      READ(10,1001)(PER(I),I=1,NPER)
C
      CALL RE2ANZ
      CALL BOUMED
c      open(7,file='darrt.res')
C
      DO 1 I=1,NPER
      PERIOD=PER(I)
      OMMI=OMMI0/PERIOD
      NEH=N2*(2*M2-MD+1)
      MEH=2*M2-MD+4
C
      IHPOL=1
      CALL BOUFLD(IHPOL)
      CALL KOEF3
      WRITE(11,2000)PERIOD,IHPOL
      CALL GAUSSR
      CALL SURFLD(IHPOL)
C
      IHPOL=2
      CALL BOUFLD(IHPOL)
      CALL KOEF3
      WRITE(11,2000)PERIOD,IHPOL
      CALL GAUSSR
      CALL SURFLD(IHPOL)
C
      CALL SURFCE
C
    1 CONTINUE
C
        write(*,*)char(7)
      STOP
C
 1000 FORMAT(I5)
 1001 FORMAT(5F10.4)
 1002 FORMAT(I5,4F10.2)
 2000 FORMAT(1H1,9X,'ITERATIVE SOLUTION FOR THE PERIOD ',F10.2,',POLARIZ
     /ATION ',I2/10X,'==================================================
     /========================================')
 2001 FORMAT(1H1,27X,'DIRECT MAGNETOTELLURIC MODELLING IN TWO-DIMENSIONA
     /L ANISOTROPIC STRUCTURES'/27X,'===================================
     /========================================='//)
C
      END
      SUBROUTINE RE2ANZ
C     ==================================================================
C
C-----------------------------------------------------------------------
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20)
      DIMENSION LAN(NCMAX),RON(NCMAX,3),SGN(NCMAX,3),
     &          USTR(NCMAX),UDIP(NCMAX),U3(NCMAX)
      character rad80*150
C
      COMMON /CAP/ OMMI,PERIOD
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /MED/ IRG(NMAX-1,MMAX-1),SG(NCMAX,3,3)
      COMMON /DIS/ SYY(NMAX),SZZ(MMAX)
C-----------------------------------------------------------------------
C
      READ(10,1000)N,M,MD
 1000 FORMAT(3I5)
      N1=N-1
      N2=N-2
      M1=M-1
      M2=M-2
      NE=N2*M2
      MH=M-MD+1
      MH1=MH-1
      MH2=MH-2
      NH=N2*MH2
      MD1=MD-1
      MDP1=MD+1
      NEH=NE+NH
C
      READ(10,1001)(SY(J),J=1,N1)
      READ(10,1001)(SZ(K),K=1,M1)
 1001 FORMAT(10F7.2)
      SYY(1)=0.
      DO 20 J=2,N
      SYY(J)=SYY(J-1)+SY(J-1)
   20 CONTINUE
      SZZ(MD)=0.
      DO 21 K=MD1,1,-1
      SZZ(K)=SZZ(K+1)-SZ(K)
   21 CONTINUE
      DO 22 K=MDP1,M
      SZZ(K)=SZZ(K-1)+SZ(K-1)
   22 CONTINUE
C
      DO 1 K=1,M1
c-----------------------------------------------------------------------
c     Possible symbols for conductivity map:
c      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 =NCMAX
c      1  2  3  4  5  6  7  8  9  :  ;  <  =  >  ?  @  A  B  C  D
c-----------------------------------------------------------------------
      read(10,1005)rad80
 1005 format(a150)
      do 50 j=1,n1
      irg(j,k)=ichar(rad80(j:j))-48
   50 continue
c      READ(10,1002)(IRG(J,K),J=1,N1)
 1002 FORMAT(150I1)
    1 CONTINUE
C
      READ(10,1003)NRO
 1003 FORMAT(I5)
      DO 2 I=1,NRO
      READ(10,1004)LAN(I),RON(I,1),RON(I,2),RON(I,3),
     /            USTR(I),UDIP(I),U3(I)
 1004 FORMAT(I5,6F10.2)
C
      IF(LAN(I).EQ.1)GOTO 3
      RON(I,2)=RON(I,1)
      RON(I,3)=RON(I,1)
      IF(RON(I,1).LE.0.)GOTO 4
      SGN(I,1)=1./RON(I,1)
      GOTO 5
    4 SGN(I,1)=0.
    5 SGN(I,2)=SGN(I,1)
      SGN(I,3)=SGN(I,1)
      GOTO 6
C
    3 DO 8 J=1,3
      IF(RON(I,J).LE.0.)GOTO 7
      SGN(I,J)=1./RON(I,J)
      GOTO 8
    7 SGN(I,J)=0.
    8 CONTINUE
    6 RSTR=PI*USTR(I)/180.
      RDIP=PI*UDIP(I)/180.
      R3=PI*U3(I)/180.
      SPS=SIN(RSTR)
      CPS=COS(RSTR)
      STH=SIN(RDIP)
      CTH=COS(RDIP)
      SFI=SIN(R3)
      CFI=COS(R3)
      POM1=SGN(I,1)*CFI*CFI+SGN(I,2)*SFI*SFI
      POM2=SGN(I,1)*SFI*SFI+SGN(I,2)*CFI*CFI
      POM3=(SGN(I,1)-SGN(I,2))*SFI*CFI
      C2PS=CPS*CPS
      S2PS=SPS*SPS
      C2TH=CTH*CTH
      S2TH=STH*STH
      CSPS=CPS*SPS
      CSTH=CTH*STH
C
      SG(I,1,1)=POM1*C2PS+POM2*S2PS*C2TH-2.*POM3*CTH*CSPS
     /         +SGN(I,3)*S2TH*S2PS
      SG(I,1,2)=POM1*CSPS-POM2*C2TH*CSPS+POM3*CTH*(C2PS-S2PS)
     /         -SGN(I,3)*S2TH*CSPS
      SG(I,1,3)=-POM2*CSTH*SPS+POM3*STH*CPS+SGN(I,3)*CSTH*SPS
      SG(I,2,1)=SG(I,1,2)
      SG(I,2,2)=POM1*S2PS+POM2*C2PS*C2TH+2.*POM3*CTH*CSPS
     /         +SGN(I,3)*S2TH*C2PS
      SG(I,2,3)=POM2*CSTH*CPS+POM3*STH*SPS-SGN(I,3)*CSTH*CPS
      SG(I,3,1)=SG(I,1,3)
      SG(I,3,2)=SG(I,2,3)
      SG(I,3,3)=POM2*S2TH+SGN(I,3)*C2TH
C
    2 CONTINUE
C
      WRITE(11,2000)
      WRITE(11,2001)N1
      WRITE(11,2002)(SY(J),J=1,N1)
      WRITE(11,2003)MD1
      WRITE(11,2002)(SZ(K),K=1,MD1)
      WRITE(11,2004)MH1
      WRITE(11,2002)(SZ(K),K=MD,M1)
      WRITE(11,2005)
      DO 10 K=1,M1
      do 51 j=1,n1
      rad80(j:j)=char(irg(j,k)+48)
   51 continue
      if(n1.lt.150)then
        do 52 j=n,150
        rad80(j:j)=' '
   52   continue
      endif
      write(11,2013)rad80
c      WRITE(11,2006)(IRG(J,K),J=1,N1)
   10 CONTINUE
      WRITE(11,2007)NRO
      WRITE(11,2009)
      WRITE(11,2008)
      DO 11 I=2,NRO
      WAVEL=SQRT(10.*RON(I,1))
      SKIND=WAVEL/(2.*PI)
      WRITE(11,2010)I,LAN(I),RON(I,1),SGN(I,1),USTR(I),
     /             SG(I,1,1),SG(I,1,2),SG(I,1,3),WAVEL,SKIND
      WAVEL=SQRT(10.*RON(I,2))
      SKIND=WAVEL/(2.*PI)
      WRITE(11,2011)RON(I,2),SGN(I,2),UDIP(I),
     /             SG(I,2,1),SG(I,2,2),SG(I,2,3),WAVEL,SKIND
      WAVEL=SQRT(10.*RON(I,3))
      SKIND=WAVEL/(2.*PI)
      WRITE(11,2012)RON(I,3),SGN(I,3),U3(I),
     /             SG(I,3,1),SG(I,3,2),SG(I,3,3),WAVEL,SKIND
      WRITE(11,2008)
   11 CONTINUE
C
      RETURN
C
 2000 FORMAT(1H1,1X,'DESCRIPTION OF THE MODEL AND THE MESH PARAMETERS'/2
     /X,'------------------------------------------------'/)
 2001 FORMAT(/2X,'MESH SPACINGS IN KM IN HORIZONTAL DIRECTION (',I3,')')
 2002 FORMAT(2X,10F7.2)
 2003 FORMAT(/2X,'MESH SPACINGS IN KM IN VERTICAL DIRECTION IN THE AIR (
     /',I3,')')
 2004 FORMAT(/2X,'MESH SPACINGS IN KM IN VERTICAL DIRECTION IN THE EARTH
     / (',I3,')')
 2005 FORMAT(/2X,'CONDUCTIVITY MAP OF THE STRUCTURE'/2X,'123456789012345
     /678901234567890123456789012345678901234567890'/)
 2006 FORMAT(2X,60I1)
 2007 FORMAT(1H1,1X,'ELECTRICAL PARAMETERS OF THE SUBDOMAINS (',I3,')'/2
     /X,'===============================================================
     /================================================')
 2008 FORMAT(2X,'-------------------------------------------------------
     /--------------------------------------------------------')
 2009 FORMAT(2X,'INDEX',1X,'ANIZ',4X,'RES1',7X,'COND1',5X,'STRIKE',7X,'S
     /GXX',8X,'SGXY',8X,'SGXZ',6X,'WAVEL1',6X,'SKIND1'/3X,'===',10X,'RES
     /2',7X,'COND2',7X,'DIP',8X,'SGYX',8X,'SGYY',8X,'SGYZ',6X,'WAVEL2',6
     /X,'SKIND2'/16X,'RES3',7X,'COND3',5X,'ANGLE3',7X,'SGZX',8X,'SGZY',8
     /X,'SGZZ',6X,'WAVEL3',6X,'SKIND3')
 2010 FORMAT(3X,I2,4X,I2,2X,F10.2,2X,F10.6,2X,F6.1,4X,3(F10.6,2X),2X,2(F
     /10.2,2X))
 2011 FORMAT(3X,'===',7X,F10.2,2X,F10.6,2X,F6.1,4X,3(F10.6,2X),2X,2(F10.
     /2,2X))
 2012 FORMAT(13X,F10.2,2X,F10.6,2X,F6.1,4X,3(F10.6,2X),2X,2(F10.2,2X))
 2013 format(a150)
C
      END
      SUBROUTINE BOUMED
C     ==================================================================
C
C-----------------------------------------------------------------------
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20,NBMAX=10)
      COMMON /BOM/ NLBL,HBL(NBMAX),SGBL(NBMAX,3,3),
     &             NLBR,HBR(NBMAX),SGBR(NBMAX,3,3)
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /MED/ IRG(NMAX-1,MMAX-1),SG(NCMAX,3,3)
C-----------------------------------------------------------------------
C
      MDP1=MD+1
C
      NLBL=1
      HBL(NLBL)=0.
      JJ=IRG(1,MD)
      DO 1 J1=1,3
      DO 2 J2=1,3
      SGBL(NLBL,J1,J2)=SG(JJ,J1,J2)
    2 CONTINUE
    1 CONTINUE
      DO 3 I=MDP1,M1
      IDEL=IRG(1,I)-IRG(1,I-1)
      IF(IDEL)4,5,4
    5 HBL(NLBL)=HBL(NLBL)+SZ(I-1)
      GOTO 3
    4 HBL(NLBL)=HBL(NLBL)+SZ(I-1)
      NLBL=NLBL+1
      HBL(NLBL)=0.
      JJ=IRG(1,I)
      DO 6 J1=1,3
      DO 7 J2=1,3
      SGBL(NLBL,J1,J2)=SG(JJ,J1,J2)
    7 CONTINUE
    6 CONTINUE
    3 CONTINUE
C
      NLBR=1
      HBR(NLBR)=0.
      JJ=IRG(N1,MD)
      DO 8 J1=1,3
      DO 9 J2=1,3
      SGBR(NLBR,J1,J2)=SG(JJ,J1,J2)
    9 CONTINUE
    8 CONTINUE
      DO 10 I=MDP1,M1
      IDEL=IRG(N1,I)-IRG(N1,I-1)
      IF(IDEL)11,12,11
   12 HBR(NLBR)=HBR(NLBR)+SZ(I-1)
      GOTO 10
   11 HBR(NLBR)=HBR(NLBR)+SZ(I-1)
      NLBR=NLBR+1
      HBR(NLBR)=0.
      JJ=IRG(N1,I)
      DO 13 J1=1,3
      DO 14 J2=1,3
      SGBR(NLBR,J1,J2)=SG(JJ,J1,J2)
   14 CONTINUE
   13 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE BOUFLD(IHPOL)
C     ==================================================================
C
C-----------------------------------------------------------------------
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20,NBMAX=10)
      COMPLEX EX,EY,EZ,HX,HY,HX0,HY0,ZXX,ZXY,ZYX,ZYY
      COMPLEX BOULE,BOURE,BOUUE,BOUDE,BOULH,BOURH,BOUUH,BOUDH
      COMPLEX BEXL,BEYL,BHXL,BHYL,BEXR,BEYR,BHXR,BHYR
      DIMENSION EX(MMAX),EY(MMAX),EZ(MMAX),HX(MMAX),HY(MMAX)
C
      COMMON /BOM/ NLBL,HBL(NBMAX),SGBL(NBMAX,3,3),
     &             NLBR,HBR(NBMAX),SGBR(NBMAX,3,3)
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /CAP/ OMMI,PERIOD
      COMMON /BUE/ BOULE(MMAX),BOURE(MMAX),BOUUE(NMAX),BOUDE(NMAX)
      COMMON /BUH/ BOULH(MMAX),BOURH(MMAX),BOUUH(NMAX),BOUDH(NMAX)
      COMMON /DIS/ SYY(NMAX),SZZ(MMAX)
      COMMON /BOS/ BEXL(2),BEYL(2),BHXL(2),BHYL(2),
     /             BEXR(2),BEYR(2),BHXR(2),BHYR(2)
C-----------------------------------------------------------------------
C
      T=PERIOD
      MD1=MD-1
      MM=M-MD1
      IF(IHPOL.EQ.2)GOTO 1
      HX0=0.
      HY0=1.
c      hx0=cmplx(2.,1.)
c      hy0=cmplx(2.4,0.)
      GOTO 2
    1 CONTINUE
      HX0=-1.
      HY0=0.
    2 CONTINUE
C
      CALL EM1ANZ(NLBL,HBL,SGBL,T,M,SZZ,HX0,HY0,EX,EY,EZ,HX,HY,
     /            ZXX,ZXY,ZYX,ZYY,IE)
      DO 3 K=1,M
      BOULE(K)=EX(K)
    3 CONTINUE
      DO 4 K=MD,M
      KK=K-MD1
      BOULH(KK)=HX(K)
    4 CONTINUE
      BEXL(IHPOL)=EX(MD)
      BEYL(IHPOL)=EY(MD)
      BHXL(IHPOL)=HX(MD)
      BHYL(IHPOL)=HY(MD)
C
      CALL EM1ANZ(NLBR,HBR,SGBR,T,M,SZZ,HX0,HY0,EX,EY,EZ,HX,HY,
     /            ZXX,ZXY,ZYX,ZYY,IE)
      DO 5 K=1,M
      BOURE(K)=EX(K)
    5 CONTINUE
      DO 6 K=MD,M
      KK=K-MD1
      BOURH(KK)=HX(K)
    6 CONTINUE
      BEXR(IHPOL)=EX(MD)
      BEYR(IHPOL)=EY(MD)
      BHXR(IHPOL)=HX(MD)
      BHYR(IHPOL)=HY(MD)
C
      DO 7 J=1,N
      DIRE=SYY(J)/SYY(N)
      BOUUE(J)=BOULE(1)+DIRE*(BOURE(1)-BOULE(1))
      BOUUH(J)=BOULH(1)+DIRE*(BOURH(1)-BOULH(1))
      BOUDE(J)=BOULE(M)+DIRE*(BOURE(M)-BOULE(M))
      BOUDH(J)=BOULH(MM)+DIRE*(BOURH(MM)-BOULH(MM))
    7 CONTINUE
C
      RETURN
      END
      SUBROUTINE EM1ANZ(NL,H,SG,T,NZS,ZS,HX0,HY0,EX,EY,EZ,HX,HY,
     /                  ZXX,ZXY,ZYX,ZYY,IE)
C     ==================================================================
C
C     COMPUTES 1D ELECTROMAGNETIC FIELD AND SURFACE IMPEDANCE TENSOR
C     FOR AN ANIZOTROPIC LAYERED MEDIUM EXCITED BY A HOMOGENEOUS PLANE
C     ELECTROMAGNETIC WAVE WITH THE PERIOD T.
C
C     INPUT  = NL - NUMBER OF LAYERS (NL.LE.10)
C     *****    H  - ARRAY OF NL VALUES OF THE THICKNESSES OF THE
C                   LAYERS IN KM, H(NL) MAY BE AN ARBITRARY REAL
C              SG - ARRAY OF NL*3*3 VALUES OF THE COMPONENTS OF THE
C                   CONDUCTIVITY TENSOR IN S/M, ORGANIZED AS FOLLOWS
C                   SG(IL,1,1),SG(IL,1,2),SG(IL,1,3)=SGXX,SGXY,SGXZ
C                   SG(IL,2,1),SG(IL,2,2),SG(IL,2,3)=SGYX,SGYY,SGYZ
C                   SG(IL,3,1),SG(IL,3,2),SG(IL,3,3)=SGZX,SGZY,SGZZ
C              T - PERIOD OF THE ELMG FIELD IN SEC
C              NZS - NUMBER OF DEPTH AND HIGHT LEVELS AT WHICH THE
C                    ELMG FIELD IS TO BE COMPUTED (NZS.LE.60)
C              ZS - ARRAY OF NZS VALUES OF VERTICAL COORDINATES OF THE
C                   DEPTH LEVELS (WHEN ZS(IZ).GE.0.) OR HIGHT LEVELS
C                   (WHEN ZS(IZ).LT.0.) IN KM AT WHICH THE COMPONENTS
C                   OF THE ELMG FIELD ARE TO BE COMPUTED
C              HX0,HY0 - VALUES OF THE NORMALIZATION CONSTANTS FOR THE
C                        MAGNETIC FIELD IN X AND Y DIRECTION AT THE
C                        SURFACE (IN A/M)
C     OUTPUT = EX,EY - ARRAYS OF NZS VALUES OF THE ELECTRIC COMPONENTS
C     ******           AT THE INDIVIDUAL LEVELS ZS(IZ) (IN V/M)
C              HX,HY - ARRAYS OF NZS VALUES OF THE MAGNETIC COMPONENTS
C                      AT THE INDIVIDUAL LEVELS ZS(IZ) (IN A/M)
C              EZ - ARRAY OF NZS VALUES OF THE VERTICAL ELECTRIC
C                   COMPONENT AT THE INDIVIDUAL LEVELS ZS(IZ) (IN A/M)
C              ZXX,ZXY,ZYX,ZYY - COMPONENTS OF THE IMPEDANCE TENSOR
C                                AT THE SURFACE OF THE EARTH (IN OHMS)
C              IE - ERROR INDICATION, IF IE.EQ.1 THEN THE CONDUCTIVITY
C                   TENSOR IS NON-SYMMETRIC AND NO COMPUTATIONS ARE
C                   CARRIED OUT. IF IE.EQ.0 THE COMPUTATIONS HAVE BEEN
C                   PERFORMED WITHOUT ERROR
C
C     ------------------------------------------------------------------
C
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(NBMAX=10,MMAX=59)
      COMPLEX*16 SOUC,WN1,WN2,WN1D,WN2D,WT1,WT2,Q1D,Q2ID,Q1,Q2I,Q1DQ2,DQ
      COMPLEX*16 S,S1,S2,C1,C2,A11,A12,A21,A22,B11,B12,B21,B22,DETA
      COMPLEX*16 ZXXD,ZXYD,ZYXD,ZYYD,T320,T340,T420,T440,GAMD,GAM1,GAM2
      COMPLEX*16 PB1,PB2,PB3,PB4,PL1,PL2,PL3,PL4,F1,F2,EXD,EYD,HXD,HYD
      COMPLEX*16 EZD
      COMPLEX*16 SH,CH,X
      complex*8 ex,ey,hx,hy,ez,hx0,hy0,zxx,zxy,zyx,zyy
C
      DIMENSION H(NBMAX),SG(NBMAX,3,3),ZS(MMAX)
      DIMENSION EXD(MMAX),EYD(MMAX),HXD(MMAX),HYD(MMAX),EZD(MMAX)
      dimension ex(mmax),ey(mmax),hx(mmax),hy(mmax),ez(mmax)
      DIMENSION WN1D(NBMAX),WN2D(NBMAX),Q1D(NBMAX),Q2ID(NBMAX)
      DIMENSION S(4,4),SOUC(NBMAX,4,4)
C
C     COMPLEX HYPERBOLIC FUNCTIONS DEFINED
C     ------------------------------------
      SH(X)=0.5*(CDEXP(X)-CDEXP(-X))
      CH(X)=0.5*(CDEXP(X)+CDEXP(-X))
C
C     FUNDAMENTAL CONSTANTS AND COMMON PARAMETERS
C     -------------------------------------------
      OMMI=8.E-7*PI*PI/T
      SOMMI2=2.E-3*PI/SQRT(10.*T)
C
C     CHARACTERISTICS OF THE BASEMENT HALFSPACE
C     -----------------------------------------
      L=NL
      AXX=SG(L,1,1)-SG(L,1,3)*SG(L,3,1)/SG(L,3,3)
      AXY=SG(L,1,2)-SG(L,1,3)*SG(L,3,2)/SG(L,3,3)
      AYX=SG(L,2,1)-SG(L,2,3)*SG(L,3,1)/SG(L,3,3)
      AYY=SG(L,2,2)-SG(L,2,3)*SG(L,3,2)/SG(L,3,3)
      IF(AXY.NE.AYX)GOTO 100
      ADA=AXX+AYY
      ADD=AXX-AYY
      ANP=AXY*AYX
      WNZ=SQRT(ADD*ADD+4.*ANP)
      IF(AXX.LT.AYY)WNZ=-WNZ
      WNZ1=SQRT(0.5*(ADA+WNZ))
      WNZ2=SQRT(0.5*(ADA-WNZ))
      WN1=(1.-IC)*SOMMI2*WNZ1
      WN2=(1.-IC)*SOMMI2*WNZ2
      WN1D(L)=IC*WN1/OMMI
      WN2D(L)=IC*WN2/OMMI
      IF(AYX.EQ.0.)GOTO 40
      Q1D(L)=2.*AYX/(ADD+WNZ)
      Q2ID(L)=0.5*(ADD-WNZ)/AYX
      Q1=IC*Q1D(L)/OMMI
      Q2I=-IC*Q2ID(L)*OMMI
      GOTO 41
   40 Q1D(L)=0.
      Q2ID(L)=0.
      Q1=0.
      Q2I=0.
   41 Q1DQ2=Q1D(L)*Q2ID(L)
      DQ=1.-Q1DQ2
      DO 1 I=1,4
      DO 2 J=1,4
      SOUC(L,I,J)=0.
    2 CONTINUE
      SOUC(L,I,I)=1.
    1 CONTINUE
C
C     IMPEDANCE TENSOR AT THE TOP OF THE BASEMENT
C     -------------------------------------------
      ZXXD=(1./WN1D(L)-1./WN2D(L))*Q2ID(L)/DQ
      ZXYD=(1./WN1D(L)-Q1DQ2/WN2D(L))/DQ
      ZYXD=-(1./WN2D(L)-Q1DQ2/WN1D(L))/DQ
      ZYYD=(1./WN1D(L)-1./WN2D(L))*Q1D(L)/DQ
      zxx=zxxd
      zxy=zxyd
      zyx=zyxd
      zyy=zyyd
C
      IF(NL.EQ.1)GOTO 3
      NL1=NL-1
C
C     CHARACTERISTICS OF THE INDIVIDUAL LAYERS
C     ----------------------------------------
      DO 4 L=NL1,1,-1
      AXX=SG(L,1,1)-SG(L,1,3)*SG(L,3,1)/SG(L,3,3)
      AXY=SG(L,1,2)-SG(L,1,3)*SG(L,3,2)/SG(L,3,3)
      AYX=SG(L,2,1)-SG(L,2,3)*SG(L,3,1)/SG(L,3,3)
      AYY=SG(L,2,2)-SG(L,2,3)*SG(L,3,2)/SG(L,3,3)
      IF(AXY.NE.AYX)GOTO 100
      ADA=AXX+AYY
      ADD=AXX-AYY
      ANP=AXY*AYX
      WNZ=SQRT(ADD*ADD+4.*ANP)
      IF(AXX.LT.AYY)WNZ=-WNZ
      WNZ1=SQRT(0.5*(ADA+WNZ))
      WNZ2=SQRT(0.5*(ADA-WNZ))
      WN1=(1.-IC)*SOMMI2*WNZ1
      WN2=(1.-IC)*SOMMI2*WNZ2
      WN1D(L)=IC*WN1/OMMI
      WN2D(L)=IC*WN2/OMMI
      WT1=1.E+3*WN1*H(L)
      WT2=1.E+3*WN2*H(L)
      IF(AYX.EQ.0.)GOTO 5
      Q1D(L)=2.*AYX/(ADD+WNZ)
      Q2ID(L)=0.5*(ADD-WNZ)/AYX
      Q1=IC*Q1D(L)/OMMI
      Q2I=-IC*Q2ID(L)*OMMI
      GOTO 6
    5 Q1D(L)=0.
      Q2ID(L)=0.
      Q1=0.
      Q2I=0.
    6 Q1DQ2=Q1D(L)*Q2ID(L)
      DQ=1.-Q1DQ2
      S1=SH(WT1)
      S2=SH(WT2)
      C1=CH(WT1)
      C2=CH(WT2)
      S(1,1)=C1-Q1DQ2*C2
      S(1,2)=-Q2ID(L)*(C1-C2)
      S(1,3)=Q2ID(L)*(S1/WN1D(L)-S2/WN2D(L))
      S(1,4)=S1/WN1D(L)-Q1DQ2*S2/WN2D(L)
      S(2,1)=Q1D(L)*(C1-C2)
      S(2,2)=-Q1DQ2*C1+C2
      S(2,3)=Q1DQ2*S1/WN1D(L)-S2/WN2D(L)
      S(2,4)=Q1D(L)*(S1/WN1D(L)-S2/WN2D(L))
      S(3,1)=-Q1D(L)*(WN1D(L)*S1-WN2D(L)*S2)
      S(3,2)=Q1DQ2*WN1D(L)*S1-WN2D(L)*S2
      S(3,3)=-Q1DQ2*C1+C2
      S(3,4)=-Q1D(L)*(C1-C2)
      S(4,1)=WN1D(L)*S1-Q1DQ2*WN2D(L)*S2
      S(4,2)=-Q2ID(L)*(WN1D(L)*S1-WN2D(L)*S2)
      S(4,3)=Q2ID(L)*(C1-C2)
      S(4,4)=C1-Q1DQ2*C2
      DO 7 I=1,4
      DO 8 J=1,4
      S(I,J)=S(I,J)/DQ
    8 CONTINUE
    7 CONTINUE
C
C     IMPEDANCE TENSOR AT THE TOP OF THE L-TH LAYER
C     ---------------------------------------------
      A11=S(3,1)*ZXXD+S(3,2)*ZYXD+S(3,3)
      A12=S(3,1)*ZXYD+S(3,2)*ZYYD+S(3,4)
      A21=S(4,1)*ZXXD+S(4,2)*ZYXD+S(4,3)
      A22=S(4,1)*ZXYD+S(4,2)*ZYYD+S(4,4)
      DETA=A11*A22-A12*A21
      B11=A22/DETA
      B12=-A12/DETA
      B21=-A21/DETA
      B22=A11/DETA
      A11=S(1,1)*ZXXD+S(1,2)*ZYXD+S(1,3)
      A12=S(1,1)*ZXYD+S(1,2)*ZYYD+S(1,4)
      A21=S(2,1)*ZXXD+S(2,2)*ZYXD+S(2,3)
      A22=S(2,1)*ZXYD+S(2,2)*ZYYD+S(2,4)
      ZXXD=A11*B11+A12*B21
      ZXYD=A11*B12+A12*B22
      ZYXD=A21*B11+A22*B21
      ZYYD=A21*B12+A22*B22
      zxx=zxxd
      zxy=zxyd
      zyx=zyxd
      zyy=zyyd
C
      DO 9 I=1,4
      DO 10 J=1,4
      SOUC(L,I,J)=0.
      DO 11 K=1,4
      SOUC(L,I,J)=SOUC(L,I,J)+S(I,K)*SOUC(L+1,K,J)
   11 CONTINUE
   10 CONTINUE
    9 CONTINUE
    4 CONTINUE
C
    3 CONTINUE
C
C     NORMALIZATION CONSTANTS IN THE BASEMENT
C     ---------------------------------------
      T320=SOUC(1,3,1)+Q1D(NL)*SOUC(1,3,2)-
     /     WN1D(NL)*Q1D(NL)*SOUC(1,3,3)+WN1D(NL)*SOUC(1,3,4)
      T340=Q2ID(NL)*SOUC(1,3,1)+SOUC(1,3,2)-
     /     WN2D(NL)*SOUC(1,3,3)+WN2D(NL)*Q2ID(NL)*SOUC(1,3,4)
      T420=SOUC(1,4,1)+Q1D(NL)*SOUC(1,4,2)-
     /     WN1D(NL)*Q1D(NL)*SOUC(1,4,3)+WN1D(NL)*SOUC(1,4,4)
      T440=Q2ID(NL)*SOUC(1,4,1)+SOUC(1,4,2)-
     /     WN2D(NL)*SOUC(1,4,3)+WN2D(NL)*Q2ID(NL)*SOUC(1,4,4)
      GAMD=T320*T440-T340*T420
      GAM1=(T440*HX0-T340*HY0)/GAMD
      GAM2=(-T420*HX0+T320*HY0)/GAMD
C
      PB1=GAM1+Q2ID(NL)*GAM2
      PB2=Q1D(NL)*GAM1+GAM2
      PB3=-WN1D(NL)*Q1D(NL)*GAM1-WN2D(NL)*GAM2
      PB4=WN1D(NL)*GAM1+WN2D(NL)*Q2ID(NL)*GAM2
C
C     ELECTROMAGNETIC FIELD AT THE LEVELS ZS(IZ),IZ=1,...,NZS
C     -------------------------------------------------------
      ZB=0.
      IF(NL.EQ.1)GOTO 20
      DO 21 I=1,NL1
      ZB=ZB+H(I)
   21 CONTINUE
   20 IZ=NZS
      IL=NL
   23 IF(ZS(IZ).LT.ZB)GOTO 22
C
C     FIELD WITHIN THE BASEMENT
C     -------------------------
      WN1=-IC*OMMI*WN1D(IL)
      WN2=-IC*OMMI*WN2D(IL)
      F1=CDEXP(-1.E+3*WN1*(ZS(IZ)-ZB))
      F2=CDEXP(-1.E+3*WN2*(ZS(IZ)-ZB))
      EXD(IZ)=GAM1*F1+Q2ID(IL)*GAM2*F2
      EYD(IZ)=Q1D(IL)*GAM1*F1+GAM2*F2
      EZD(IZ)=-(SG(IL,3,1)*EXD(IZ)+SG(IL,3,2)*EYD(IZ))/SG(IL,3,3)
      HXD(IZ)=-WN1D(IL)*Q1D(IL)*GAM1*F1-WN2D(IL)*GAM2*F2
      HYD(IZ)=WN1D(IL)*GAM1*F1+WN2D(IL)*Q2ID(IL)*GAM2*F2
      ex(iz)=exd(iz)
      ey(iz)=eyd(iz)
      ez(iz)=ezd(iz)
      hx(iz)=hxd(iz)
      hy(iz)=hyd(iz)
C
      IF(IZ.EQ.1)GOTO 30
      IZ=IZ-1
      GOTO 23
   22 IL=IL-1
      IF(IL.EQ.0)GOTO 24
      ZB=ZB-H(IL)
      IF(ZS(IZ).LT.ZB)GOTO 22
C
C     FIELD WITHIN THE IL-TH LAYER
C     ----------------------------
      IL1=IL+1
      PL1=SOUC(IL1,1,1)*PB1+SOUC(IL1,1,2)*PB2+
     /    SOUC(IL1,1,3)*PB3+SOUC(IL1,1,4)*PB4
      PL2=SOUC(IL1,2,1)*PB1+SOUC(IL1,2,2)*PB2+
     /    SOUC(IL1,2,3)*PB3+SOUC(IL1,2,4)*PB4
      PL3=SOUC(IL1,3,1)*PB1+SOUC(IL1,3,2)*PB2+
     /    SOUC(IL1,3,3)*PB3+SOUC(IL1,3,4)*PB4
      PL4=SOUC(IL1,4,1)*PB1+SOUC(IL1,4,2)*PB2+
     /    SOUC(IL1,4,3)*PB3+SOUC(IL1,4,4)*PB4
      WN1=-IC*OMMI*WN1D(IL)
      WN2=-IC*OMMI*WN2D(IL)
      Q1=IC*Q1D(IL)/OMMI
      Q2I=-IC*Q2ID(IL)*OMMI
      Q1DQ2=Q1D(IL)*Q2ID(IL)
      DQ=1.-Q1DQ2
C
   27 DZ=ZB+H(IL)-ZS(IZ)
      WT1=1.E+3*WN1*DZ
      WT2=1.E+3*WN2*DZ
      S1=SH(WT1)
      S2=SH(WT2)
      C1=CH(WT1)
      C2=CH(WT2)
      S(1,1)=C1-Q1DQ2*C2
      S(1,2)=-Q2ID(IL)*(C1-C2)
      S(1,3)=Q2ID(IL)*(S1/WN1D(IL)-S2/WN2D(IL))
      S(1,4)=S1/WN1D(IL)-Q1DQ2*S2/WN2D(IL)
      S(2,1)=Q1D(IL)*(C1-C2)
      S(2,2)=-Q1DQ2*C1+C2
      S(2,3)=Q1DQ2*S1/WN1D(IL)-S2/WN2D(IL)
      S(2,4)=Q1D(IL)*(S1/WN1D(IL)-S2/WN2D(IL))
      S(3,1)=-Q1D(IL)*(WN1D(IL)*S1-WN2D(IL)*S2)
      S(3,2)=Q1DQ2*WN1D(IL)*S1-WN2D(IL)*S2
      S(3,3)=-Q1DQ2*C1+C2
      S(3,4)=-Q1D(IL)*(C1-C2)
      S(4,1)=WN1D(IL)*S1-Q1DQ2*WN2D(IL)*S2
      S(4,2)=-Q2ID(IL)*(WN1D(IL)*S1-WN2D(IL)*S2)
      S(4,3)=Q2ID(IL)*(C1-C2)
      S(4,4)=C1-Q1DQ2*C2
      DO 25 I=1,4
      DO 26 J=1,4
      S(I,J)=S(I,J)/DQ
   26 CONTINUE
   25 CONTINUE
      EXD(IZ)=S(1,1)*PL1+S(1,2)*PL2+S(1,3)*PL3+S(1,4)*PL4
      EYD(IZ)=S(2,1)*PL1+S(2,2)*PL2+S(2,3)*PL3+S(2,4)*PL4
      EZD(IZ)=-(SG(IL,3,1)*EXD(IZ)+SG(IL,3,2)*EYD(IZ))/SG(IL,3,3)
      HXD(IZ)=S(3,1)*PL1+S(3,2)*PL2+S(3,3)*PL3+S(3,4)*PL4
      HYD(IZ)=S(4,1)*PL1+S(4,2)*PL2+S(4,3)*PL3+S(4,4)*PL4
      ex(iz)=exd(iz)
      ey(iz)=eyd(iz)
      ez(iz)=ezd(iz)
      hx(iz)=hxd(iz)
      hy(iz)=hyd(iz)
C
      IF(IZ.EQ.1)GOTO 30
      IZ=IZ-1
      IF(ZS(IZ).LT.ZB)GOTO 22
      GOTO 27
C
C     FIELD WITHIN THE AIR LAYER
C     --------------------------
   24 IL1=IL+1
      PL1=SOUC(IL1,1,1)*PB1+SOUC(IL1,1,2)*PB2+
     /    SOUC(IL1,1,3)*PB3+SOUC(IL1,1,4)*PB4
      PL2=SOUC(IL1,2,1)*PB1+SOUC(IL1,2,2)*PB2+
     /    SOUC(IL1,2,3)*PB3+SOUC(IL1,2,4)*PB4
      PL3=SOUC(IL1,3,1)*PB1+SOUC(IL1,3,2)*PB2+
     /    SOUC(IL1,3,3)*PB3+SOUC(IL1,3,4)*PB4
      PL4=SOUC(IL1,4,1)*PB1+SOUC(IL1,4,2)*PB2+
     /    SOUC(IL1,4,3)*PB3+SOUC(IL1,4,4)*PB4
C
   28 S(1,1)=1.
      S(1,2)=0.
      S(1,3)=0.
      S(1,4)=-1.E+3*IC*OMMI*(ZB-ZS(IZ))
      S(2,1)=0.
      S(2,2)=1.
      S(2,3)=-S(1,4)
      S(2,4)=0.
      S(3,1)=0.
      S(3,2)=0.
      S(3,3)=1.
      S(3,4)=0.
      S(4,1)=0.
      S(4,2)=0.
      S(4,3)=0.
      S(4,4)=1.
      EXD(IZ)=S(1,1)*PL1+S(1,2)*PL2+S(1,3)*PL3+S(1,4)*PL4
      EYD(IZ)=S(2,1)*PL1+S(2,2)*PL2+S(2,3)*PL3+S(2,4)*PL4
      EZD(IZ)=0.
      HXD(IZ)=S(3,1)*PL1+S(3,2)*PL2+S(3,3)*PL3+S(3,4)*PL4
      HYD(IZ)=S(4,1)*PL1+S(4,2)*PL2+S(4,3)*PL3+S(4,4)*PL4
      ex(iz)=exd(iz)
      ey(iz)=eyd(iz)
      ez(iz)=ezd(iz)
      hx(iz)=hxd(iz)
      hy(iz)=hyd(iz)
C
      IF(IZ.EQ.1)GOTO 30
      IZ=IZ-1
      GOTO 28
C
C     SUBROUTINE FINISHED WITHOUT ERROR
C     ---------------------------------
   30 IE=0
      RETURN
C
C     NON-SYMMETRIC CONDUCTIVITY TENSOR - ERROR MESSAGE
C     -------------------------------------------------
  100 IE=1
      RETURN
C
      END
      SUBROUTINE KOEF3
C     ==================================================================
C
C-----------------------------------------------------------------------
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20,NPMAX=30,NBMAX=10)
      PARAMETER(NEHMAX=2*(NMAX-2)*(MMAX-2),MEHMAX=2*(MMAX-2)+3)
      COMPLEX CE
      COMPLEX*16 APOM,BPOM
      COMPLEX BOULE,BOURE,BOUUE,BOUDE,BOULH,BOURH,BOUUH,BOUDH
      COMPLEX AHES
C
      COMMON /CAP/ OMMI,PERIOD
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /NHS/ NEH
      COMMON /RHS/ BPOM(NEHMAX)
      COMMON /BUE/ BOULE(MMAX),BOURE(MMAX),BOUUE(NMAX),BOUDE(NMAX)
      COMMON /BUH/ BOULH(MMAX),BOURH(MMAX),BOUUH(NMAX),BOUDH(NMAX)
      COMMON /MED/ IRG(NMAX-1,MMAX-1),SG(NCMAX,3,3)
      COMMON /SCO/ AHES(NMAX-2)
      COMMON /MAT/ APOM(NEHMAX,MEHMAX)
C-----------------------------------------------------------------------
C
      MD1=MD-1
      MD2=MD-2
      MH2=M1-MD
      MH1=MH2+1
      MH=MH1+1
      MEH1=M2+MH2+1
      MEH3=MEH1+2
      OMMII=1./OMMI
C
      JK=0
      DO 1 J=1,N2
      HYM=1.E+3*SY(J)
      HYP=1.E+3*SY(J+1)
C
      DO 2 K=1,MD2
      HZM=1.E+3*SZ(K)
      HZP=1.E+3*SZ(K+1)
      JK=JK+1
C
      DO 40 L=1,MEH3
      APOM(JK,L)=0.
   40 CONTINUE
C
      PE=0.5*(HZM+HZP)/HYM
      QE=0.5*(HYM+HYP)/HZM
      RE=0.5*(HYM+HYP)/HZP
      SE=0.5*(HZM+HZP)/HYP
      TE=-(PE+QE+RE+SE)
      CE=TE
C
      APOM(JK,1)=CE
      APOM(JK,2)=RE
      APOM(JK,MEH1)=SE  
C
      BPOM(JK)=0.
C
      IF(J.GT.1)GOTO 11
      BPOM(JK)=BPOM(JK)-PE*BOULE(K+1)
C
   11 IF(J.LT.N2)GOTO 12
      APOM(JK,MEH1)=0.
      BPOM(JK)=BPOM(JK)-SE*BOURE(K+1)
C
   12 IF(K.GT.1)GOTO 34
      BPOM(JK)=BPOM(JK)-QE*BOUUE(J)
C
   34 CONTINUE
      DO 38 L=1,MEH3
      APOM(JK,L)=-IC*OMMII*APOM(JK,L)
   38 CONTINUE
      BPOM(JK)=-IC*OMMII*BPOM(JK)
C
    2 CONTINUE
C
      DO 3 K=MD1,M2
      HZM=1.E+3*SZ(K)
      HZP=1.E+3*SZ(K+1)
      KN2=K-MD1
      JK=JK+1
C
      DO 41 L=1,MEH3
      APOM(JK,L)=0.
   41 CONTINUE
C
      IF(K.EQ.MD1)GOTO 21
C
      IG=IRG(J,K)
      IGMM=IG
      DD=SG(IG,2,2)*SG(IG,3,3)-SG(IG,2,3)*SG(IG,3,2)
      SBYMM=SG(IG,2,3)/DD
      SBZMM=SG(IG,2,2)/DD
      SCYMM=SG(IG,3,3)/DD
      SCZMM=SG(IG,3,2)/DD
      SAYMM=SG(IG,3,1)*SBYMM-SG(IG,2,1)*SCYMM
      SAZMM=SG(IG,3,1)*SBZMM-SG(IG,2,1)*SCZMM
C
      IG=IRG(J+1,K)
      IGPM=IG
      DD=SG(IG,2,2)*SG(IG,3,3)-SG(IG,2,3)*SG(IG,3,2)
      SBYPM=SG(IG,2,3)/DD
      SBZPM=SG(IG,2,2)/DD
      SCYPM=SG(IG,3,3)/DD
      SCZPM=SG(IG,3,2)/DD
      SAYPM=SG(IG,3,1)*SBYPM-SG(IG,2,1)*SCYPM
      SAZPM=SG(IG,3,1)*SBZPM-SG(IG,2,1)*SCZPM
C
      GOTO 22
C
   21 IG=IRG(J,K)
      IGMM=IG
      SBYMM=0.
      SBZMM=0.
      SCYMM=0.
      SCZMM=0.
      SAYMM=0.
      SAZMM=0.
C
      IG=IRG(J+1,K)
      IGPM=IG
      SBYPM=0.
      SBZPM=0.
      SCYPM=0.
      SCZPM=0.
      SAYPM=0.
      SAZPM=0.
C
   22 IG=IRG(J,K+1)
      IGMP=IG
      DD=SG(IG,2,2)*SG(IG,3,3)-SG(IG,2,3)*SG(IG,3,2)
      SBYMP=SG(IG,2,3)/DD
      SBZMP=SG(IG,2,2)/DD
      SCYMP=SG(IG,3,3)/DD
      SCZMP=SG(IG,3,2)/DD
      SAYMP=SG(IG,3,1)*SBYMP-SG(IG,2,1)*SCYMP
      SAZMP=SG(IG,3,1)*SBZMP-SG(IG,2,1)*SCZMP
C
      IG=IRG(J+1,K+1)
      IGPP=IG
      DD=SG(IG,2,2)*SG(IG,3,3)-SG(IG,2,3)*SG(IG,3,2)
      SBYPP=SG(IG,2,3)/DD
      SBZPP=SG(IG,2,2)/DD
      SCYPP=SG(IG,3,3)/DD
      SCZPP=SG(IG,3,2)/DD
      SAYPP=SG(IG,3,1)*SBYPP-SG(IG,2,1)*SCYPP
      SAZPP=SG(IG,3,1)*SBZPP-SG(IG,2,1)*SCZPP
C
      PE=0.5*(HZM+HZP)/HYM
      QE=0.5*(HYM+HYP)/HZM
      RE=0.5*(HYM+HYP)/HZP
      SE=0.5*(HZM+HZP)/HYP
      TE=-(PE+QE+RE+SE)
      VE=0.25*OMMI*
     /   ((SG(IGMM,1,1)+SG(IGMM,1,2)*SAYMM-SG(IGMM,1,3)*SAZMM)*HYM*HZM
     /   +(SG(IGPM,1,1)+SG(IGPM,1,2)*SAYPM-SG(IGPM,1,3)*SAZPM)*HYP*HZM
     /   +(SG(IGMP,1,1)+SG(IGMP,1,2)*SAYMP-SG(IGMP,1,3)*SAZMP)*HYM*HZP
     /   +(SG(IGPP,1,1)+SG(IGPP,1,2)*SAYPP-SG(IGPP,1,3)*SAZPP)*HYP*HZP)
      CE=CMPLX(TE,VE)
C
      IF(K.EQ.MD1)GOTO 23
C
      SWYE=-0.25*OMMI*((SG(IGMM,1,2)*SBYMM-SG(IGMM,1,3)*SBZMM)*HZM
     /                +(SG(IGMP,1,2)*SBYMP-SG(IGMP,1,3)*SBZMP)*HZP)
      SXYE=0.25*OMMI*((SG(IGPM,1,2)*SBYPM-SG(IGPM,1,3)*SBZPM)*HZM
     /               +(SG(IGPP,1,2)*SBYPP-SG(IGPP,1,3)*SBZPP)*HZP)
      SWZE=-0.25*OMMI*((SG(IGMM,1,2)*SCYMM-SG(IGMM,1,3)*SCZMM)*HYM
     /                +(SG(IGPM,1,2)*SCYPM-SG(IGPM,1,3)*SCZPM)*HYP)
      GOTO 24
C
   23 SWYE=-0.25*OMMI*(SG(IGMP,1,2)*SBYMP-SG(IGMP,1,3)*SBZMP)*HZP
      SXYE=0.25*OMMI*(SG(IGPP,1,2)*SBYPP-SG(IGPP,1,3)*SBZPP)*HZP
      SWZE=0.
C
   24 SXZE=0.25*OMMI*((SG(IGMP,1,2)*SCYMP-SG(IGMP,1,3)*SCZMP)*HYM
     /               +(SG(IGPP,1,2)*SCYPP-SG(IGPP,1,3)*SCZPP)*HYP)
      SZE=-(SWYE+SWZE+SXYE+SXZE)
C
      APOM(JK,1)=CE
      IF(K.EQ.MD1)THEN
        APOM(JK,2)=RE
      ELSE IF(K.GT.MD1)THEN
        APOM(JK,3)=RE
      ENDIF
      APOM(JK,MEH1)=SE
C
      IF(K.EQ.MD1)THEN
        APOM(JK,3)=CMPLX(0.,SXZE)
      ELSE IF(K.GT.MD1)THEN
        APOM(JK,2)=CMPLX(0.,SZE)
        APOM(JK,4)=CMPLX(0.,SXZE)
        APOM(JK,MEH1+1)=CMPLX(0.,SXYE)
      ENDIF
C
      BPOM(JK)=0. 
C
      IF(K.GT.MD1)GOTO 51
      BPOM(JK)=BPOM(JK)-IC*SZE*BOUUH(J+1)
C
   51 IF(K.GT.MD)GOTO 52
      BPOM(JK)=BPOM(JK)-IC*SWZE*BOUUH(J+1)
C
   52 IF(J.GT.1)GOTO 25
      BPOM(JK)=BPOM(JK)-PE*BOULE(K+1)
CORRECTED ON November 6, 1999 by J.Pek
C     BPOM(JK)=BPOM(JK)-SWYE*BOULH(KN2+1)
	BPOM(JK)=BPOM(JK)-IC*SWYE*BOULH(KN2+1)
C
   25 IF(J.LT.N2)GOTO 26
      APOM(JK,MEH1)=0.
      BPOM(JK)=BPOM(JK)-SE*BOURE(K+1)
      IF(K.GT.MD1)THEN
        APOM(JK,MEH1+1)=0.
      ENDIF
      BPOM(JK)=BPOM(JK)-IC*SXYE*BOURH(KN2+1)
C
   26 IF(K.LT.M2)GOTO 36
      APOM(JK,3)=0.
      BPOM(JK)=BPOM(JK)-RE*BOUDE(J+1)
      APOM(JK,4)=0.
      BPOM(JK)=BPOM(JK)-IC*SXZE*BOUDH(J+1)
C
   36 CONTINUE
      DO 37 L=1,MEH3
      APOM(JK,L)=-IC*OMMII*APOM(JK,L)
   37 CONTINUE
      BPOM(JK)=-IC*OMMII*BPOM(JK)
C
   27 CONTINUE
C
      IF(K.EQ.MD1)GOTO 3
C
      JK=JK+1
C
      DO 42 L=1,MEH3
      APOM(JK,L)=0.
   42 CONTINUE
C
      PMH=0.25*(SCZMM+SBYMM)
      PCH=0.5*(SBZMM*HZM+SBZMP*HZP)/HYM
      PH=PCH
      PPH=-0.25*(SCZMP+SBYMP)
      QCH=0.5*(SCYMM*HYM+SCYPM*HYP)/HZM
      QH=QCH
      RCH=0.5*(SCYMP*HYM+SCYPP*HYP)/HZP
      RH=RCH
      SMH=-0.25*(SCZPM+SBYPM)
      SCH=0.5*(SBZPM*HZM+SBZPP*HZP)/HYP
      SH=SCH
      SPH=0.25*(SCZPP+SBYPP)
      TCH=-(PCH+QCH+RCH+SCH)
      TH=TCH+0.25*(SCZMP+SCZPM-SCZMM-SCZPP+SBYMP+SBYPM-SBYMM-SBYPP)
      VH=0.25*OMMI*(HYM+HYP)*(HZM+HZP)
C
      APOM(JK,1)=CMPLX(TH,VH)
      APOM(JK,3)=RH
      APOM(JK,MEH1-2)=SMH
      APOM(JK,MEH1)=SH
      APOM(JK,MEH1+2)=SPH  
C
      APOM(JK,2)=-OMMII*SXZE
      APOM(JK,MEH1-1)=-OMMII*SXYE 
      IF(K.EQ.MD)AHES(J)=-OMMII*SWZE
C
      BPOM(JK)=0.
C
      IF(J.GT.1)GOTO 28
      BPOM(JK)=BPOM(JK)-PMH*BOULH(KN2)-PH*BOULH(KN2+1)-PPH*BOULH(KN2+2)
      BPOM(JK)=BPOM(JK)+OMMII*SWYE*BOULE(K+1)
C
   28 IF(J.LT.N2)GOTO 29
      APOM(JK,MEH1-2)=0.
      APOM(JK,MEH1)=0.
      APOM(JK,MEH1+2)=0.
      APOM(JK,MEH1-1)=0.
      BPOM(JK)=BPOM(JK)-SMH*BOURH(KN2)-SH*BOURH(KN2+1)-SPH*BOURH(KN2+2)
      BPOM(JK)=BPOM(JK)+OMMII*SXYE*BOURE(K+1)
C
   29 IF(K.GT.MD)GOTO 30
      BPOM(JK)=BPOM(JK)-QH*BOUUH(J+1)
      IF(J.GT.1)BPOM(JK)=BPOM(JK)-PMH*BOUUH(J)
      IF(J.LT.N2)BPOM(JK)=BPOM(JK)-SMH*BOUUH(J+2)
C
   30 IF(K.LT.M2)GOTO 35
      APOM(JK,2)=0.
      APOM(JK,3)=0.
      APOM(JK,MEH1+2)=0.
      BPOM(JK)=BPOM(JK)-RH*BOUDH(J+1)
      IF(J.GT.1)BPOM(JK)=BPOM(JK)-PPH*BOUDH(J)
      IF(J.LT.N2)BPOM(JK)=BPOM(JK)-SPH*BOUDH(J+2)
      BPOM(JK)=BPOM(JK)+OMMII*SXZE*BOUDE(J+1)
C
   35 CONTINUE
C
    3 CONTINUE
C
    1 CONTINUE
C
      RETURN
      END
      SUBROUTINE SURFLD(IHPOL)
C     ==================================================================
C
C-----------------------------------------------------------------------
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20)
      PARAMETER(NEHMAX=2*(NMAX-2)*(MMAX-2),MEHMAX=2*(MMAX-2)+3)
      COMPLEX SEX,SHX,SEY,SHY,SHZ,SEZ
      COMPLEX BEXL,BEXR,BEYL,BEYR,BHXL,BHXR,BHYL,BHYR
      COMPLEX BOULE,BOURE,BOUUE,BOUDE,BOULH,BOURH,BOUUH,BOUDH
      COMPLEX PEM,PEC,PEP,PEL,PER,PHM,PHC,PHP
      COMPLEX DEDZ,DEDY,DHDZ
      COMPLEX ENOR,EXN,EYN
      COMPLEX*16 BPOM
C
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /NHS/ NEH
      COMMON /RHS/ BPOM(NEHMAX)
      COMMON /CAP/ OMMI,PERIOD
      COMMON /BUE/ BOULE(MMAX),BOURE(MMAX),BOUUE(NMAX),BOUDE(NMAX)
      COMMON /BUH/ BOULH(MMAX),BOURH(MMAX),BOUUH(NMAX),BOUDH(NMAX)
      COMMON /MED/ IRG(NMAX-1,MMAX-1),SG(NCMAX,3,3)
      COMMON /SFD/ SEX(2,NMAX),SEY(2,NMAX,2),SEZ(2,NMAX,2),
     &             SHX(2,NMAX),SHY(2,NMAX),SHZ(2,NMAX)
      COMMON /BOS/ BEXL(2),BEYL(2),BHXL(2),BHYL(2),
     /             BEXR(2),BEYR(2),BHXR(2),BHYR(2)
      COMMON /DIS/ SYY(NMAX),SZZ(MMAX)
C-----------------------------------------------------------------------
C
      OMMII=1./OMMI
      MH2=M1-MD
      MEH2=2*M2-MD+1
C
      SEX(IHPOL,1)=BEXL(IHPOL)
      SEX(IHPOL,N)=BEXR(IHPOL)
      SHX(IHPOL,1)=BHXL(IHPOL)
      SHX(IHPOL,N)=BHXR(IHPOL)
      SEY(IHPOL,1,1)=BEYL(IHPOL)
      SEY(IHPOL,1,2)=SEY(IHPOL,1,1)
      SEY(IHPOL,N,1)=BEYR(IHPOL)
      SEY(IHPOL,N,2)=SEY(IHPOL,N,1)
      SHY(IHPOL,1)=BHYL(IHPOL)
      SHY(IHPOL,N)=BHYR(IHPOL)
      SHZ(IHPOL,1)=0.
      SHZ(IHPOL,N)=0.
C
      IEH=MD-1
      DO 1 J=1,N2
      J1=J+1
      SEX(IHPOL,J1)=BPOM(IEH)
c	write(11,*)ihpol,j1,sex(ihpol,j1)
      SHX(IHPOL,J1)=BOUUH(J1)
      PEM=BPOM(IEH-1)
      PEC=BPOM(IEH)
      PEP=BPOM(IEH+1)
      PHM=BOUUH(J1)
      PHC=BPOM(IEH+2)
      PHP=BPOM(IEH+4)
      IF(J.GT.1)GOTO 2
      PEL=BOULE(MD)
      GOTO 3
    2 PEL=BPOM(IEH-MEH2)
    3 IF(J.LT.N2)GOTO 4
      PER=BOURE(MD)
      GOTO 5
    4 PER=BPOM(IEH+MEH2)
    5 CONTINUE
      IM=1
      IS=2
      HYM=1.E+3*SY(J)
      HYP=1.E+3*SY(J1)
      HZM=1.E+3*SZ(MD-1)
      HZC=1.E+3*SZ(MD)
      HZP=1.E+3*SZ(MD+1)
C
      CALL DER3PO(PEM,PEC,PEP,HZM,HZC,IS,DEDZ)
      CALL DER3PO(PEL,PEC,PER,HYM,HYP,IS,DEDY)
      CALL DER3PO(PHM,PHC,PHP,HZC,HZP,IM,DHDZ)
C
      SHY(IHPOL,J1)=-IC*OMMII*DEDZ
      SHZ(IHPOL,J1)=IC*OMMII*DEDY
      S1=SG(IRG(J,MD),2,2)
      S2=SG(IRG(J,MD),2,1)
c      SEY(IHPOL,J1,1)=(DHDZ-S2*SEX(IHPOL,J1))/S1
      S1D=SG(IRG(J,MD),2,2)*SG(IRG(J,MD),3,3)-
     -    SG(IRG(J,MD),3,2)*SG(IRG(J,MD),2,3)
      S2K=SG(IRG(J,MD),2,3)*SG(IRG(J,MD),3,1)-
     -    SG(IRG(J,MD),3,3)*SG(IRG(J,MD),2,1)
      S2L=SG(IRG(J,MD),3,2)*SG(IRG(J,MD),2,1)-
     -    SG(IRG(J,MD),2,2)*SG(IRG(J,MD),3,1)
      S3=SG(IRG(J,MD),3,3)
      SEY(IHPOL,J1,1)=(S3*DHDZ+S2K*SEX(IHPOL,J1))/S1D
      SEZ(IHPOL,J1,1)=(-SG(IRG(J,MD),3,2)*DHDZ+S2L*SEX(IHPOL,J1))/S1D
      S1=SG(IRG(J1,MD),2,2)
      S2=SG(IRG(J1,MD),2,1)
c      SEY(IHPOL,J1,2)=(DHDZ-S2*SEX(IHPOL,J1))/S1
      S1D=SG(IRG(J1,MD),2,2)*SG(IRG(J1,MD),3,3)-
     -    SG(IRG(J1,MD),3,2)*SG(IRG(J1,MD),2,3)
      S2K=SG(IRG(J1,MD),2,3)*SG(IRG(J1,MD),3,1)-
     -    SG(IRG(J1,MD),3,3)*SG(IRG(J1,MD),2,1)
      S2L=SG(IRG(J1,MD),3,2)*SG(IRG(J1,MD),2,1)-
     -    SG(IRG(J1,MD),2,2)*SG(IRG(J1,MD),3,1)
      S3=SG(IRG(J1,MD),3,3)
      SEY(IHPOL,J1,2)=(S3*DHDZ+S2K*SEX(IHPOL,J1))/S1D
      SEZ(IHPOL,J1,2)=(-SG(IRG(J1,MD),3,2)*DHDZ+S2L*SEX(IHPOL,J1))/S1D
C
      IF(J.GE.N2)GOTO 1
      IEH=IEH+MEH2
    1 CONTINUE
C
      WRITE(11,2000)IHPOL
      WRITE(11,2005)
C
      IF(IHPOL.EQ.2)GOTO 6
      ENOR=SEX(1,1)
      GOTO 7
    6 ENOR=SEY(2,1,1)
    7 CONTINUE
C
      DO 9 J=1,N
C
      IF(J.EQ.1.OR.J.EQ.N)GOTO 8
      IF(IRG(J-1,MD).EQ.IRG(J,MD))GOTO 8
C
      WRITE(11,2001)J,SYY(J),SEX(IHPOL,J),SEY(IHPOL,J,1),
     /             SHX(IHPOL,J),SHY(IHPOL,J),SHZ(IHPOL,J)
      EXN=SEX(IHPOL,J)/ENOR
      EYN=SEY(IHPOL,J,1)/ENOR
      WRITE(11,2002)EXN,EYN,SEZ(IHPOL,J,1)
C
      WRITE(11,2003)J,SYY(J),SEX(IHPOL,J),SEY(IHPOL,J,2),
     /             SHX(IHPOL,J),SHY(IHPOL,J),SHZ(IHPOL,J)
      EYN=SEY(IHPOL,J,2)/ENOR
      WRITE(11,2002)EXN,EYN,SEZ(IHPOL,J,2)
C
      GOTO 9
C
    8 WRITE(11,2004)J,SYY(J),SEX(IHPOL,J),SEY(IHPOL,J,1),
     /             SHX(IHPOL,J),SHY(IHPOL,J),SHZ(IHPOL,J)
      EXN=SEX(IHPOL,J)/ENOR
      EYN=SEY(IHPOL,J,1)/ENOR
      WRITE(11,2002)EXN,EYN,SEZ(IHPOL,J,1)
C
    9 CONTINUE
      WRITE(11,2005)
C
      RETURN
C
 2000 FORMAT(1H1,1X,'SURFACE VALUES OF THE ELECTROMAGNETIC COMPONENTS FO
     /R THE POLARISATION ',I2/2X,'======================================
     /==================================================================
     /========================'/2X,'POINT',2X,'DISTANCE',6X,'RE(EX)',5X,
     /'IM(EX)',5X,'RE(EY)',5X,'IM(EY)',5X,'RE(HX)',5X,'IM(HX)',5X,'RE(HY
     /)',5X,'IM(HY)',5X,'RE(HZ)',5X,'IM(HZ)'/21X,'RE(EX/EN)',2X,'IM(EX/E
     /N)',2X,'RE(EY/EN)',2X,'IM(EY/EN)',49X,'RE(EZ)',5X,'IM(EZ)')
 2001 FORMAT(2X,I3,'-',2X,F10.2,2X,10E11.3)
 2002 FORMAT(20X,4E11.3,44X,2E11.3)
 2003 FORMAT(2X,I3,'+',2X,F10.2,2X,10E11.3)
 2004 FORMAT(2X,I3,3X,F10.2,2X,10E11.3)
 2005 FORMAT(2X,'-------------------------------------------------------
     /------------------------------------------------------------------
     /-------')
C
      END
      SUBROUTINE SURFCE
C     ==================================================================
C
C-----------------------------------------------------------------------
      COMPLEX IC
      PARAMETER(IC=(0.,1.),PI=3.141592653589793)
      PARAMETER(NMAX=151,MMAX=59,NCMAX=20)
      COMPLEX SEX,SEY,SHX,SHY,SHZ,SEZ
      COMPLEX DET,D,ZXX,ZXY,ZYXL,ZYXR,ZYYL,ZYYR
      COMPLEX WX,WY,TR,TI
      complex zpom
C
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
      COMMON /CAP/ OMMI,PERIOD
      COMMON /MED/ IRG(NMAX-1,MMAX-1),SG(NCMAX,3,3)
      COMMON /SFD/ SEX(2,NMAX),SEY(2,NMAX,2),SEZ(2,NMAX,2),
     &             SHX(2,NMAX),SHY(2,NMAX),SHZ(2,NMAX)
      COMMON /DIS/ SYY(NMAX),SZZ(MMAX)
C-----------------------------------------------------------------------
C
      OMMII=1./OMMI
C
      WRITE(11,2004)
      WRITE(11,2005)
      WRITE(11,2006)
      WRITE(11,2005)
      WRITE(11,2007)
C
      DO 1 J=1,N
c      zpom=sey(1,j,1)/sey(1,1,1)
c      z1=real(zpom)
c      z2=aimag(zpom)
c      z3=cabs(zpom)
c      z4=phase(zpom)
c      zpom=sex(1,j)/sex(1,1)
c      z3=cabs(zpom)
c      z4=phase(zpom)
c      write(7,7000)syy(j),z1,z2,z3,z4
c 7000 format(1x,5f11.4)
c      zpom=shz(1,j)/shy(1,j)
c      z1=cabs(zpom)
c      z2=phase(zpom)
c      zpom=shz(1,j)/shy(1,1)
c      z3=cabs(zpom)
c      z4=phase(zpom)
c      write(7,7000)syy(j),z1,z2,z3,z4
C
      DET=SHX(1,J)*SHY(2,J)-SHX(2,J)*SHY(1,J)
C
      D=SEX(1,J)*SHY(2,J)-SEX(2,J)*SHY(1,J)
      ZXX=D/DET
C
      D=SEX(2,J)*SHX(1,J)-SEX(1,J)*SHX(2,J)
      ZXY=D/DET
      RAXY=OMMII*CABS(ZXY)*CABS(ZXY)
      RAXYL=ALOG10(RAXY)
      PAXY=PHASE(ZXY)
C
      D=SEY(1,J,1)*SHY(2,J)-SEY(2,J,1)*SHY(1,J)
      ZYXL=D/DET
      RAYXL=OMMII*CABS(ZYXL)*CABS(ZYXL)
      RAYXLL=ALOG10(RAYXL)
      PAYXL=PHASE(ZYXL)
      D=SEY(1,J,2)*SHY(2,J)-SEY(2,J,2)*SHY(1,J)
      ZYXR=D/DET
      RAYXR=OMMII*CABS(ZYXR)*CABS(ZYXR)
      RAYXRL=ALOG10(RAYXR)
      PAYXR=PHASE(ZYXR)
C
      D=SEY(2,J,1)*SHX(1,J)-SEY(1,J,1)*SHX(2,J)
      ZYYL=D/DET
      D=SEY(2,J,2)*SHX(1,J)-SEY(1,J,2)*SHX(2,J)
      ZYYR=D/DET
C
      CALL MTPARM(ZXX,ZXY,ZYXL,ZYYL,PSWIFL,PANIZL,PSKEWL)
      CALL MTPARM(ZXX,ZXY,ZYXR,ZYYR,PSWIFR,PANIZR,PSKEWR)
C
      D=SHZ(1,J)*SHY(2,J)-SHZ(2,J)*SHY(1,J)
      WX=D/DET
C
      D=SHZ(2,J)*SHX(1,J)-SHZ(1,J)*SHX(2,J)
      WY=D/DET
C
      TR=CMPLX(REAL(WX),REAL(WY))
      TRA=CABS(TR)
      TRP=PHASE(TR)
      TI=CMPLX(AIMAG(WX),AIMAG(WY))
      TIA=CABS(TI)
      TIP=PHASE(TI)
C
      PERL=ALOG10(PERIOD)
      PERL2=0.5*PERL
C
      IF(J.EQ.1.OR.J.EQ.N)GOTO 2
      IF(IRG(J-1,MD).EQ.IRG(J,MD))GOTO 2
C
      WRITE(11,2001)J,SYY(J),PERIOD,PERL,PERL2
      WRITE(11,2000)ZXX,ZXY,ZYXL,ZYYL,RAXY,RAXYL,PAXY,RAYXL,RAYXLL,
     /             PAYXL,PSWIFL,PANIZL,PSKEWL,WX,WY,TRA,TRP,TIA,TIP
c      write(7,7002)syy(j),raxy,paxy,rayxl,payxl
c 7002 format(5f12.3)
c      raxx=ommii*cabs(zxx)*cabs(zxx)
c      paxx=phase(zxx)
c      rayyl=ommii*cabs(zyyl)*cabs(zyyl)
c      payyl=phase(zyyl)
c      write(7,7002)syy(j),raxx,paxx,rayyl,payyl
      WRITE(11,2002)J,SYY(J),PERIOD,PERL,PERL2
      WRITE(11,2000)ZXX,ZXY,ZYXR,ZYYR,RAXY,RAXYL,PAXY,RAYXR,RAYXRL,
     /             PAYXR,PSWIFR,PANIZR,PSKEWR,WX,WY,TRA,TRP,TIA,TIP
c      write(7,7002)syy(j),raxy,paxy,rayxr,payxr
c      rayyr=ommii*cabs(zyyr)*cabs(zyyr)
c      payyr=phase(zyyr)
c      write(7,7002)syy(j),raxx,paxx,rayyr,payyr
      GOTO 1
C
    2 WRITE(11,2003)J,SYY(J),PERIOD,PERL,PERL2
      WRITE(11,2000)ZXX,ZXY,ZYXL,ZYYL,RAXY,RAXYL,PAXY,RAYXL,RAYXLL,
     /             PAYXL,PSWIFL,PANIZL,PSKEWL,WX,WY,TRA,TRP,TIA,TIP
c      write(7,7001)j,zxx,zxy,zyxl,zyyl,wx,wy
 7001 format(i10/8f10.6/4f10.6)
c      write(7,7002)syy(j),raxy,paxy,rayxl,payxl
c      raxx=ommii*cabs(zxx)*cabs(zxx)
c      paxx=phase(zxx)
c      rayyl=ommii*cabs(zyyl)*cabs(zyyl)
c      payyl=phase(zyyl)
c      write(7,7002)syy(j),raxx,paxx,rayyl,payyl
C
    1 CONTINUE
C
      RETURN
C
 2000 FORMAT(2X,8E11.3/2X,2F11.4,F11.1,2F11.4,F11.1/4X,F9.2,2F11.4/2X,5F
     /11.4,F11.1,F11.4,F11.1/)
 2001 FORMAT(2X,'POINT=',I3,'-',4X,'DISTANCE=',F10.2,5X,'PERIOD=',F10.2,
     /2X,'LOG(PER)=',F7.3,2X,'LOG(SQRT(PER))=',F7.3/2X,'----------------
     /------------------------------------------------------------------
     /---------------')
 2002 FORMAT(2X,'POINT=',I3,'+',4X,'DISTANCE=',F10.2,5X,'PERIOD=',F10.2,
     /2X,'LOG(PER)=',F7.3,2X,'LOG(SQRT(PER))=',F7.3/2X,'----------------
     /------------------------------------------------------------------
     /---------------')
 2003 FORMAT(2X,'POINT=',I3,5X,'DISTANCE=',F10.2,5X,'PERIOD=',F10.2,2X,'
     /LOG(PER)=',F7.3,2X,'LOG(SQRT(PER))=',F7.3/2X,'--------------------
     /------------------------------------------------------------------
     /-----------')
 2004 FORMAT(1H1,1X,'SURFACE VALUES OF THE GEOELECTRIC FUNCTIONS')
 2005 FORMAT(2X,'=======================================================
     /==========================================')
 2006 FORMAT(5X,'RE(ZXX)',4X,'IM(ZXX)',4X,'RE(ZXY)',4X,'IM(ZXY)',4X,'RE(
     /ZYX)',4X,'IM(ZYX)',4X,'RE(ZYY)',4X,'IM(ZYY)'/6X,'ROAXY',3X,'LOG(RO
     /AXY)',2X,'PH(ZXY)',6X,'ROAYX',3X,'LOG(ROAYX)',2X,'PH(ZYX)'/5X,'SWI
     /FT',7X,'ANIZ',7X,'SKEW'/5X,'RE(WX)',5X,'IM(WX)',5X,'RE(WY)',5X,'IM
     /(WY)',7X,'TR',7X,'PH(TR)',7X,'TI',7X,'PH(TI)')
 2007 FORMAT(/)
C
      END
      SUBROUTINE DER3PO(Y1,Y2,Y3,H12,H23,IY,DER)
C     ==================================================================
C
C-----------------------------------------------------------------------
      COMPLEX Y1,Y2,Y3,D1,D2,DER
C-----------------------------------------------------------------------
C
      D1=(Y2-Y1)/H12
      D2=(Y3-Y2)/H23
      D12=H12/(H12+H23)
      D21=H23/(H12+H23)
C
      IF(IY-2)1,2,3
C
    1 DER=(D12+D12+D21)*D1-D12*D2
      RETURN
C
    2 DER=D21*D1+D12*D2
      RETURN
C
    3 DER=-D21*D1+(D12+D21+D21)*D2
      RETURN
C
      END
      FUNCTION PHASE(Z)
C     ==================================================================
C
C-----------------------------------------------------------------------
      PARAMETER(PII=3.141592653589793)
      COMPLEX Z
C-----------------------------------------------------------------------
C
      ZR=REAL(Z)
      ZI=AIMAG(Z)
C
      IF(ZR)1,2,3
C
    1 ZNAZI=1.
      IF(ZI.LT.0.)ZNAZI=-1.
      PHASE=ATAN(ZI/ZR)+PII*ZNAZI
      GOTO 4
C
    2 IF(ZI)10,20,30
   10 PHASE=-PII/2.
      GOTO 4
   20 PHASE=0.
      GOTO 4
   30 PHASE=PII/2.
      GOTO 4
C
    3 PHASE=ATAN(ZI/ZR)
    4 PHASE=180.*PHASE/PII
C
      RETURN
C
      END
      SUBROUTINE GAUSSR
C     ==================================================================
C
C-----------------------------------------------------------------------
      PARAMETER(NMAX=151,MMAX=59)
      PARAMETER(NEHMAX=2*(NMAX-2)*(MMAX-2),MEHMAX=2*(MMAX-2)+3)
      COMPLEX*16 APOM,BPOM
      COMPLEX*16 C
      REAL*8 AIK,ANK
C
      COMMON /NHS/ NEH
      COMMON /RHS/ BPOM(NEHMAX)
      COMMON /MAT/ APOM(NEHMAX,MEHMAX)
      COMMON /MES/ N,M,MD,N1,N2,M1,M2,SY(NMAX-1),SZ(MMAX-1)
C-----------------------------------------------------------------------
C
      WRITE(*,'(a12)')' GAUSS in ->'
      NK=NEH
      M1STOR=M1
      MH2=M1-MD
      MEH=2*M2-MD+4
      M1=MEH
C
      NK1=NK-1
      NKM=NK-M1+2
      MG=M1
      NKM1=NKM-1
        N9=NK/10
        J9M=-1
      DO 210 I=1,NK1
      I1=I-1
      IF(I.GE.NKM)MG=MG-1
      ME=2
      DO 220 K=2,MG
      AIK=CDABS(APOM(I,K))
      IF(AIK.EQ.0.)GOTO 220
      C=APOM(I,K)/APOM(I,1)
      IK=I1+K
      J=0
      DO 230 L=ME,MG
      J=J+1
  230 APOM(IK,J)=APOM(IK,J)-C*APOM(I,L)
      BPOM(IK)=BPOM(IK)-C*BPOM(I)
  220 ME=ME+1
        J9=I/N9
        IF(J9.GT.J9M)THEN
          WRITE(*,'(i3)')J9
          J9M=J9
        ENDIF
  210 CONTINUE
C
      NE=NK+1
  310 NE=NE-1
      IF(NE)350,350,320
  320 L=NE
      DO 340 K=2,M1
      L=L+1
      IF(L.GT.NK)GOTO 360
      ANK=CDABS(APOM(NE,K))
      IF(ANK)330,340,330
  330 BPOM(NE)=BPOM(NE)-APOM(NE,K)*BPOM(L)
  340 CONTINUE
  360 BPOM(NE)=BPOM(NE)/APOM(NE,1)
      GOTO 310
  350 CONTINUE
C
      M1=M1STOR
C
      WRITE(*,'(a13)')' -> GAUSS out'
C
      RETURN
      END
      SUBROUTINE MTPARM(ZXX,ZXY,ZYX,ZYY,PSWIF,PANIZ,PSKEW)
C     ==================================================================
C
C-----------------------------------------------------------------------
      PARAMETER(PI=3.141592653589793)
      COMPLEX ZXX,ZXY,ZYX,ZYY,D1,D2,S1,S2
C-----------------------------------------------------------------------
C
      D1=ZXX-ZYY
      D2=ZXY-ZYX
      S1=ZXX+ZYY
      S2=ZXY+ZYX
C
      IF(CABS(D1).EQ.CABS(S2))THEN
        PSWIF=0.
        PANIZ=AMAX1(CABS(ZXY)/CABS(ZYX),CABS(ZYX)/CABS(ZXY))
      ELSE
        PSWIF=0.25*ATAN(2.*REAL(D1*CONJG(S2))/
     /                  (CABS(D1)**2.-CABS(S2)**2.))
        SWIF2=2.*PSWIF
        CO2=COS(SWIF2)
        SI2=SIN(SWIF2)
        S2SWI1=CABS(S2*CO2-D1*SI2)
        S2SWI2=CABS(-S2*SI2+D1*CO2)
        IF(S2SWI2.GT.S2SWI1)PSWIF=PSWIF+0.25*PI
        IF(PSWIF.GT.0.5*PI)PSWIF=PSWIF-PI
        IF(PSWIF.LE.-0.5*PI)PSWIF=PSWIF+PI
        SWIF2=2.*PSWIF
        CO2=COS(SWIF2)
        SI2=SIN(SWIF2)
        S2SWI1=CABS(D2+S2*CO2-D1*SI2)
        S2SWI2=CABS(-D2+S2*CO2-D1*SI2)
        PANIZ=AMAX1(S2SWI1/S2SWI2,S2SWI2/S2SWI1) 
      ENDIF
      PSWIF=180.*PSWIF/PI
C
      PSKEW=CABS(S1)/CABS(D2)
C
      RETURN
      END
