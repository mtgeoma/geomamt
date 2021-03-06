      DOUBLE PRECISION FUNCTION S07AAF(X,IFAIL)
C     MARK 5A REVISED  -  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     TAN(X)
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S07AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 FL1, FL2, SGN, T, TBPI, THETA,
     *                                 XN, Y
      INTEGER                          N
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MOD, SIGN, DBLE, INT, AINT
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C08   DATA FL1,FL2,TBPI
C08  A/1.0D5,1.0D-6,6.3661977D-1/
C12   DATA FL1,FL2,TBPI
C12  A/1.0D7,1.0D-9,6.36619772368D-1/
C14   DATA FL1,FL2,TBPI
C14  A/1.0D11,1.0D-12,6.3661977236758D-1/
      DATA FL1,FL2,TBPI
     A/1.0D13,1.0D-14,6.366197723675813D-1/
C18   DATA FL1,FL2,TBPI
C18  A/1.0D15,1.0D-16,6.36619772367581343D-1/
C     .. Executable Statements ..
C
C     ERROR 1 TEST
      IF (ABS(X).GT.FL1) GO TO 20
C
C     RANGE REDUCTION
      SGN = 1.0D0
      IF (X.LT.0.0D0) SGN = -SGN
      XN = AINT(X*TBPI+SGN*0.5D0)
      THETA = X - XN/TBPI
      N = INT(MOD(ABS(XN),2.0D0))
C
C     ERROR 2 TEST
      IF (N.EQ.1 .AND. ABS(THETA).LT.ABS(X)*FL2) GO TO 40
C
C     EXPANSION ARGUMENT
      T = (TBPI*THETA)**2*8.0D0 - 1.0D0
C
C      * EXPANSION (0007) *
C
C     EXPANSION (0007) EVALUATED AS Y(T)  --PRECISION 08E
C08   Y = (((((((+1.1781590D-6)*T+8.2048199D-6)*T+5.5077467D-5)
C08  *    *T+3.8561954D-4)*T+2.7010438D-3)*T+1.8924537D-2)
C08  *    *T+1.3386247D-1)*T + 1.1173014D+0
C
C     EXPANSION (0007) EVALUATED AS Y(T)  --PRECISION 12E
C12   Y = ((((((((((+3.48825813626D-9)*T+2.42925841226D-8)
C12  *    *T+1.60455379111D-7)*T+1.12350071991D-6)
C12  *    *T+7.87409839245D-6)*T+5.51184605550D-5)
C12  *    *T+3.85828287625D-4)*T+2.70103244160D-3)
C12  *    *T+1.89244955238D-2)*T+1.33862473672D-1)*T +
C12  *     1.11730141051D+0
C
C     EXPANSION (0007) EVALUATED AS Y(T)  --PRECISION 14E
C14   Y = ((((((((((((+7.1924759283372D-11)*T+5.0089133229408D-10)
C14  *    *T+3.2724838584129D-9)*T+2.2915132958756D-8)
C14  *    *T+1.6069812517338D-7)*T+1.1248781710765D-6)
C14  *    *T+7.8739725241234D-6)*T+5.5117857920090D-5)
C14  *    *T+3.8582831712490D-4)*T+2.7010325492100D-3)
C14  *    *T+1.8924495521240D-2)*T+1.3386247366625D-1)*T +
C14  *     1.1173014105142D+0
C
C     EXPANSION (0007) EVALUATED AS Y(T)  --PRECISION 16E
      Y = ((((((((((((((+1.483024132942717D-12)
     *    *T+1.032793075951471D-11)*T+6.673417481807284D-11)
     *    *T+4.673255573256614D-10)*T+3.279620912052723D-9)
     *    *T+2.295709017746642D-8)*T+1.606932590004467D-7)
     *    *T+1.124852996745234D-6)*T+7.873974227283909D-6)
     *    *T+5.511786526260323D-5)*T+3.858283168410358D-4)
     *    *T+2.701032548292210D-3)*T+1.892449552125797D-2)
     *    *T+1.338624736662861D-1)*T + 1.117301410514158D+0
C
C     EXPANSION (0007) EVALUATED AS Y(T)  --PRECISION 18E
C18   Y = (((((((((((((((+2.12952684337496064D-13)
C18  *    *T+1.48302413294271693D-12)*T+9.52935819324910387D-12)
C18  *    *T+6.67341748180728381D-11)*T+4.68523416175059801D-10)
C18  *    *T+3.27962091205272343D-9)*T+2.29561751464009084D-8)
C18  *    *T+1.60693259000446697D-7)*T+1.12485337107612437D-6)
C18  *    *T+7.87397422728390854D-6)*T+5.51178651839937393D-5)
C18  *    *T+3.85828316841035795D-4)*T+2.70103254829948866D-3)
C18  *    *T+1.89244955212579706D-2)*T+1.33862473666285911D-1)*T +
C18  *     1.11730141051415794D+0
C
C
      S07AAF = THETA*Y
      IFAIL = 0
      IF (N.EQ.1) S07AAF = -1.0D0/S07AAF
      RETURN
C
C     ERROR 1 EXIT
   20 S07AAF = 0.0D0
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
C
C     ERROR 2 EXIT
   40 IF (THETA.EQ.0.0D0) THETA = -1.0D0
      S07AAF = -SIGN(1.0D0/(X*FL2),THETA)
      IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      RETURN
C
      END
