      SUBROUTINE F06QMF( UPLO, PIVOT, DIRECT, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        UPLO, PIVOT, DIRECT
      INTEGER            N, K1, K2, LDA
C     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * ), A( LDA, * )
C     ..
C
C  F06QMF performs the transformation
C
C     A := P*A*P',
C
C  where A is an n by n symmetric matrix and  P is an orthogonal matrix,
C  consisting of a sequence of plane rotations, applied in planes  k1 to
C  k2,  determined  by the  parameters  PIVOT  and  DIRECT  as  follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ). The two by two plane rotation part of P( k ), R( k ),
C  is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ),
C              ( -s( k )  c( k )   )
C
C  The parameter  UPLO  determines whether the upper or lower triangular
C  part of A is referenced as follows.
C
C     When  UPLO = 'U' or 'u'  then the leading  n by n upper triangular
C     part of  A must contain the upper triangular part of the symmetric
C     matrix  and the  strictly  lower  triangular  part  of  A  is  not
C     referenced.
C
C     When  UPLO = 'L' or 'l'  then the leading  n by n lower triangular
C     part of  A must contain the lower triangular part of the symmetric
C     matrix  and the  strictly  upper  triangular  part  of  A  is  not
C     referenced.
C
C  If  n or k1 are less than unity, or  k2 is not greater than k1, or k2
C  is greater than n, then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-December-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C     .. External Subroutines ..
      EXTERNAL           F06BHF
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C     .. Local Scalars ..
      INTEGER            I, I1, I2, J
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( UPLO.EQ.'L' ).OR.( UPLO.EQ.'l' ) )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
C
C              First use  F06BHF  to perform the  two by two  similarity
C              transformations  to those  elements  that get  changed by
C              P( j ), j = k1, ..., k2 - 1,  both from the  left and the
C              right.  See the comment in  F06BHF  for an explanation of
C              the calling sequence.
C
C
C              Next apply the right hand rotations, columnwise.
C
               DO 20 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     CALL F06BHF( A( J + 1, J + 1 ), A( J + 1, J ),
     $                            A( J, J ), CTEMP, -STEMP )
                     DO 10 I = J + 2, N
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
C
C              Now apply the left hand rotations, columnwise.
C
               DO 40 J = 1, K2 - 2
                  I1 = MAX( K1, J + 1 )
                  AIJ = A( I1, J )
                  DO 30 I = I1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   30             CONTINUE
                  A( K2, J ) = AIJ
   40          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
C
C              First apply the left hand rotations, columnwise.
C
               DO 60 J = K2 - 2, 1, -1
                  AIJ = A( K2, J )
                  I2 = MAX( K1, J + 1 )
                  DO 50 I = K2 - 1, I2, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( I2, J ) = AIJ
   60          CONTINUE
C
C              Next use  F06BHF  to perform  the  two by two  similarity
C              transformations  to those  elements  that get  changed by
C              P( j ), j = k2 - 1, ..., k1,  both from the  left and the
C              right.
C
               DO 70 J = K2 - 1, K1, -1
                  CALL F06BHF( A( J + 1, J + 1 ), A( J + 1, J ),
     $                         A( J, J ), C( J ), -S( J ) )
   70          CONTINUE
C
C              Now  apply  the  right hand rotations,  again columnwise.
C
               DO 90 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 80 I = N, J + 2, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   80                CONTINUE
                  END IF
   90          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 110 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
C
C                    First  use   F06BHF  to  perform  the  two  by  two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( J, J ), A( J, K1 ), A( K1, K1 ),
     $                            CTEMP, -STEMP )
C
C                    Next  apply the  right hand rotations,  columnwise.
C
                     DO 100 I = J + 1, N
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  100                CONTINUE
                  END IF
  110          CONTINUE
C
C              Now apply the left hand rotations, columnwise.
C
               DO 130 J = 1, K1 - 1
                  TEMP = A( K1, J )
                  DO 120 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
  120             CONTINUE
                  A( K1, J ) = TEMP
  130          CONTINUE
               DO 150 J = K1 + 1, K2 - 1
                  TEMP = A( J, K1 )
                  DO 140 I = J, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
  140             CONTINUE
                  A( J, K1 ) = TEMP
  150          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 180 J = K2, K1 + 1, -1
C
C                 First  apply   the   immediately  required  left  hand
C                 rotations,  columnwise.
C
                  TEMP = A( J, K1 )
                  DO 160 I = K2 - 1, J, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = -S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP + S( I )*AIJ
  160             CONTINUE
                  A( J, K1 ) = TEMP
C
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
C
C                    Next  use   F06BHF  to  perform   the  two  by  two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( J, J ), A( J, K1 ), A( K1, K1 ),
     $                            CTEMP, -STEMP )
C
C                    Next  apply the  right hand rotations,  columnwise.
C
                     DO 170 I = J + 1, N
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
C
C              Now apply the remaining left hand rotations,  columnwise.
C
               DO 200 J = 1, K1 - 1
                  TEMP = A( K1, J )
                  DO 190 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
  190             CONTINUE
                  A( K1, J ) = TEMP
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
C
C           N.B.  When the bottom element is  the pivot,  and the matrix
C           is stored in the lower triangle of A, some rotations have to
C           be performed rowwise.  This may be disadvantageous on  paged
C           machines.
C
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 230 J = K1, K2 - 1
C
C                 First  apply   the   immediately  required   left hand
C                 rotations.   These  have  to  be  performed   rowwise.
C
                  TEMP = A( K2, J )
                  DO 210 I = K1, J - 1
                     AIJ = A( J, I )
                     A( J, I ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  210             CONTINUE
                  A( K2, J ) = TEMP
C
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
C
C                    Next  use   F06BHF   to  perform  the  two  by  two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( K2, K2 ), A( K2, J ), A( J, J ),
     $                            CTEMP, -STEMP )
C
C                    Next  apply the  right hand rotations,  columnwise.
C
                     DO 220 I = K2 + 1, N
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP + STEMP*A( I, K2 )
                        A( I, K2 ) = -STEMP*TEMP + CTEMP*A( I, K2 )
  220                CONTINUE
                  END IF
  230          CONTINUE
C
C              Now apply the remaining left hand rotations,  columnwise.
C
               DO 250 J = 1, K2 - 1
                  TEMP = A( K2, J )
                  I1 = MAX( K1, J + 1 )
                  DO 240 I = I1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  240             CONTINUE
                  A( K2, J ) = TEMP
  250          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 280 J = K2 - 1, K1, -1
C
C                 First  apply   the   immediately  required  left  hand
C                 rotations,  columnwise.
C
                  TEMP = A( K2, J )
                  DO 260 I = K2 - 1, J + 1, -1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  260             CONTINUE
                  A( K2, J ) = TEMP
C
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
C
C                    Next  use   F06BHF   to  perform  the  two  by  two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( K2, K2 ), A( K2, J ), A( J, J ),
     $                            CTEMP, -STEMP )
C
C                    Next  apply the  right hand rotations,  columnwise.
C
                     DO 270 I = K2 + 1, N
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP + STEMP*A( I, K2 )
                        A( I, K2 ) = -STEMP*TEMP + CTEMP*A( I, K2 )
  270                CONTINUE
                  END IF
  280          CONTINUE
C
C              Now apply the remaining left hand rotations,  columnwise.
C
               DO 300 J = 1, K1 - 1
                  TEMP = A( K2, J )
                  DO 290 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  290             CONTINUE
                  A( K2, J ) = TEMP
  300          CONTINUE
C
C              Now apply the remaining right hand rotations.  These have
C              to be performed rowwise.
C
               DO 320 J = K2 - 2, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 310 I = K2 - 1, J + 1, -1
                        TEMP = A( K2, I )
                        A( K2, I ) = -STEMP*A( I, J ) + CTEMP*TEMP
                        A( I, J ) = CTEMP*A( I, J ) + STEMP*TEMP
  310                CONTINUE
                  END IF
  320          CONTINUE
            END IF
         END IF
      ELSE IF( ( UPLO.EQ.'U' ).OR.( UPLO.EQ.'u' ) )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
C
C              First apply the left hand rotations, columnwise.
C
               DO 340 J = K1 + 2, N
                  AIJ = A( K1, J )
                  I2 = MIN( K2 - 1, J - 2 )
                  DO 330 I = K1, I2
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
  330             CONTINUE
                  A( I2 + 1, J ) = AIJ
  340          CONTINUE
C
C              Now  apply the  right hand rotations,  again  columnwise.
C
               DO 360 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 350 I = 1, J - 1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  350                CONTINUE
                  END IF
C
C                 Next use  F06BHF  to perform the two by two similarity
C                 transformations  to those elements that get changed by
C                 P( j ), j = k1, ..., k2 - 1,  both  from the  left and
C                 the right.
C
                  CALL F06BHF( A( J, J ), A( J, J + 1 ),
     $                         A( J + 1, J + 1 ), C( J ), S( J ) )
  360          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
C
C              First apply the right hand rotations, columnwise.
C
               DO 380 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 370 I = 1, J - 1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  370                CONTINUE
                  END IF
  380          CONTINUE
C
C              Next use  F06BHF  to  perform the  two by two  similarity
C              transformations  to those  elements  that get  changed by
C              P( j ), j = k2 - 1, ..., k1,  both from the  left and the
C              right.
C
               DO 390 J = K2 - 1, K1, -1
                  CALL F06BHF( A( J, J ), A( J, J + 1 ),
     $                         A( J + 1, J + 1 ), C( J ), S( J ) )
  390          CONTINUE
C
C              Now apply the left hand rotations, columnwise.
C
               DO 410 J = N, K1 + 2, -1
                  I1 = MIN( K2 - 1, J - 2 )
                  AIJ = A( I1 + 1, J )
                  DO 400 I = I1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
  400             CONTINUE
                  A( K1, J ) = AIJ
  410          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
C
C           N.B.  When the top element is the  pivot,  and the matrix is
C           stored in the upper triangle of A, some rotations have to be
C           performed  rowwise.  This  may be  disadvantageous on  paged
C           machines.
C
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 440 J = K1 + 1, K2
C
C                 First  apply   the   immediately  required   left hand
C                 rotations,  columnwise.
C
                  TEMP = A( K1, J )
                  DO 420 I = K1, J - 2
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = -S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP + S( I )*AIJ
  420             CONTINUE
                  A( K1, J ) = TEMP
C
C                 Apply the  right hand rotations that can be  performed
C                 columnwise, and use  F06BHF  to perform the two by two
C                 similarity transformations  to those elements that get
C                 changed by  P( j ), j = k2 - 1, ..., k1  both from the
C                 left and the right.
C
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 430 I = 1, K1 - 1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  430                CONTINUE
C
                     CALL F06BHF( A( K1, K1 ), A( K1, J ), A( J, J ),
     $                            CTEMP, STEMP )
C
                  END IF
  440          CONTINUE
C
C              Now apply the remaining left hand rotations,  columnwise.
C
               DO 460 J = K2 + 1, N
                  TEMP = A( K1, J )
                  DO 450 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
  450             CONTINUE
                  A( K1, J ) = TEMP
  460          CONTINUE
C
C              Apply  the remaining  right hand rotations  that  must be
C              performed rowwise.
C
               DO 480 J = K1 + 2, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 470 I = K1 + 1, J - 1
                        TEMP = A( K1, I )
                        A( K1, I ) = CTEMP*TEMP + STEMP*A( I, J )
                        A( I, J ) = -STEMP*TEMP + CTEMP*A( I, J )
  470                CONTINUE
                  END IF
  480          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 520 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
C
C                    First  use   F06BHF  to  perform  the  two  by  two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( K1, K1 ), A( K1, J ), A( J, J ),
     $                            CTEMP, STEMP )
C
C                    Then apply the right hand rotations that have to be
C                    performed rowwise.
C
                     DO 490 I = J - 1, K1 + 1, -1
                        TEMP = A( K1, I )
                        A( K1, I ) = CTEMP*TEMP + STEMP*A( I, J )
                        A( I, J ) = -STEMP*TEMP + CTEMP*A( I, J )
  490                CONTINUE
C
C                    Apply    the   remaining    right  hand  rotations,
C                    columnwise.
C
                     DO 500 I = 1, K1 - 1
                        TEMP = A( I, J )
                        A( I, J ) = -STEMP*A( I, K1 ) + CTEMP*TEMP
                        A( I, K1 ) = CTEMP*A( I, K1 ) + STEMP*TEMP
  500                CONTINUE
                  END IF
C
C                 Apply  some  of the  left hand rotations,  columnwise.
C
                  TEMP = A( K1, J )
                  DO 510 I = J - 2, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = -S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP + S( I )*AIJ
  510             CONTINUE
                  A( K1, J ) = TEMP
  520          CONTINUE
C
C              Now apply the remaining left hand rotations,  columnwise.
C
               DO 540 J = K2 + 1, N
                  TEMP = A( K1, J )
                  DO 530 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = -S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP + S( I )*AIJ
  530             CONTINUE
                  A( K1, J ) = TEMP
  540          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 570 J = K1, K2 - 1
C
C                 First  apply   the  immediately  required   left  hand
C                 rotations, columnwise.
C
                  TEMP = A( J, K2 )
                  DO 550 I = K1, J - 1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  550             CONTINUE
                  A( J, K2 ) = TEMP
C
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
C
C                    Next  use   F06BHF   to  perform  the  two  by  two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( J, J ), A( J, K2 ), A( K2, K2 ),
     $                            CTEMP, STEMP )
C
C                    Next  apply the  right hand rotations,  columnwise.
C
                     DO 560 I = J - 1, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP + STEMP*A( I, K2 )
                        A( I, K2 ) = -STEMP*TEMP + CTEMP*A( I, K2 )
  560                CONTINUE
                  END IF
  570          CONTINUE
C
C              Now apply the remaining left hand rotations,  columnwise.
C
               DO 590 J = K2 + 1, N
                  TEMP = A( K2, J )
                  DO 580 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  580             CONTINUE
                  A( K2, J ) = TEMP
  590          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 610 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
C
C                    First  use   F06BHF  to  perform  the  two  by two
C                    similarity transformation  to those  elements  that
C                    get  changed  both  from  the  left  and the right.
C
                     CALL F06BHF( A( J, J ), A( J, K2 ), A( K2, K2 ),
     $                            CTEMP, STEMP )
C
C                    Next  apply the  right hand rotations,  columnwise.
C
                     DO 600 I = J - 1, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP + STEMP*A( I, K2 )
                        A( I, K2 ) = -STEMP*TEMP + CTEMP*A( I, K2 )
  600                CONTINUE
                  END IF
  610          CONTINUE
C
C              Now apply the left hand rotations, columnwise.
C
               DO 630 J = N, K2 + 1, -1
                  TEMP = A( K2, J )
                  DO 620 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  620             CONTINUE
                  A( K2, J ) = TEMP
  630          CONTINUE
               DO 650 J = K2 - 1, K1 + 1, -1
                  TEMP = A( J, K2 )
                  DO 640 I = J - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = C( I )*AIJ + S( I )*TEMP
                     TEMP = -S( I )*AIJ + C( I )*TEMP
  640             CONTINUE
                  A( J, K2 ) = TEMP
  650          CONTINUE
            END IF
         END IF
      END IF
      RETURN
C
C     End of F06QMF. ( SSYSRC )
C
      END
