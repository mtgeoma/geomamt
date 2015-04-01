      SUBROUTINE F02XUV(SHIFT,M,N,D,E,C,S,WANTLT,CL,SL,WANTRT,CR,SR)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  F02XUV performs a QR step on a bidiagonal matrix A.
C
C  When  SHIFT = 'N' or 'n'  ( Non-zero ),  then the parameters  c and s
C  must define the  initial right hand  plane rotation  to be used,  but
C  when  SHIFT = 'Z' or 'z'  ( Zero ), then the initial right hand plane
C  is determined by F02XUV so that
C
C     ( d( m )  0 ) := ( d( m )  e( m ) )*( c  -s )
C                                         ( s   c )
C
C  and c and s need not be supplied.
C
C  In either case,  plane rotations are then performed  alternately from
C  the  left and right  in order to  recover  the  bidiagonal form.  The
C  bidiagonal matrix  must be  supplied  with the  diagonal elements  in
C  d( m ), d( m + 1 ), ..., d( n )  and the  super-diagonal elements  in
C  e( m ), e( m + 1 ), ..., e( n - 1 ). When  m  is greater than  unity,
C  the assumption is that   e( m - 1 ) = 0.0.
C
C  Thus F02XUV performs the transformation
C
C     A := Q*A*P',
C
C  where  Q and P  are  orthogonal matrices.  P and Q  each consist of a
C  sequence of plane rotations
C
C     P = P( n - 1 )*...*P( m + 1 )*P( m )
C
C     Q = Q( n - 1 )*...*Q( m + 1 )*Q( m ),
C
C  where Q( k ) and P( k ) are each plane rotations for the ( k, k + 1 )
C  plane. The two by two part of the plane rotation matrix  Q( k )  will
C  be of the form
C
C     Q2 = (  cl( k )  sl( k ) ),
C          ( -sl( k )  cl( k ) )
C
C  and  if  WANTLT  is supplied as  true  then  cl( k ) and sl( k )  are
C  returned in the  corresponding elements of the arrays  CL and SL.  If
C  WANTLT  is supplied as false then  CL and SL  are not referenced. The
C  two by two  part of the plane rotation matrix  P( k )  will be of the
C  form
C
C     P2 = (  cr( k )  sr( k ) ),
C          ( -sr( k )  cr( k ) )
C
C  and  if  WANTRT  is supplied as  true  then  cr( k ) and sr( k )  are
C  returned in the corresponding elements of the arrays  CR and SR. Note
C  that  cr( m ) = c  and  sr( m ) = s. If  WANTRT  is supplied as false
C  then  CR and SR  are not referenced.
C
C  If  m.le.0  or  n.le.m  then an immediate return is effected.
C
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 10-August-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, S
      INTEGER           M, N
      LOGICAL           WANTLT, WANTRT
      CHARACTER*1       SHIFT
C     .. Array Arguments ..
      DOUBLE PRECISION  CL(*), CR(*), D(*), E(*), SL(*), SR(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CO, CS, EI, SI, SN, TEMP
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Executable Statements ..
      IF ((M.GT.0) .AND. (M.LT.N)) THEN
         IF ((SHIFT.EQ.'N') .OR. (SHIFT.EQ.'n')) THEN
            IF (WANTRT) THEN
               CR(M) = C
               SR(M) = S
            END IF
            TEMP = C*D(M) + S*E(M)
            E(M) = C*E(M) - S*D(M)
            D(M) = TEMP
            TEMP = S*D(M+1)
            D(M+1) = C*D(M+1)
            DO 20 I = M, N - 2
               CALL F06BAF(D(I),TEMP,CS,SN)
               IF (WANTLT) THEN
                  CL(I) = CS
                  SL(I) = SN
               END IF
               TEMP = CS*E(I) + SN*D(I+1)
               D(I+1) = CS*D(I+1) - SN*E(I)
               E(I) = TEMP
               TEMP = SN*E(I+1)
               E(I+1) = CS*E(I+1)
               CALL F06BAF(E(I),TEMP,CS,SN)
               IF (WANTRT) THEN
                  CR(I+1) = CS
                  SR(I+1) = SN
               END IF
               TEMP = CS*D(I+1) + SN*E(I+1)
               E(I+1) = CS*E(I+1) - SN*D(I+1)
               D(I+1) = TEMP
               TEMP = SN*D(I+2)
               D(I+2) = CS*D(I+2)
   20       CONTINUE
            CALL F06BAF(D(N-1),TEMP,CS,SN)
            IF (WANTLT) THEN
               CL(N-1) = CS
               SL(N-1) = SN
            END IF
            TEMP = CS*E(N-1) + SN*D(N)
            D(N) = CS*D(N) - SN*E(N-1)
            E(N-1) = TEMP
         ELSE
            CALL F06BAF(D(M),E(M),C,S)
            IF (WANTRT) THEN
               CR(M) = C
               SR(M) = S
            END IF
            TEMP = S*D(M+1)
            D(M+1) = C*D(M+1)
            DO 40 I = M, N - 2
               CALL F06BAF(D(I),TEMP,CO,SI)
               IF (WANTLT) THEN
                  CL(I) = CO
                  SL(I) = SI
               END IF
               EI = D(I+1)
               TEMP = E(I+1)
               CALL F06BAF(EI,TEMP,CS,SN)
               IF (WANTRT) THEN
                  CR(I+1) = CS
                  SR(I+1) = SN
               END IF
               E(I) = SI*EI
               D(I+1) = CO*EI
               TEMP = SN*D(I+2)
               D(I+2) = CS*D(I+2)
   40       CONTINUE
            CALL F06BAF(D(N-1),TEMP,CS,SN)
            IF (WANTLT) THEN
               CL(N-1) = CS
               SL(N-1) = SN
            END IF
            E(N-1) = SN*D(N)
            D(N) = CS*D(N)
         END IF
      END IF
      RETURN
C
C     End of F02XUV. ( SBIQR1 )
C
      END
