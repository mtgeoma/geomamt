c*******************************************************
c
c    estimated measured impedance creator (PRODUCT FACTORIZATION)
c
c********************************************************
       subroutine estim_imp(ztemp,gamma1,gamma2,A,B,theta)

C********************************************
C   Ztemp - estmated MEASURED IMPEDANCE TENSOR
C   C     - SCATTERING MATRIX
C   ROT   - ROTATION TO MEASURED CO-ORDINATES
C   Z2    - REGIONAL IMPEDANCE
C   ROTT  - TRANSPOSE OF ROTATION MATRIX  
C***********************************************
       complex*16 c(2,2),rot(2,2),z2(2,2),ROTT(2,2)
       complex*16 ztemp1(2,2),ztemp2(2,2),A,B,ztemp(2,2)
       real*8     gamma1,gamma2,theta

c**************************************************
c 
c  twist and scaling  -  c matrix
c
c*************************************************


c --Note negative of twist and regional azimuth are taken
c   as rotations must be counter-clockwise (Z2 is theta
c   clockwise from Zm, therefore Zm is theta counter-clockwise
c   from Z2).

        GAMMA2 = -GAMMA2
        THETA = -THETA       


        C(1,1) = DCMPLX(1.D0,0.D0)+GAMMA1*GAMMA2
        C(1,2) = DCMPLX(GAMMA1+GAMMA2,0.D0)
        C(2,1) = DCMPLX(GAMMA1-GAMMA2,0.D0)
        C(2,2) = DCMPLX(1.D0,0.D0)-GAMMA1*GAMMA2

        ROT(1,1) = DCMPLX( DCOS(THETA),0.D0)
        ROT(2,1) = DCMPLX(-DSIN(THETA),0.D0)
        ROT(1,2) = DCMPLX( DSIN(THETA),0.D0)
        ROT(2,2) = DCMPLX( DCOS(THETA),0.D0)


        CALL MAT_MULTIPLY(2,ROT,C,Ztemp1)



        Z2(1,1)=DCMPLX(0.D0,0.D0)
        Z2(2,2)=DCMPLX(0.D0,0.D0)


        Z2(1,2)= A
        Z2(2,1)= -B


        CALL MAT_MULTIPLY(2,Ztemp1,Z2,Ztemp2)


        ROTT(1,1)=ROT(1,1)
        ROTT(2,2)=ROT(2,2)
        ROTT(1,2)=ROT(2,1)
        ROTT(2,1)=ROT(1,2)

        CALL MAT_MULTIPLY(2,Ztemp2,ROTT,ztemp)

        GAMMA2 = -GAMMA2
        THETA = -THETA       

        RETURN 

        END
