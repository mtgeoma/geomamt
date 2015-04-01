c*************************************************
c
c  this subroutine multiplies matrix mat1 by mat2
c  with product in mat3
c
c***********************************************


       subroutine mat_multiply(n,mat1,mat2,mat3)
       complex*16 mat1(n,n),mat2(n,n),mat3(n,n)

       do j = 1, n
         do i = 1, n
         mat3(i,j) = dcmplx(0.,0.)
           do k = 1, n
             mat3(i,j)=mat3(i,j) + mat1(i,k)*mat2(k,j)
           end do
         end do
       end do

       return
       end
