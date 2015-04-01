       subroutine sdmsort(z)
       complex z(*),ztemp(15)
       integer iswitch(15),iconj(8)
       data iswitch/10,14,15,9,13,6,7,11,4,1,8,12,5,2,3/
       data iconj/4,5,7,8,9,11,12,13/

       do 10 i = 1,15
10     ztemp(i) = z(i)
       do 20 i = 1,15
20     z(i) = ztemp(iswitch(i))
       do 30 i = 1,9
30     z(iconj(i)) = conjg(z(iconj(i)))
       return
       end
