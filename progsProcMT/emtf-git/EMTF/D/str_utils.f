      subroutine strlen(string,i1,i2)
      character string*(*)
      do i=len(string),1,-1
         if(string(i:i).ne.' ' .and. string(i:i).ne.'	')then
            i2=i
            exit
         endif
      enddo
      do i=1,len(string)
         if(string(i:i).ne.' ' .and. string(i:i).ne.'	')then
            i1=i
            return
         endif
      enddo
      return
      end

      subroutine strtolower(string)
      character string*(*)
      integer del
      del=iachar('a') - iachar('A')
      do i=1,len(string)
         if(string(i:i).ge.'A' .and. string(i:i).le.'Z')then
            string(i:i)=achar(iachar(string(i:i))+del)
         endif
      enddo
      return
      end
