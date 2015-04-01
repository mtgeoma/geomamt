        subroutine shiftdown (string)
        character*(*) string
        integer table(256), ns
c converts uppercase letters in string to lowercase

C       character table (256) /
        data table /
     1  1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
     2 17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
     3 33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,
     4 49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,
     5  97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,112,
     6 113,114,115,116,117,118,119,120,121,122, 91, 92, 93, 94, 95, 96,
     7  97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,112,
     8 113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,
     9 128*0/

        include 'ctrlblk.inc'

        if( iprint.ge.30 ) then
           write(*,*)'shiftdown: entered with string=',string
        end if

        call trunc( string, ns )
        do 1 i = 1, ns
           string(i:i) = char(table(ichar(string(i:i))))
    1   continue

        return
        end
