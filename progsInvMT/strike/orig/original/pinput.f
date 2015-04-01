      subroutine pinput
C ***************************************************************
C *                                                             *
C *  This is a simple routine for prompting input, without the  *
C *  trouble of formatting.  A prompt text is written, then     *
C *  the current value of the variable being prompted for.      *
C *  If a carriage return is the user response, the variable    *
C *  is left at the default value (the value on entry to the    *
C *  subroutine.  If something is actually entered, the         *
C *  variable is given this value.  By calling this routine     *
C *  under different names, REAL, INTEGER, or CHARACTER         *
C *  variables can be prompted for and obtained.                *
C *                                                             *
C *        CALL RIN(PROMPT,RVALUE) will print the text in       *
C *             character variable PROMPT and read a real       *
C *             value RVALUE.                                   *
C *        CALL RIN8(PROMPT,RVALU8) will print the text in      *
C *             character variable PROMPT and read a real*8     *
C *             value RVALU8.                                   *
C *        CALL IIN(PROMPT,IVALUE) does the same except that    *
C *             it returns an integer value IVALUE.             *
C *        CALL IIN2(PROMPT,IVALU2) does the same except that   *
C *             it returns an integer*2 value IVALU2.           *
C *        CALL LIN(PROMPT,LVALUE) does the same except that    *
C *             it returns a logical value LVALUE.              *
C *        CALL TIN(PROMPT,TVALUE) does the same except that    *
C *             it returns a character value TVALUE.            *
C *        CALL NIN(PROMPT) merely prompts without reading      *
C *             any value.                                      *
C *                                                             *
C *  The character variables PROMPT and TVALUE can be of any    *
C *  length up to 80 for PROMPT and 132 for TVALUE.             *
C *  Writing is done on FORTRAN unit OUTPUT,                    *
C *  Reading from FORTRAN unit INPUT                            *
C *  (these default to "*" if set to zero )                     *
C *                                                             *
C ***************************************************************
 
      implicit none
 
      character*(*) prompt,tvalue
      character*80 code
      character*1 ans
      character*133 response
      logical lvalue
      real rvalue
      real*8 rvalu8
      integer max, mode, length, input, output, ivalue, iloc, lprompt
      integer*2 ivalu2
 
      include 'ctrlblk.inc'

      common / iodev / input, output
 
      entry rin(prompt, rvalue)
        mode=1
        max=80
        go to 1
 
      entry rin8(prompt, rvalu8)
        mode=12
        max=80
        go to 1
 
      entry iin(prompt, ivalue)
        mode=2
        max=80
        go to 1
 
      entry iin2(prompt, ivalu2)
        mode=22
        max=80
        go to 1
 
      entry lin(prompt, lvalue)
        mode=3
        max=80
        go to 1
 
      entry tin(prompt,tvalue)
        mode=4
        max = min0 (132, len (tvalue))
        go to 1
 
      entry nin(prompt)
        mode=5
        go to 1

      entry yin(prompt,ans)
        mode=6
        max = 80
        goto 1



    1 if( mode .eq. 1 ) then
        write (code, '(g12.5)', err=99) rvalue
        if (rvalue .eq. 0.) code = '0.0'
        call trunc(code, iloc)
 
      else if( mode .eq. 12 ) then
        write (code, '(g12.5)', err=99) rvalu8
        if (rvalu8 .eq. 0.d0) code = '0.d0'
        call trunc(code, iloc)
 
      else if( mode .eq. 2 ) then
        write (code, '(i12)', err=99) ivalue
        call trunc(code, iloc)
 
      else if( mode .eq. 22 ) then
        write (code, '(i12)', err=99) ivalu2
        call trunc(code, iloc)
 
      else if( mode .eq. 3 ) then
        write (code, '(l1)', err=99) lvalue
        call trunc(code, iloc)
 
      else if( mode .eq. 4 ) then
        call trunc(tvalue, iloc)
        if (iloc .eq. 0) then
          code = ' '
          iloc = 1
        else
          code = tvalue
        end if
 
      else if( mode .eq. 5 ) then
        if( iprint.ge.0 ) then
          lprompt = len(prompt)
          if( output .ne. 0 ) then
            write( output,603) prompt(:lprompt)
          else
            write(      *,603) prompt(:lprompt)
          end if
        end if
        return
 
      else if( mode.eq.6 ) then
        code = ans
        call shiftup( code )
        if( code.ne.'Y' .and. code.ne.'N' ) then
          write(*,*)'pinput: erroneous input to yin: ans = ', ans
          return
        endif
        iloc = 1

      end if
 
 
      if( iprint.ge.0 ) then
        lprompt = len(prompt)
        if( lprompt .gt. 0 ) then
          if( output .ne. 0 ) then
            write( output,600) prompt(:lprompt), code(:iloc)
          else
            write(      *,600) prompt(:lprompt), code(:iloc)
          end if
        else
          if( output .ne. 0 ) then
            write ( output, 601) code(:iloc)
          else
            write (      *, 601) code(:iloc)
          end if
        end if
      end if
 
      if( input.gt.0 ) then
        read( input,'(a)' ) response
      else
        read(     *,'(a)' ) response
      end if
      call trunc( response, iloc )
      length = iloc
 
10    format (a)
 
      if( length.eq.0 ) return
 
      if( mode.eq.4 .and. response.eq.' ') then
        tvalue = ' '
        return
      endif
 
      call trunc (response, iloc)
      if (iloc .eq. 0) return
 
c      if(mode.eq.1) read (response(:iloc),'(g12.0)',err=13) rvalue
c      if(mode.eq.12)read (response(:iloc),'(g12.0)',err=13) rvalu8
c      if(mode.eq.2) read (response(:iloc), '(bn,i12)', err=13) ivalue
c      if(mode.eq.22)read (response(:iloc), '(bn,i6)', err=13) ivalu2
c      if(mode.eq.3) read (response(:iloc), '(l1)', err=13) lvalue
      if(mode.eq.1) read(response(:iloc),*,err=13) rvalue
      if(mode.eq.12) read(response(:iloc),*,err=13) rvalu8
      if(mode.eq.2) read(response(:iloc),*,err=13) ivalue
      if(mode.eq.22) read(response(:iloc),*,err=13) ivalu2
      if(mode.eq.3) read(response(:iloc),*,err=13) lvalue

       if(mode.eq.4) tvalue = response
       if(mode.eq.6) then
        call shiftup(response)
        if( response.ne.'Y' .and. response.ne.'N' ) then
          goto 13
        endif
        ans = response
      endif


      return
 
   13 write(6,*) '***** Format error: try again *****'
      go to 1
 
99    write (6, *) 'error printing default value.  pinput aborted.'
 
  600 format( $, a,'  [default: ', a, ' ] >', $ )
  601 format ('+  [default: ', a, ' ]  >', $ )
  602 format( a, $ )
  603 format( a )
 
      END
