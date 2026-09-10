!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_String_Functions </br></br>
!#
!# **Brief**: Routines to manipulate Strings - http://fortranwiki.org/fortran/show/String_Functions</br></br>
!# 
!# **Author**: David Frank </br>
!#
!# **Version**: 2.0.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>??-??-???? - David Frank - version: 1.0.0 </li>
!#  <li>07-08-2019 - Denis Eiras - version: 2.0.0 </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time</li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> </br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> </br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---

module Mod_String_Functions  ! by David Frank  dave_frank@hotmail.com
  implicit none            ! http://home.earthlink.net/~dave_gemini/strings.f90


  interface Copy
    !# Copies (generic) char array to string or string to char array
    !# ---
    !# @info
    !# **Brief:** Copies (generic) char array to string or string to char array. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    module procedure copy_a2s, copy_s2a
  end interface Copy


contains


  pure function Copy_a2s(a)  result (s)
    !# Copies char array to string
    !# ---
    !# @info
    !# **Brief:** Copies char array to string. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character, intent(in) :: a(:)
    character(size(a)) :: s
    integer :: i
    do i = 1, size(a)
      s(i:i) = a(i)
    end do
  end function Copy_a2s


  pure function Copy_s2a(s)  result (a)
    !# Copies s(1:Clen(s)) to char array
    !# ---
    !# @info
    !# **Brief:** Copies s(1:Clen(s)) to char array. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*), intent(in) :: s
    character :: a(len(s))
    integer :: i
    do i = 1, len(s)
      a(i) = s(i:i)
    end do
  end function Copy_s2a


  pure integer function Clen(s)
    !# Returns same result as len unless last non-blank char = null
    !# ---
    !# @info
    !# **Brief:** Returns same result as len unless last non-blank char = null. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*), intent(in) :: s       ! last non-blank char is null
    integer :: i
    Clen = len(s)
    i = len_trim(s)
    if (s(i:i) == char(0)) Clen = i - 1  ! len of C string
  end function Clen


  pure integer function Clen_trim(s) 
    !# Returns same result as len_trim unless last non-blank char = null
    !# ---
    !# @info
    !# **Brief:** Returns same result as len_trim unless last non-blank char = null. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*), intent(in) :: s       ! last char non-blank is null, if true:
    integer :: i                       ! then len of C string is returned, note:
    ! Ctrim is only user of this function
    i = len_trim(s) ; Clen_trim = i
    if (s(i:i) == char(0)) Clen_trim = Clen(s)   ! len of C string
  end function Clen_trim


  function Ctrim(s1)  result(s2)
    !# Returns same result as trim unless last non-blank char = null
    !# ---
    !# @info
    !# **Brief:** Returns same result as trim unless last non-blank char = null. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*), intent(in) :: s1     ! last non-blank char is null in which
    character(Clen_trim(s1)) :: s2     ! case trailing blanks prior to null
    s2 = s1                            ! are output
  end function Ctrim


  integer function Count_Items(s1)
    !# Counts items in string or C string that are blank or comma separated
    !# ---
    !# @info
    !# **Brief:** Counts items in string or C string that are blank or comma separated. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s1
    character(Clen(s1)) :: s
    integer :: i, k

    s = s1                            ! remove possible last char null
    k = 0  ; if (s /= ' ') k = 1      ! string has at least 1 item
    do i = 1, len_trim(s) - 1
      if (s(i:i) /= ' '.and.s(i:i) /= ',' &
        .and.s(i + 1:i + 1) == ' '.or.s(i + 1:i + 1) == ',') k = k + 1
    end do
    Count_Items = k
  end function Count_Items


  function Reduce_Blanks(s)  result (outs)
    !# Reduces blanks in string to 1 blank between items, last char not blank
    !# ---
    !# @info
    !# **Brief:** Reduces blanks in string to 1 blank between items, last char not blank. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s
    character(len_trim(s)) :: outs
    integer :: i, k, n

    n = 0  ; k = len_trim(s)          ! k=index last non-blank (may be null)
    do i = 1, k - 1                      ! dont process last char yet
      n = n + 1 ; outs(n:n) = s(i:i)
      if (s(i:i + 1) == '  ') n = n - 1  ! backup/discard consecutive output blank
    end do
    n = n + 1  ; outs(n:n) = s(k:k)    ! last non-blank char output (may be null)
    if (n < k) outs(n + 1:) = ' '       ! pad trailing blanks
  end function Reduce_Blanks


  function Replace_Text (s, text, rep)  result(outs)
    !# Replaces text in all occurances in string with replacement string
    !# ---
    !# @info
    !# **Brief:** Replaces text in all occurances in string with replacement string. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s, text, rep
    character(len(s) + 100) :: outs     ! provide outs with extra 100 char len
    integer :: i, nt, nr

    outs = s ; nt = len_trim(text) ; nr = len_trim(rep)
    do
      i = index(outs, text(:nt)) ; if (i == 0) EXIT
      outs = outs(:i - 1) // rep(:nr) // outs(i + nt:)
    end do
  end function Replace_Text


  function Spack (s, ex)  result (outs)
    !# Spack pack string's chars == extract string's chars
    !# ---
    !# @info
    !# **Brief:** Spack pack string's chars == extract string's chars. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s, ex
    character(len(s)) :: outs
    character :: aex(len(ex))   ! array of ex chars to extract
    integer :: i, n

    n = 0  ;  aex = Copy(ex)
    do i = 1, len(s)
      if (.not.any(s(i:i) == aex)) cycle   ! dont pack char
      n = n + 1 ; outs(n:n) = s(i:i)
    end do
    outs(n + 1:) = ' '     ! pad with trailing blanks
  end function Spack


  integer function Tally (s, text)
    !# Talles occurances in string of text arg
    !# ---
    !# @info
    !# **Brief:** Talles occurances in string of text arg. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s, text
    integer :: i, nt

    Tally = 0 ; nt = len_trim(text)
    do i = 1, len(s) - nt + 1
      if (s(i:i + nt - 1) == text(:nt)) Tally = Tally + 1
    end do
  end function Tally


  function Translate(s1, codes)  result (s2)
    !# Translates text arg via indexed code table
    !# ---
    !# @info
    !# **Brief:** Translates text arg via indexed code table. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s1, codes(2)
    character(len(s1)) :: s2
    character :: ch
    integer :: i, j

    do i = 1, len(s1)
      ch = s1(i:i)
      j = index(codes(1), ch) ; if (j > 0) ch = codes(2)(j:j)
      s2(i:i) = ch
    end do
  end function Translate


  function Upper(s1)  result (s2)
    !# Returns in upper case the text arg
    !# ---
    !# @info
    !# **Brief:** Returns in upper case the text arg. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s1
    character(len(s1)) :: s2
    character :: ch
    integer, parameter :: DUC = ichar('A') - ichar('a')
    integer :: i

    do i = 1, len(s1)
      ch = s1(i:i)
      if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch) + DUC)
      s2(i:i) = ch
    end do
  end function Upper


  function Lower(s1)  result (s2)
    !# Returns in lower case the text arg
    !# ---
    !# @info
    !# **Brief:** Returns in lower case the text arg. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    character(*) :: s1
    character(len(s1)) :: s2
    character :: ch
    integer, parameter :: DUC = ichar('A') - ichar('a')
    integer :: i

    do i = 1, len(s1)
      ch = s1(i:i)
      if (ch >= 'A'.and.ch <= 'Z') ch = char(ichar(ch) - DUC)
      s2(i:i) = ch
    end do
  end function Lower


  function intToStr(intParam) result (intString)
    !# Converts integer to string
    !# ---
    !# @info
    !# **Brief:** Converts integer to string. </br>
    !# **Authors**: </br>
    !# &bull; David Frank </br>
    !# **Date**: aug/2019 </br>
    !# @endinfo
    implicit none
    integer :: intParam
    character(4) :: intString
    character(4) :: formatString

    if(intParam < 10) then
      formatString = "(I1)"
    elseif (intParam < 100) then
      formatString = "(I2)"
    elseif (intParam < 1000) then
      formatString = "(I3)"
    else
      formatString = "(I4)"
    endif
    write(intString, fmt = formatString) intParam
  end function intToStr

end module Mod_String_Functions
