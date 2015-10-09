module m_tools
   
   ! Tools module
   ! ===================

   ! List of routines in order of appearance
   !
   ! type         name           description
   ! -----------  -------------  -----------------------------------------------
   ! function     f_ictoi        converts a string into its corresp. integer
   ! function     f_lval         converts a symbol into its corresp. quant. num.
   ! function     f_symb         converts a quant. num. into its corresp. symbol
   ! subroutine   s_del          delete the leading space of a given string
   ! subroutine   s_strsh        shifts a string left to a given starting point
   ! subroutine   s_decomp       decomposes a string of electron replacements
   ! ---------------------------------------------------------------------------
   
   use m_parameters
   
   implicit none

contains
   
   function f_ictoi(str)

      ! Convert character string into its corresponding integer

      implicit none
      
      ! function declaration
      integer(i4) :: f_ictoi
      ! dummy arguments
      character(*), intent(in) :: str
      ! local variables
      integer(i4) :: n
      
      n = ichar(str(1:1)) - ichar('0')
      
      if(len(str) > 1) then
         if( str(2:2) /= ' ' ) then
            n = n*10 + ichar(str(2:2)) - ichar('0')
         end if
      else
         write(*,*) 'Error in f_ictoi: string has length < 2. Program stop!'
         stop
      end if
      
      f_ictoi = n

   end function f_ictoi

   function f_lval(symbol)

      ! convert the symbol into its corresponding quantum number

      implicit none
      
      ! function declaration
      integer(i4) :: f_lval
      ! dummy arguments
      character,intent(in) :: symbol
      ! local variables
      integer(i4)          :: locate
      character(40), save  :: set = 'spdfghiklmnopqrstuvwSPDFGHIKLMNOPQRSTUVW'

      locate = index(set,symbol)
      if( locate <= 20_i4 ) then
         f_lval = locate - 1_i4
      else
         f_lval = locate - 21_i4
      end if

   end function f_lval

   function f_symb(l)

      !  quantum number into its corresponding symbol
      
      implicit none
      
      ! function declaration
      character :: f_symb
      ! dummy arguments
      integer(i4), intent(in) :: l
      ! local variables
      character(len=20), save :: set = 'SPDFGHIKLMNOPQRSTUVW'
 
      f_symb = set(l+1:l+1)

   end function f_symb

   subroutine s_del(str)

      ! Delete the leading space of the string

      implicit none
     
      ! dummy arguments
      character(*), intent(inout) :: str
      
      ! local variables
      integer(i4)   :: i, length
      character(72) :: temp
 
      length = len(str)
      i      = 0_i4

      do while ( str(i+1_i4:i+1_i4)==' ' )
         i = i + 1_i4
         if ( i>=length ) exit
      end do

      temp = str(i+1_i4:)
      str  = temp

   end subroutine s_del

   subroutine s_strsh(str,i)
      
      ! Shift the string left
      
      implicit none

      ! dummy arguments
      integer(i4),  intent(in)    :: i
      character(*), intent(inout) :: str
      ! local variables
      character(72) :: temp
      
      temp = str(i:)
      str  = temp

   end subroutine s_strsh

   subroutine s_decomp(string,e_lbl,q_numb,n_repl)
      
      ! Decompose the string of Replacements
      ! ====================================
      
      use m_globals

      implicit none

      ! input  :
      !     string = string to be decomposed
      !     e_lbl  = electron label array
      !                 where e_lbl(i)(1=1)  ---  blank
      !                       e_lbl(i)(2=2)  ---  n-symbol
      !                       e_lbl(i)(3=3)  ---  L-symbol
      !     q_numb =  occupation number
      !                       0 (empty) <= q_numb(i)  <= 2(2L(i)+1) (full)
      ! output  :
      !     n_repl =  number of e_lbl to be replaced

      ! dummy arguments
      character(60),                 intent(inout) :: string
      character(3), dimension(nels), intent(out)   :: e_lbl
      integer(i4),  dimension(nels), intent(out)   :: q_numb
      integer(i4),                   intent(out)   :: n_repl
      
      ! local variables
      character(3) :: ch3
      integer :: i, j, k, left, n, right
      
      do i = 1, 5
         if ( string(:5)=='     ' ) then
            n_repl = i - 1
            return
         end if
         call s_del(string)
 
         left  = index(string,'(')
         right = index(string,')')
 
         ! ... If the Replacement is like 2s.2p = 3s.3p
         if ( left==0 ) then
            k = 1
            n = ichar(string(3:3))
            if ( n>=ord0 .and. n<=ord9 ) then
               j = 3
            else
               j = 2
            end if
            ch3 = string(:j)
 
         ! ... If the Replacement is like 2p(2) = 3p(2)
         else
            ch3 = string(left+1:right-1)
            k = f_ictoi(ch3)
            ch3 = string(:left-1)
         end if
 
         ! ... Convert uppercase to lowercase, and assign value to e_lbl(i),q_numb(i)
         n = ichar(ch3(2:2))
         if ( n>=ordua .and. n<=orduz ) ch3(2:2) = char(n-ordua+ordla)
         e_lbl(i) = ch3
         q_numb(i) = k
         if ( left==0 ) then
            call s_strsh(string,(j+2))
         else
            call s_strsh(string,(right+1))
         end if
      end do

   end subroutine s_decomp

end module m_tools
