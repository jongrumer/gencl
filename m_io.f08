module m_io

   ! Input/Output module
   ! ===================

   ! List of routines in order of appearance
   !
   ! type         name           description
   ! -----------  -------------  -----------------------------------------------
   ! subroutine   s_help         prints help info if requested during runtime
   ! subroutine   s_input        process and checks input
   ! subroutine   s_print
   ! ---------------------------------------------------------------------------
   
   use m_parameters
   implicit none
   
contains

   subroutine s_help
      
      ! Input format help information

      implicit none
      
      ! local variables
      integer(kind=i4) :: i, j
      character(10)    :: str
 
 100  write (0,10001)
10001 format (//,5x,'This program prompts for each required input.  The user'/5x&
            & ,'should enter data or a RETURN after a question (?) mark')
      write (0,*)
      write (0,*) '     Example 1 :'
      write (0,*) '    --------------------'
      write (0,*) '                 Header  ?  S II ground state'
      write (0,*) '          Closed shells  ?   1s 2s 2p'
      write (0,*) '          Reference Set  ?  3s(2) 3p(3)'
      write (0,*) '                      2  ?  RETURN'
      write (0,*) '             Active Set  ?  3s,3p'
      write (0,*) 'Type of set generation   ?  0'
      write (0,*) '            Replacement  ?  3s(2) = 3d(2)'
      write (0,*) '                      2  ?  3s = 3d'
      write (0,*) '                      3  ?  3s.3p = 4s.3d'
      write (0,*) '                      4  ?  <RETURN>'
      write (0,*) '             Final Term  ?  4S'
      write (0,*) '                      2  ?  RETURN'
      write (0,10002)
10002 format (/5x,'Header and Closed Shells cannot exceed 72 characters and'/5x,&
             &'will be copied to the output file. The electrons are'/5x,        &
             &'separated by blanks in the Closed Shells.')
      write (0,10003)
10003 format (5x,'Press RETURN for more... ')
      read (5,'(A)') str
      do
 
         write (0,10004)
10004    format (///////5x,                                                     &
                &'Items are separated by a blank in the Reference Set, by a'/5x,&
                &'comma or a blank in the Active set, and by a period or a'/5x, &
                &'blank in Replacements.'//5x,                                  &
                &'Reference Set, Replacement, and Final Term are three sets'/5x,&
                &'of input, each with 0 to 10 members.  Each member must be'/5x,&
                &'entered on a separate line.')
 
         write (0,*) '       PRINT RETURN to terminate the input set.'
         write (0,*) '       PRINT RETURN if the set is empty.'
         write (0,*)
         write (0,*) '     Example 2 :'
         write (0,*) '    --------------------'
         write (0,*)
         write (0,*) '         Reference Set  ?  2s(1) 2p(2) 3s(1)'
         write (0,*) '                     2  ?  RETURN'
         write (0,*) '            Active Set  ?  2s,2p,3s'
         write (0,*) 'Type of set generation  ?  0'
         write (0,*) '           Replacement  ?  RETURN '
         write (0,*) '            Final Term  ?  RETURN'
         write (0,*)
         write (0,*) '     Where the Replacement and the Final Term are empty.'
         write (0,*)
         write (0,*)
         write (0,10010)
         read (5,'(A)') str
         i = index(str,'B')
         j = index(str,'b')
         if ( i/=0 .or. j/=0 ) go to 100
         do
 
            write (0,10005)
10005       format (///////5x,                                                  &
                  & 'By inputing "s" or "d" or "sd" you can compute the config-'&
                  & /5x,'urations from the Virtual Set, where  S means Single'/5&
                  & x,'Replacement, D means Double Replacement, SD means Single'&
                  & /5x,'and Double Replacement.'//5x,                          &
                  & 'GENCL will give you prompts for the Virtual Set automati-'/&
                  & 5x,                                                         &
                  & 'cally, then you need to specify the range of shells that'/5&
                  & x,'are to participate in the replacements. For instance, a'/&
                  & 5x,                                                         &
                  & 'response of 2 to "From which shell" and of 3 to "To which'/&
                  & 5x,'shell" implies that shells 2 and 3 participates in the'/&
                  & 5x,'replacements, and shell 1 does not enter into any',1x,  &
                  & 'replacements.')
            write (0,*)
            write (0,*) '     Example 3  :'
            write (0,*) '    -------------------------'
            write (0,*) '                ...'
            write (0,*) '         Reference Set  ?  2s(1) 2p(1) 3s(1)'
            write (0,*) '            Active Set  ?  RETURN'
            write (0,*) '           Replacement  ?  sd'
            write (0,*) '           Virtual Set  ?  3p,3d,4s'
            write (0,*) '      From which shell  ?  2'
            write (0,*) '        To which shell  ?  3'
            write (0,*) '            Final Term  ?  RETURN'
            write (0,*)
            write (0,*)
            write (0,*)
            write (0,10010)
            read (5,'(A)') str
            i = index(str,'B')
            j = index(str,'b')
            if ( i/=0 .or. j/=0 ) exit
            do
               write (0,10006)
10006          format (///////5x,                                               &
                 &'After terminating an input line, if you find the previous'   &
                 & /5x,                                                         &
                 &'input to be wrong, type ''B'' or ''b'' to go back one step.' &
                 & /5x,                                                         &
                 &'For example, before inputing Active Set, if you find the'/5x &
                 & ,'wrong spelling in the Header, type "B" and GENCL will '/5x,&
                 &'prompt for the Header again.'//)
               write (0,*) '     Example 4 :'
               write (0,*) '    ----------------------'
               write (0,*) '                Header  ?  OXYYEN'
               write (0,*) '            Active Set  ?  B'
               write (0,*) '                Header  ?  OXYGEN '
               write (0,*) '            Active Set  ?  2s '
               write (0,*)
               write (0,*) '     Then the following prompts will continue.'
               write (0,*)
               write (0,*)
               write (0,10010)
               read (5,'(A)') str
               i = index(str,'B')
               j = index(str,'b')
               if ( i/=0 .or. j/=0 ) exit
               do
 
                  write (0,10007)
10007             format (///////5x,                                            &
                     & 'Example 4 shows the procedure for going back four steps'&
                     & /5x,'to correct the Closed Shells from  5s  to 1s  2s.'//&
                     & )
                  write (0,*) '     Example 5 :'
                  write (0,*) '    -----------------------'
                  write (0,*)
                  write (0,*) '            Active Set  ?  5s'
                  write (0,*) 'Type of set generation  ?  0'
                  write (0,*) '         Reference Set  ?  2s(1) 2p(2) 3s(1)'
                  write (0,*) '                     2  ?  2P(4)'
                  write (0,*) '                     3  ?  b '
                  write (0,*) '                     2  ?  b '
                  write (0,*) '         Reference Set  ?  b '
                  write (0,*) '         Closed Shells  ?    1s  2s '
                  write (0,*) '         Reference Set  ?  '
                  write (0,*)
                  write (0,*)
!
                  write (0,*) '   Then reenter the data for the Reference Set', &
                             &' and continue the input.'
                  write (0,*)
                  write (0,10010)
                  read (5,'(A)') str
                  i = index(str,'B')
                  j = index(str,'b')
                  if ( i/=0 .or. j/=0 ) exit
 
                  write (0,10008)
10008             format (///////5x,                                            &
                   &'When the following error conditions are detected,',        &
                   &' a message is given.'/10x,                                 &
                   &' 1). The parentheses are not matched'/10x,                 &
                   &' 2). The number of electrons in a shell is more than FULL' &
                   & //15x,'For each member of the Reference Set ,'/10x,        &
                   &' 3). The number of electrons is not the same '/10x,        &
                   &' 4). The parity is not the same '//10x,                    &
                   &' 5). The number of couplings generated by a configuration' &
                   & /10x,'     is more than 500'///////)
                  write (0,10009)
10009             format (5x,'Press ''b'' for BACK or RETURN to begin the',1x,  &
                         &'program.')
                  read (5,'(A)') str
                  i = index(str,'B')
                  j = index(str,'b')
                  if ( i==0 .and. j==0 ) then
                     do i = 1, 30
                        write (0,*)
                     end do
                     go to 99999
                  end if
               end do
            end do
         end do
      end do
10010 format (5x,'Press ''b'' for BACK or RETURN for more... ')
99999 continue

   end subroutine s_help
   
   subroutine s_input(maxset,nset,set,mark,sflag,dflag,*,*)
      
      ! Process the input set and check the input error
      
      use m_tools, only: s_del

      implicit none
      
      ! dummy arguments
      integer(i4),                 intent(in)    :: maxset, mark
      integer(i4),                 intent(inout) :: nset, sflag, dflag
      character(60), dimension(*), intent(inout) :: set
      ! local variables
      character(60) :: ch1, temp
      integer(i4)   :: i, j

      !     maxset   =  maximum number of input elements
      !     nset     =  number of members in the set
      !     set      =  character array with nset elements
      !     *1       =  return label if input is 'B'
      !     *2       =  return label if the set is empty
      !     mark     =  1 if input is replacement, 0 otherwise

      nset = 0
      set(1) = '    '
      do while ( nset<maxset )
         read (5,'(A60)') temp
 
         ! ... If input is 'B' or 'b', go back one step
         i = index(temp,'B')
         j = index(temp,'b')
         if ( i/=0 .or. j/=0 ) then
            if ( nset==0 ) then
               return 1
            else
               nset = nset - 1
               go to 50
            end if
         end if
 
         ! ... Go to for next input if the input is empty
         !     Return if the input is finished
         if ( temp(1:5)=='     ' ) then
            if ( nset==0 ) then
               return 2
            else
               return
            end if
         end if
         call s_del(temp)
 
         ! ... If Replacement is 's' or 'd' or 'sd', set single and
         !     double flag for Virtual Set
         if ( mark/=0 ) then
            ch1 = temp(:1)
            if ( ch1=='S' .or. ch1=='s' ) sflag = 1
            i = index(temp,'SD')
            j = index(temp,'sd')
            if ( ch1=='D' .or. ch1=='d' .or. i/=0 .or. j/=0 ) dflag = 1
            if ( sflag/=0 .or. dflag/=0 ) then
               nset = 1
               set(1) = temp
               return
            end if
         end if
 
         ! ... READ the input and delete the repeated member
         do i = 1, nset
            if ( set(i)==temp ) then
               write (0,*) '     You give a repeated input!'
               go to 50
            end if
         end do
         nset = nset + 1
         set(nset) = temp
 50      write (0,10001,advance='no') nset + 1, ': '
10001    format (t24,i10,a2)
      end do

   end subroutine s_input

   subroutine s_print(e_lbl,q_numb,n_shells,nc,mark)
      
      ! Print out the values of couplings
      ! =================================

      !   input :
      !        e_lbl =  electron labels
      !                 where e_lbl(i)(1=1)  ---  blank
      !                       e_lbl(i)(2=2)  ---  n-symbol
      !                       e_lbl(i)(3=3)  ---  L-symbol
      !       q_numb =  occupation number
      !                       0 (empty) <= q_numb(i)  <= 2(2L(i)+1) (full)
      !     n_shells =  number of shells
      !                       0  <=  n_shells  <=  5
      !           nc =  number of couplings
      !


      use m_parameters
      use m_globals
      implicit none

      ! dummy arguments
      integer(kind=i4), intent(in)                  :: n_shells, mark, nc
      integer(kind=i4), dimension(nels), intent(in) :: q_numb
      character(3),     dimension(nels), intent(in) :: e_lbl
      ! local variables
      integer   :: i, j, k
      character :: ch1
      character(3), dimension(ncoupl) :: couple
 
      ! mark  =  1 for reference set
      !       =  2 for reference set, active set
      !       =  3 for reference set, replacement
      !       =  4 for reference set, virtual set

      ! ... print the input of header
      j = index(header,'          ')
      write (6,10001) header(1:j)
10001 format (' '/t10,'-------------       ',a,'  --------',/)
 
      ! ... print the input of Closed Shells
      j = index(shells,'     ')
      write (6,10002) shells(:j)
10002 format ('          Closed Shells  :  ',a/)
      if ( mark/=1 .and. mark/=2 ) then
 
         ! ... print the input of Reference Set
         write (6,10003) ref(nf)
10003    format ('          Reference set  :  ',a60/)
         if ( mark/=3 ) then
 
            ! ... print the input of Virtual Set
            if ( mark==4 ) ch1 = 'S'
            if ( mark==5 ) ch1 = 'D'
            if ( mark==4 .or. mark==5 ) then
               write (6,10004) virtul
10004          format ('            Virtual Set  :  ',a60/)
               write (6,10005) ch1, repl(nels)
10005          format ('          ',a1,'-Replacement  :  ',a60/)
               go to 100
            end if
 
            ! ... print the input of Active Set
            write (6,10006) act
10006       format ('             Active Set  :  ',a60/)
            if ( mark==2 ) go to 100
         end if
 
         ! ... print the input of Replacement
         write (6,10007) repl(nr)
10007    format ('            Replacement  :  ',a60/)
      end if
 
      ! ... print the new configuration by Replacement
 100  if ( e_lbl(1)(1:1)==' ' ) then
         k = 1
      else
         k = 2
      end if
      write (6,10008) (e_lbl(j),q_numb(j),j=1,n_shells)
10008 format ('          Configuration  :',5(2x,a3,'(',i2,')')/)
 
      ! ... print couplings generated from the configuration
      k = 2*n_shells - 1
      if ( n_shells<=3 ) then
         read (file3(1),10013) (couple(j),j=1,k)
         write (6,10009) (couple(j),j=1,k)
10009    format ('        Their couplings  :  ',9(5x,a3))
         do i = 2, nc
            read (file3(i),10013) (couple(j),j=1,k)
            write (6,10010) (couple(j),j=1,k)
10010       format (t29,9(5x,a3))
         end do
      else
         write (6,10011)
10011    format ('         Their couplings  :'/)
         do i = 1, nc
            read (file3(i),10013) (couple(j),j=1,k)
            write (6,10012) (couple(j),j=1,k)
10012       format (t5,9(5x,a3))
         end do
      end if

10013 format (9(a3))

   end subroutine s_print

end module m_io
