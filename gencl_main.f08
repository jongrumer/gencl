program gencl_main
      
      use m_parameters
      use m_globals
      use m_tools
      use m_io
      use m_generate

      implicit none

      ! all function nams start with f_ and subroutines with s_

      ! local variables
      
      ! left      =  position of the first '(' in the string
      ! right     =  position of the first ')' in the string
      ! sumq      =  total number of q(i) for each configuration
      ! nl        =  total quantum number for each reference set
      ! const     =  total number of qi for the first reference set
      ! elv       =  eli for virtual set
      ! pl        =  parity of qi for reference set
      ! llmin     =  minimum angular coupling of the pair
      ! llmax     =  maximum angular coupling of the pair
      ! str       =  string to be packed as output for replacement
      ! elv       =  electron label for virtual set
      ! mv        =  number of electron for virtual set
      ! vl        =  parity for each qi in virtual set

      ! header    =  csf list header
      ! shells    =  csf list closed shell(s)
      ! act       =  csf list active set of orbital(s)
      ! ref       =  csf list reference configuration(s)
      ! repl      =  csf list number of replacements
      ! finals    =  csf list final term(s)
      ! nref      =  number of reference sets
      ! nrpl      =  number of replacements
      ! nftm      =  number of final terms
      ! nels      =  number of electrons

      ! PROCESS THE REFERENCE SET
      !
      ! m         =  maximin number of eli for reference set
      ! els       =  storage of eli for different reference set
      ! qs        =  storage of qi for each reference set
      ! ms        =  storage of m for each reference set
      ! strl      =  temporary string
      
      ! PROCESS THE ACTIVE SET - the member of the input is limited to nels
      !
      ! nq        =  number of configurations
      ! ela       =  eli for active set
      ! qa        =  qi for active set
      ! ma        =  number of eli for active set
      ! rl        =  L-value for each shell
      ! f(i)      =  maximum allowed number of electrons for each shell
      ! h(i)      =  number of unassigned electrons or ``holes'' for each shell
      
      ! rl        =  parity of eli of reference set
      ! vl        =  parity of eli of virtual set
      ! pl        =  pairty for the pair qi and qj in reference set

      real(kind=r8)     :: s, sumq
      integer(kind=i4)  :: i, j, k, l
      integer(kind=i4)  :: dflag, iflag, sflag, itype, left, right
      integer(kind=i4)  :: ll1, ll2, llmax, llmin, lmax, lmin, lvir, nvir
      integer(kind=i4)  :: minq1, mv, n, n1, n2, nc, nl, nrepl, pl
      character(len=1)  :: ch1
      character(len=2)  :: ch2
      character(len=3)  :: ch3
      character(len=72) :: strl, strr
      character(len=72) :: temp
      character(len=3),  dimension(ncoupl) :: couple
      character(len=60), dimension(nels)   :: finals
      character(len=3),  dimension(nels)   :: elb, elc, elv
      integer(kind=i4),  dimension(nels)   :: f, h, vl

      integer(kind=i4)  :: q1, q2, q3, q4, q5      ! occupation numbers (1-5)
      integer(kind=i4)  :: q6, q7, q8, q9, q10     ! occupation numbers (6-10)
      integer(kind=i4)  :: q11, q12, q13, q14, q15 ! occupation numbers (11-15)

      nq    = 0 ! init number of configurations
      iflag = 0 ! and some flags
      sflag = 0 !
      dflag = 0 !

      ! ... Print program header

      write (0,*)
      write (0,*) ' -----------------------------------------------------------'
      write (0,*) '                                                            '
      write (0,*) '                        G E N C L 9 0                       '
      write (0,*) '                        =============                       '
      write (0,*) '            A code to GENerate Configuration Lists          '
      write (0,*) '                                                            '
      write (0,*) '      Modern Fortran version with support for old (ATSP)    '
      write (0,*) '          and new (ATSP2K) style coupling formatting        '
      write (0,*) '                                                            '
      write (0,*) ' -----------------------------------------------------------'
      write (0,*)

      ! ... if the user needs to obtain information about input format,
      !     call subroutine s_help
      
      write (0,'(a59)', advance='no') ' - Type "h" if you need Help about the input format [h/*]: '
      read (5,'(a)') strl
      i = index(strl,'H')
      j = index(strl,'h')
      if ( i/=0 .or. j/=0 ) call s_help

      ! ... read and analysis the input data

      ! ... get old or new style coupling formatting 
      !     n => pre atsp2k style, y => atsp2k style (denser)
 50   write (0,'(a52)', advance='no') ' - New style coupling formatting (ATSP2K)? [*/n/b]: '
      read (5,'(a1)') answ_coupl_fmt_new

      i = index(answ_coupl_fmt_new,'B')
      j = index(answ_coupl_fmt_new,'b')
      if ( i/=0 .or. j/=0 ) go to 50
      if ( answ_coupl_fmt_new == 'N') answ_coupl_fmt_new = 'n'
      if ( answ_coupl_fmt_new == 'Y') answ_coupl_fmt_new = 'y'

      ! ... get header string
 100  write (0,'(a21)',advance='no') ' - Give list header: '
      read (5,'(a)') header
 
      !ch1 = ch3(2:2) ! CHECK THIS! IS CH3 REALLY DEFINED?
      
 200  do
         ! ... Closed Shells - check if it satisfies the FORMAT(18(1X,A3))
         write (0,'(a23)',advance='no') ' - Give closed shells: '
         read (5,'(a)') temp

         i = index(temp,'B')
         j = index(temp,'b')
         if ( i/=0 .or. j/=0 ) go to 100

         n      = 2
         shells = ' '
         do
            if ( index(temp,'     ')>2 ) then
               call s_del(temp)
               j = ichar(temp(1:1))
               if ( j<ord0 .or. j>ord9 ) then
                  write (0,*) '     INPUT ERROR !'
                  exit
               end if
               shells(n:n+2) = temp(:3)
               n = n + 4
               call s_strsh(temp,4)
               cycle
            end if
            go to 300
         end do
      end do
 300  do
 
         ! ... Input Reference Set, and check the input error

         write (0,'(a35)',advance='no') ' - Give reference configuration 1: '
         
         call s_input(nels,nref,ref,iflag,sflag,dflag,*200,*300)

         do i = 1, nref
            sumq = 0
            nl   = 0
            temp = ref(i)
            
            do
               if ( temp(1:5)/='     ' ) then
                  
                  call s_del(temp)
                  
                  ! ... error if the input has unmatched parenthesis
                  left  = index(temp,'(')
                  right = index(temp(:7),')')
                  
                  if ( left==0 .or. right==0 ) then
                     write (0,*) '       Unmatched parenthesis!',               &
                                &' Please input again .'
                     go to 400
                  end if
                  
                  ! ... error if number of electron is more than FULL
                  ch3 = temp(left+1:right-1)
                  n   = f_ictoi(ch3)
                  l   = f_lval(temp(2:2))
                  
                  if ( n>l*4+2 ) then
                     write (0,*)                                                &
                     &'     Number of electrons in a shell is  more than FULL !'
                     go to 400
                  end if
                  sumq = sumq + n
                  nl   = nl + n*l
                  call s_strsh(temp,right+1)
                  cycle
               end if
 
               ! ... The first Reference Set is defined as the standard value
               if ( i==1 ) then
                  const     = int(sumq)
                  parityval = mod(nl,2)
               else

                  ! ... Error if members of the reference set have different
                  !     electron numbers or parity
                  if ( int(sumq)/=const ) then
                     write (0,*) '     Total number of electrons is wrong!'
                     go to 400
                  end if
                  if ( mod(nl,2)/=parityval ) then
                     write (0,*) '       Parity is wrong!'
                     go to 400
                  end if
               end if
               exit
            end do
         end do

         exit
 400  end do
 
      ! ... Input Active Set
 500  write (0,'(a20)',advance='no') ' - Give active set: '
      read (5,'(A)') act
      
      i = index(act,'B')
      j = index(act,'b')

      if ( i/=0 .or. j/=0 ) go to 300
      if ( act/='   ' ) then
         write (0,'(a29)',advance='no') ' - Give set generation type: '
         read *, itype
      end if

 600  do
         ! ... Input Replacement
         write (0,'(a22)',advance='no') ' - Give replacements: '
         iflag = 1

         call s_input(nels,nrepl,repl,iflag,sflag,dflag,*500,*800)
 
         ! ... if replacement = S or D or SD, input virtual set
         if ( sflag==0 .and. dflag==0 ) exit
 
         write (0,'(a21)',advance='no') ' - Give virtual set: '
         read (5,'(a)') virtul

         temp      = virtul
         j         = index(temp,'     ')
         temp(j:j) = ','
         mv        = 0
 
         ! ... Decompose the input of Virtual Set
         do while ( temp(:5)/='     ' )
            mv = mv + 1
            if ( mv>(nels) ) then
               write (0,*) ' Virtual set too large: max = ', nels
               stop
            end if
            
            call s_del(temp)

            j      = index(temp,',')
            vl(mv) = mod(f_lval(temp(2:2)),2)
            
            ! ... Convert the input of uppercase to lowercase and assign value ELVi
            n      = ichar(temp(2:2))
            if ( n>=ordua .and. n<=orduz ) then
               temp(2:2) = char(n-ordua+ordla)
            end if
            
            elv(mv) = temp(:j-1)
            
            call s_strsh(temp,j+1)

         end do
         do
 
            write (0,'(a15)',advance='no') ' - From shell: '
            read (5,'(A)') temp

            i = index(temp,'B')
            j = index(temp,'b')
            
            if ( i/=0 .or. j/=0 ) exit
            
            call s_del(temp)
            nvir = ichar(temp(1:1)) - ord0
            
            if ( nvir<1 .or. nvir>9 ) then
               write (0,*) '    Please input digit for the shell position!'
               cycle
            end if
            
            do
               write (0,'(a13)',advance='no') ' - To shell: '
               read (5,'(a)') temp
               
               i = index(temp,'B')
               j = index(temp,'b')
               
               if ( i/=0 .or. j/=0 ) go to 700
               
               call s_del(temp)
               lvir = ichar(temp(:1)) - ord0
               
               if ( lvir<0 .or. lvir>9 ) then
                  write (0,*) '     Please input digit for the shell position!'
                  cycle
               end if
               
               go to 800
            end do
         end do
 700  end do
 
      ! ... Input Final Terms

 800  write (0,'(a21)',advance='no') ' - Give final terms: '
      
      iflag = 0
      call s_input(nels,nftm,finals,iflag,sflag,dflag,*600,*900)
      
      do i = 1, nftm
         ch2 = finals(i)(:2)
         
         ! ... Conver the input of lowercase to uppercase
         n = ichar(ch2(2:2))
         if ( n>=ordla .and. n<=ordlz ) then
            ch2(2:2) = char(n-ordla+ordua)
         end if
         ftm(i) = ch2
      end do
 
      ! ... open the file for ci.lst
 900  open (7,file=fname,status='UNKNOWN')

      ! ... print out all input data of the user

      write (6,10011)
10011 format (///////t5,'***************          I N P U T',                   &
             &'  D A T A          **********'///)
      write(6,'(a39,1x,a1)') ' New style (ATSP2K) coupling format :  ', answ_coupl_fmt_new
      write (6,10012) header
10012 format (t5,'          Header  :  ',a60/)
      write (6,10013) shells
10013 format (t5,'   Closed shells  :  ',a60/)
      write (6,10014) ref(1)
10014 format (t5,'   Reference Set  :  ',a60/)
      do i = 2, nref
         write (6,10033) i, ref(i)
      end do
      write (6,10015) act
10015 format (t5,'      Active Set  :  ',a60/)
      write (6,10016) repl(1)
10016 format (t5,'    Replacements  :  ',a60/)
      if ( sflag/=0 .or. dflag/=0 ) then
         write (6,10017) virtul
10017    format (t5, '     Virtual Set  :  ',a72/)
         write (6,*) '             From which shell  :  ', char(nvir+ord0)
         write (6,*)
         write (6,*) '               To which shell  :  ', char(lvir+ord0)
         write (6,*)
      end if
      do i = 2, nrepl
         write (6,10033) i, repl(i)
      end do
      write (6,10018) finals(1)
10018 format (t5,'     Final Terms  :  ',a60/)
      do i = 2, nftm
         write (6,10033) i, finals(i)
      end do
 
      write (7,10019) header
10019 format (' ',a60)
      write (7,'(1X,A72)') shells
 
      ! ... process the reference set

      ! m     =  maximin number of eli for reference set
      ! els   =  storage of eli for different reference set
      ! qs    =  storage of qi for each reference set
      ! ms    =  storage of m for each reference set
      ! strl  =  temporary string

      do nf = 1, nref
         m = 0
         temp = ref(nf)
         do i = 1, nels
 
            ! ...  Decompose the input for Reference Set
            if ( temp(:5)=='     ' ) then
               q(i) = 0
               exit
            end if

            call s_del(temp)
            
            left  = index(temp,'(')
            right = index(temp,')')
            m     = m + 1
 
            ! ... convert the input of uppercase to lowercase, and assign initial
            !     values for ELi and Qi
            ch3   = temp(:left-1)
            n     = ichar(ch3(2:2))
            
            if ( n>=ordua .and. n<=orduz ) then
               ch3(2:2) = char(n-ordua+ordla)
            end if
            
            el(m) = ch3
            ch3   = temp(left+1:right-1)
            q(m)  = f_ictoi(ch3)

            call s_strsh(temp,right+1)
         end do

         do j = 1, nels
            qs(j,nf) = 0
         end do
 
         ! ... Add new electrons to the first configuration
         if ( nf==1 ) then
            do j = 1, m
               els(j)   = el(j)
               qs(j,nf) = q(j)
            end do
            ms(nf) = m
         else
            do j = 1, m
               n = ms(1)
               do i = 1, n
                  if ( el(j)==els(i) ) then
                     qs(i,nf) = q(j)
                     go to 920
                  end if
               end do
               n = n + 1
               if ( n>nels ) then
                  write (0,*) ' Too many shells in reference set: max = ', nels
                  stop
               end if
               ms(1)    = n
               els(n)   = el(j)
               qs(n,nf) = q(j)
 920        end do
            ms(nf) = n
         end if
 
         ! ... generate all couplings for the reference set
         !     max = -5  means the first time to call subroutine s_coupld
         k = 0
         do i = 1, m
            if ( q(i)/=0 ) then
               k = k + 1
               elb(k) = el(i)
               qb(k)  = q(i)
            end if
         end do
 
         maxfinal = -5
         call s_coupld(elb,qb,k,elc,nc,*1100)
         if ( nc/=0 ) then
            write (6,10020)
10020       format (//                                                          &
           &'      GENERATE ALL COUPLINGS FOR EACH MEMBER  OF THE REFERENCE SET'&
          & //)
            call s_print(elc,qb,k,nc,1)
         end if
 
         ! ... Compute max and min value for the finals.
         !     Rule :  MAXFINAL = 2*|S+L| ,  MINFINAL = 2*|S-L|
         if ( nf==1 ) then
            maxfinal = 0
            minfinal = 100
         
            do i = 1, nc
               n = 2*k - 1
               read (file3(i),'(9(a3))') (couple(j),j=1,n)
               
               ch3  = couple(n)(1:1)
               ch1  = couple(n)(2:2)
               s    = real((f_ictoi(ch3)-1),r8)/2.0_r8
               l    = f_lval(ch1)
               lmin = 2_i4*abs( int(s,i4) - l )
               lmax = 2_i4*abs( int(s,i4) + l )
               
               if ( lmin<minfinal ) minfinal = lmin
               if ( lmax>maxfinal ) maxfinal = lmax

            end do
         end if
      end do
      if ( act(:5)/='     ' ) then
 
         ! ... process the active set
         !     the member of the input is limited to nels .
         !
         !     nq   =  number of configurations
         !     ela  =  eli for active set
         !     qa   =  qi for active set
         !     ma   =  number of eli for active set
         !     rl   =  L-value for each shell
         !     f    =  full value for each shell

         ma   = 0
         temp = act
 
         ! ... Decompose the input of Active Set, Convert the input of uppercase
         !     to lowercase, assign values to ela,qa,rl,f .
         do while ( ma<nels )
            if ( temp(1:5)=='     ' ) then
 
               ! ... nels is maximun value of active electrons
               do i = ma + 1, nels
                  f(i) = 0
               end do
               do
 
               ! ... Generate other possible configurations from the given set
               !     f(i) is the maximum allowed number of electrons
               !     h(i) is the number of unassigned electrons or ``holes''
               if ( itype==0 ) then
                     minq1 = 0
                  else if ( itype==1 ) then
                     minq1 = max0(min0(f(1),const)-1,0)
                  else if ( itype==2 ) then
                     minq1 = max0(min0(f(1),const)-2,0)
                  else if ( itype==3 ) then
                     minq1 = max0(min0(f(1),const)-3,0)
                  else
                     write (0,*) ' Unknown type: Re-enter'
                     read *, itype
                     cycle
                  end if
                  do q1 = min0(f(1),const), minq1, -1
                     h(1) = max0(0,const-q1)
                     do q2 = min0(f(2),h(1)), 0, -1
                        h(2) = max0(0,h(1)-q2)
                        do q3 = min0(f(3),h(2)), 0, -1
                           h(3) = max0(0,h(2)-q3)
                           do q4 = min0(f(4),h(3)), 0, -1
                              h(4) = max0(0,h(3)-q4)
                              do q5 = min0(f(5),h(4)), 0, -1
                                 h(5) = max0(0,h(4)-q5)
                                 do q6 = min0(f(6),h(5)), 0, -1
                                    h(6) = max0(0,h(5)-q6)
                                    do q7 = min0(f(7),h(6)), 0, -1
                                       h(7) = max0(0,h(6)-q7)
                                       do q8 = min0(f(8),h(7)), 0, -1
                                          h(8) = max0(0,h(7)-q8)
                                          do q9 = min0(f(9),h(8)), 0, -1
                                             h(9) = max0(0,h(8)-q9)
                                             do q10 = min0(f(10),h(9)), 0, -1
                                                h(10) = max0(0,h(9)-q10)
                                                do q11 = min0(f(11),h(10)), 0,  &
                                                 & -1
                                                   h(11) = max0(0,h(10)-q11)
                                                   do q12 = min0(f(12),h(11)),  &
                                                    & 0, -1
                                                      h(12) = max0(0,h(11)-q12)
                                                      do q13 = min0(f(13),h(12))&
                                                       & , 0, -1
                                                         h(13)                  &
                                                          & = max0(0,h(12)-q13)
                                                         do q14 = min0(f(14),   &
                                                          & h(13)), 0, -1
                                                           h(14) = max0(0,h(13)-&
                                                            & q14)
                                                           if ( h(14)<=f(15) )  &
                                                            & then
                                                           q15 = h(14)
                                                           call s_config
                                                           end if
                                                         end do
                                                      end do
                                                   end do
                                                end do
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do

                  ! ... Print Header of the output file
                  
                  if ( nq==0 ) go to 1000
                  write (6,10022)
10022             format (//                                                    &
             &'      GENERATE ALL POSSIBLE CONFIGURATIONS FROM THE   ACTIVE SET'&
             & //)
                  j = index(header,'         ')
                  write (6,10023) header(1:j)
10023             format (' '/t10,'-------------       ',a,'  --------'/)
                  write (6,10024) shells
10024             format ('          Closed Shells  :  ',a60/)
                  write (6,10025) act
10025             format ('             Active Set  :  ',a60/)
 
                  ! ... Print all configurations generated from Active Set

                  write (6,*) '         Configurations  :'
                  do i = 1, nq
                     read (file1(i),10034) (qa(j),j=1,ma)
                     if ( ela(1)(1:1)==' ' ) then
                        n = 27
                     else
                        n = 28
                     end if
 
                     k = 0
                     do j = 1, ma
                        if ( qa(j)/=0 .or. ma==2 ) then
                           k      = k + 1
                           elb(k) = ela(j)
                           qb(k)  = qa(j)
                        end if
                     end do
                     write (6,10026) (elb(j),qb(j),j=1,k)
10026                format (t28,8(1x,a3,'(',i2,')'))
                  end do
 
                  ! ... For each new configuration, generate all possible couplings
                  n1 = 0
                  n2 = 0
                  do i = 1, nq
                     read (file1(i),10034) (qa(j),j=1,ma)
 
                     ! ... omit ela(i) if corresponding qa(i)=0
                     k = 0
                     do j = 1, ma
                        if ( qa(j)/=0 .or. ma==2 ) then
                           k      = k + 1
                           elb(k) = ela(j)
                           qb(k)  = qa(j)
                        end if
                     end do
 
                     ! ... omit configurations which have more than 8 shells
                     if ( k<=8 ) then
                        call s_coupld(elb,qb,k,elc,nc,*1100)
                        if ( nc>0 ) then
                           write (6,10027)
10027                      format (//'    LIST THE COUPLINGS FOR EACH',         &
                                  &' CONFIGURATION GENERATED BY THE ACTIVE SET' &
                                  & //)
                           call s_print(elc,qb,k,nc,2)
                        else
                           n2 = n2 + 1
                        end if
                     else
                        n1 = n1 + 1
                     end if
                  end do
                  if ( n1/=0 ) print 10028, n1
10028             format (t5,'Too many occuplied shells --- ',i3,               &
                         &' configuration omitted!')
                  if ( n2/=0 ) print 10029, n2
10029             format (t5,'No final term as your selection for ',i3,         &
                         &' Active set!')
 
                  if ( nrepl/=0 ) go to 1000
                  go to 1100
               end do
            else
               call s_del(temp)
               ma = ma + 1
               n  = ichar(temp(2:2))
               
               if ( n>=ordua .and. n<=orduz ) then
                  temp(2:2) = char(n-ordua+ordla)
               end if

               ela(ma) = temp(:2)
               rl(ma)  = f_lval(temp(2:2))
               
               if ( rl(ma)>3 ) then
                  f(ma) = 2
               else
                  f(ma) = min0(const,4*rl(ma)+2)
               end if

               call s_strsh(temp,4)

            end if
         end do
         write (0,*) 'Too many electrons in the active set: MAX= 15'
         stop
      end if
 
      ! ... process the replacements
      !
      !     strl  :  string of the left side of '='
      !     strr  =  string of the right side of '='
 
 1000 if ( sflag/=0 .or. dflag/=0 ) then
 
         ! ... process the virtual set
         !
         !     rl  =  parity of eli of reference set
         !     vl  =  parity of eli of virtual set

         do nf = 1, nref
 
            ! ... set the initial values of reference set
            m = ms(nf)
            do i      = 1, nels
               ch3    = els(i)
               el(i)  = ch3
               elb(i) = ch3
               q(i)   = qs(i,nf)
               qb(i)  = q(i)
               rl(i)  = mod(f_lval(ch3(2:2)),2)
            end do
 
            if ( sflag/=0 .or. dflag==0 ) then
 
               ! ... preparation for the single replacement
               ql(1) = 1
               ml    = 1
               qr(1) = 1
               mr    = 1
               do i = nvir, lvir
                  if ( q(i)/=0 ) then
                     ell(1) = el(i)
                     if ( rl(i)>0 ) then
                        n = 1
                     else
                        n = 0
                     end if
                     do j = 1, mv
 
                     ! ... Replace Qi of Reference Set which has the same parity
                     !     with Qj of Virtual Set
                     if ( vl(j)==n ) then
                           elr(1) = elv(j)
                           call s_replace(elb,qb,k,*1002)
                           call s_coupld(elb,qb,k,elc,nc,*1100)

                           if ( nc>0 ) then
                              write (6,10030)
10030                         format (//'        FOR VIRTUAL SET, ',            &
                       &'GENERATE CONFIGURATION AND COUPLINGS FOR S-REPLACEMENT'&
                      & //)
                              nr       = nels
                              temp     = ell(1)//' = '//elr(1)
                              call s_del(temp)
                              repl(nr) = temp
                              call s_print(elc,qb,k,nc,4)

                           end if
                        end if
 1002                end do
                  end if
               end do
            end if
 
 
! ... Replace pairs of q(i) and q(j) by double virtual set
!     pl  =  pairty for the pair qi and qj in reference set

            if ( dflag==0 ) exit
            ml    = 2
            ql(1) = 1
            ql(2) = 1
            
            do i = nvir, lvir - 1
               ell(1) = el(i)
               ll1    = f_lval(ell(1)(2:2))
               
               do j = i + 1, m
                  if ( q(i)/=0 .and. q(j)/=0 ) then
                     ell(2) = el(j)
                     ll2    = f_lval(ell(2)(2:2))
                     temp   = ell(1)//'.'//ell(2)//' = '
                     llmin  = iabs(ll1-ll2)
                     llmax  = ll1 + ll2
                     pl     = mod(llmax,2)
                     call s_vpair(elv,mv,pl,llmin,llmax,temp,*1100)
                  end if
               end do
            end do
 
            ! ... Replace pairs of (Qi)=2 by Double Virtual Set
            ml    = 1
            ql(1) = 2
            do i = nvir, lvir
               if ( q(i)>1 ) then

                  ell(1) = el(i)
                  ll1    = f_lval(ell(1)(2:2))
                  llmin  = 0
                  llmax  = ll1 + ll1
                  temp   = ell(1)//'(2) = '
                  pl     = mod(llmax,2)
                  call s_vpair(elv,mv,pl,llmin,llmax,temp,*1100)
               end if
            end do
         end do
      else
         n1 = 0
         do nr = 1, nrepl
            strr = repl(nr)
            j    = index(strr,'=')
            strl = strr(1:j-1)
            call s_strsh(strr,j+1)
 
            ! ... Decompose the substring on the left of ''=''
            !     ell  =  old value of eli to be replaced
            !     ql   =  old value of qi to be replaced
            !     ml   =  number of old value of eli to be replaced

            call s_del(strl)
            call s_decomp(strl,ell,ql,ml)
 
            ! ... Decompose the substring on the right of ''=''
            !
            !     elr  =  new value of eli
            !     qr   =  new varue of qi
            !     mr   =  number of the new value of eli

            call s_del(strr)
            call s_decomp(strr,elr,qr,mr)
 
           ! ... For each Replacement, replace all Reference Sets

            do nf = 1, nref
               m = ms(nf)
               if ( m>=ml ) then
                  do i = 1, m
                     el(i) = els(i)
                     q(i) = qs(i,nf)
                  end do

                  call s_replace(elb,qb,k,*1020)
                  call s_coupld(elb,qb,k,elc,nc,*1100)

                  if ( nc>0 ) then
                     write (6,10031)
10031                format (//'      FOR EACH REPLACEMENT, GENERATE',          &
                          &' CONFIGURATIONS AND COUPLINGS FOR THE REFERENCE SET'&
                         & //)
                     call s_print(elc,qb,k,nc,3)
                  else
                     n1 = n1 + 1
                  end if
               end if
 1020       end do
         end do
         if ( n1/=0 ) write (0,10032) n1
10032    format (t5,'No Final Term as selected for ',i3,' Replacements!')
      end if
 
      ! ... Program ending

 1100 write (7,'(a1)') '*'
      close (7)
      
      write (0,*)
      write (0,*) '         OK !'
      write (0,*) '         List of configurations and their couplings'
      write (0,*) '         is in the file ', fname

10033 format (t29,i2,'  :  ',a60/)
10034 format (21(i2))
      
end program gencl_main
