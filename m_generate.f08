module m_generate
   
   ! Module with routines for generating configuration state functions
   ! =================================================================

   use m_parameters
   implicit none
   
   ! List of routines in order of appearance
   !
   ! type         name           description
   ! -----------  -------------  -----------------------------------------------
   ! subroutine   s_config       tests electron number & parity of new config.
   ! subroutine   s_coupld       generates all possible couplings
   ! subroutine   s_replace      generates electron replacements   
   ! subroutine   s_vpair        generates occ. or virt. pairs for D-replacement
   ! ---------------------------------------------------------------------------

contains

   subroutine s_coupld(e_lbl,q_numb,n_shells,e_lbl_o,nc,*)

      ! This subroutine generates all possible couplings
      ! ================================================

      !   first, compute alpha value from the given configuration,
      !   then compute beta from each value of alpha

      !   input :
      !        e_lbl =  electron labels
      !                 where e_lbl(i)(1=1)  ---  blank
      !                       e_lbl(i)(2=2)  ---  n-symbol
      !                       e_lbl(i)(3=3)  ---  L-symbol
      !       q_numb =  occupation number
      !                       0 (empty) <= q_numb(i)  <= 2(2L(i)+1) (full)
      !     n_shells =  number of shells
      !                       0  <=  n_shells  <=  5
      !   output :
      !  e_lbl_o     =  electron labels (output)     
      !       nc     =  number of couplings
      !        *     =  return label if the maximun number of couplings > nscoup
      !

      use m_globals
      use m_tools, only: f_lval, f_ictoi, f_symb

      implicit none

      ! dummy arguments
      character(3), dimension(*), intent(in)    :: e_lbl
      integer(i4),  dimension(*), intent(in)    :: q_numb
      integer(i4),                intent(in)    :: n_shells
      character(3), dimension(*), intent(inout) :: e_lbl_o
      integer(i4),                intent(inout) :: nc 

      ! local variables

      real(r8)       :: full, half, s, s1, s2, s3, s4, kr
      character(2)   :: a2, b1, ch2
      character(30)  :: calfa
      character      :: cch1, ch1
      character(3)   :: cch3, ch3
      character(11)  :: coupling_fmt ! the coupling format (the denser atsp2k style or the old style)
      character(65), save :: termq3
      character(3),  dimension(nshel,nscoup)  :: alfa
      character(3),  dimension(ncoupl,nscoup) :: couple
      integer(i4),   dimension(nels)          :: beta, posit
      integer(i4),   dimension(10), save      :: lposit
      integer(i4),   dimension(27), save      :: nterm
      character(65), dimension(27), save      :: term,  term2, term3, &
                                                 term4, term5, term6
      integer(i4) :: child, i, j, jj, k, ks, l, l1, l2, l3, l4, lmax, lmin, &
               & loca, locb, loct, mbeta, n, nalfa, nb, nbeta, nt, parent, ptr


      ! ... Number of possible terms for configurations P(1-3),
      !     D(1-5), F(1-2), G(1-2), ... M(1-2)
      data (nterm(i),i=1,27)/1, 3, 3, 1, 5, 8, 16, 16, 1, 7, 17, 47, 73, 119,   &
          & 119, 1, 9, 1, 9, 1, 9, 1, 9, 1, 9, 1, 9/

      ! ... Starting position in term table for given L
      data (lposit(i),i=1,9)/1, 4, 9, 16, 18, 20, 22, 24, 26/

      ! ... Possible terms for configurations P(1-3),D(1-5),F(1-2),G(1-2)
      data (term(i),i=1,27)/'2P1', '1S01D23P2', '2P12D34S3', '2D1',             &
           &'1S01D21G23P23F2', '2D12P32D32F32G32H34P34F3',                      &
           &'1S01D21G23P23F21S41D41F41G41I43P43D43F43G43H45D4',                 &
           &'2D12P32D32F32G32H34P34F32S52D52F52G52I54D54G56S5', '2F1',          &
           &'1S01D21G21I23P23F23H2',                                            &
           &'2P12D12D22F12F22G12G22H12H22I12K12L14S14D14F14G14I1',              &
           &'1S11S21D11D21D31D41F11G11G21G31G41H11H21I11I21I31K11L11L21N13P1',  &
           &'2P12P22P32P42D12D22D32D42D52F12F22F32F42F52F62F72G12G22G32G42G5',  &
           &'1S11S21S31S41P01D11D21D31D41D51D61F11F21F31F41G11G21G31G41G51G6',  &
           &'2S12S22P12P22P32P42P52D12D22D32D42D52D62D72F12F22F32F42F52F62F7',  &
           &'2G1', '1S01D21G21I21L23P23F23H23K2', '2H1',                        &
           &'1S01D21G21I21L23P23F23H23K2', '2I1', '1S01D21G21I21L23P23F23H23K2',&
           &'2K1', '1S01D21G21I21L23P23F23H23K2', '2L1',                        &
           &'1S01D21G21I21L23P23F23H23K2', '2M1', '1S01D21G21I21L23P23F23H23K2'/
      data (term2(i),i=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',  &
           &'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',                          &
           &'3P23P33D13D23F13F23F33F43G13G23G33H13H23H33H43I13I23K13K23L13M1',  &
           &'2G62H12H22H32H42H52H62H72I12I22I32I42I52K12K22K32K42K52L12L22L3',  &
           &'1G71G81H11H21H31H41I11I21I31I41I51I61I71K11K21K31L11L21L31L41M1',  &
           &'2F82F92FA2G12G22G32G42G52G62G72G82G92GA2H12H22H32H42H52H62H72H8',  &
           &'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',     &
           &'ERR', 'ERROR', 'ERR', 'ERROR'/
      data (term3(i),i=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',  &
           &'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', '5S05D15F15G15I1',       &
           &'2M12M22N12O04S14P14P24D14D24D34F14F24F34F44G14G24G34G44H14H24H3',  &
           &'1M21N11N21Q03P13P23P33P43P53P63D13D23D33D43D53F13F23F33F43F53F6',  &
           &'2H92I12I22I32I42I52I62I72I82I92K12K22K32K42K52K62K72L12L22L32L4',  &
           &'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',     &
           &'ERR', 'ERROR', 'ERR', 'ERROR'/
      data (term4(i),i=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',  &
           &'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', 'ERROR',                 &
           &'4I14I24I34K14K24L14M06P06F06H0',                                   &
           &'3F73F83F93G13G23G33G43G53G63G73H13H23H33H43H53H63H73H83H93I13I2',  &
           &'2L52M12M22M32M42N12N22O02Q04S14S24P14P24D14D24D34D44D54D64F14F2',  &
           &'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',     &
           &'ERR', 'ERROR', 'ERR', 'ERROR'/
      data (term5(i),i=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',  &
           &'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', 'ERROR', 'ERROR',        &
           &'3I33I43I53I63K13K23K33K43K53K63L13L23L33M13M23M33N03O07F05S05P0',  &
           &'4F34F44F54G14G24G34G44G54G64G74H14H24H34H44H54I14I24I34I44I54K1',  &
           &'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR',     &
           &'ERR', 'ERROR', 'ERR', 'ERROR'/
      data (term6(i),i=1,27)/'ERR', 'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR',  &
           &'ERROR', 'ERROR', 'ERR', 'ERROR', 'ERROR', 'ERROR', 'ERROR',        &
           &'5D15D25D35F15F25G15G25G35H15H25I15I25K05L0',                       &
           &'4K24K34L14L24L34M04N06P06D06F06G06H06I08S0', 'ERR', 'ERROR', 'ERR',&
           &'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR', 'ERROR', 'ERR',     &
           &'ERROR'/
      data termq3/'1S11D11G11I13P13F13H1'/

      ! ... define coupling printing format variable

      if ( answ_coupl_fmt_new == 'n' ) then
         ! atsp pre-2k (such as the "bookcodes") coupling formatting
         coupling_fmt = '(15(5x,a3))'
      else
         ! atsp2k coupling formatting (denser)
         coupling_fmt = '(15(1x,a3))'
      end if

      ! ... compute the possible values for alpha
      !
      !     empty  =  1 when q_numb(i)=0
      !     full   =  2*(2L+1)
      !     half   =  2L+1
      !     alfa   =  matrix of (NELS,NSCOUP)
      !     nalfa  =  Number of ALFA
      !     posit  =  array of NELS, store the position in table NTERM and TERM
      !     corresponding q_numb(i). Rule  :
      !     position = (L-1)*2+q_numb(i)        if 1 <= q_numb(i) <= HALF
      !     position = (L-1)*2+(FULL-q_numb(i))   if q_numb(i) > HALF

      nalfa = 1
      do i = 1, n_shells
         ch3 = e_lbl(i)
         ch1 = ch3(2:2)
         if ( ch3(3:3)==' ' ) then
            cch3 = ch3
            ch3 = ' '//cch3(1:2)
         end if
         e_lbl_o(i) = ch3
         cch1 = char(ichar(ch1)-ordla+ordua)
         ch1 = cch1
         l = f_lval(ch1)
         full = 4*l + 2
         k = q_numb(i)

         ! ... If shell is full, ALFA(i) = 1S0
         if ( k==0 .or. k==int(full) ) then
            ch3 = '1S0'

         ! ... If q_numb(i) = 1,  then ALFA(i)=2<L-symbol>1
         else if ( k==1 ) then
            ch3 = '2'//ch1//'1'

         ! ... Otherwise, get the possible value from array nterm and term
         else
            half = full/2.0_r8
            if ( k<=half ) then
               posit(i) = lposit(l) + k - 1
            else
               if( abs(real(int(full,i4),r8)-full) > 0._r8 ) then
                  print*, ' -- WARNING! (m_generate, row 187): real8 variable "full" should be integer valued!'
               end if
               posit(i) = lposit(l) + (int(full,i4)-k) - 1
            end if
            nalfa = nalfa*nterm(posit(i))
            ch3 = '   '
         end if
         if ( nalfa>nscoup ) then
            write (0,*) 'Array ALFA in routine COUPLD exceeded'
            stop
         end if

         ! ... calfa is a string storing the elements in one alfa
         calfa(3*i-2:3*i) = ch3
      end do

      ! ... assign values to all elements of alfa
      !
      !     nt  =
      !     loct  =  current position in the table term
      !     loca  =  current position in the matrix alfa

      nt = 1
      do i = n_shells, 1, -1
         ch3 = calfa(i*3-2:i*3)
         if ( ch3/='   ' ) then
            do j = 1, nalfa
               alfa(i,j) = ch3
            end do
         else
            loct = posit(i)
            n    = nterm(loct)
            loca = 1
            do
               if ( loca<=nalfa ) then
                  do j = 1, n

                     !CH3 = TERM(LOCT)(J*4-3:J*4-1) ! -- remove this??

                     ! ... the case f(12)
                     if ( loct==10 .and. q_numb(i)==12 ) then
                        ch3 = termq3(j*3-2:j*3)

                     !... other cases
                     else if ( j<22 ) then
                        ch3 = term(loct)(j*3-2:j*3)
                     else if ( j<43 ) then
                        jj  = j - 21
                        ch3 = term2(loct)(jj*3-2:jj*3)
                     else if ( j<64 ) then
                        jj  = j - 42
                        ch3 = term3(loct)(jj*3-2:jj*3)
                     else if ( j<85 ) then
                        jj  = j - 63
                        ch3 = term4(loct)(jj*3-2:jj*3)
                     else if ( j<106 ) then
                        jj  = j - 84
                        ch3 = term5(loct)(jj*3-2:jj*3)
                     else if ( j<120 ) then
                        jj  = j - 105
                        ch3 = term6(loct)(jj*3-2:jj*3)
                     end if
                     do k = 1, nt
                        alfa(i,loca) = ch3
                        loca         = loca + 1
                     end do
                  end do
                  cycle
               end if
               nt = nt*n
               exit
            end do
         end if
      end do

      ! ... GENERATE POSSIBLE VALUE OF BETA FROM ALFA

      nc = 0
      do nb = 1, nalfa

         ! ... There is only one coupling if n_shells = 1
         if ( n_shells==1 ) then
            beta(1) = 1
            couple(1,1) = alfa(1,nb)
            nbeta = 1
            go to 50
         end if

         ! ... Define BETA(1)=ALFA(1), then the next basic coupling steps is :
         !     S1 = (BETA(1)(1:1)-1)/2 ,   S2 = (ALFA(2)(1:1)-1)/2 ;
         !     | S1-S2 | <= BETA(2)(1:1) <= | S1+S2 |
         !     L1 = L-number of BETA(1)(2:2),   L2 = L-number of ALFA(2)(2:2) ;
         !     Symbol(| L1-L2 |) <= BETA(2)(2:2) <= Symbol(| L1+L2 |) .

         b1 = alfa(1,nb)(1:2)
         do j = 2, n_shells
            a2      = alfa(j,nb)(1:2)
            parent  = 1
            child   = 1
            beta(j) = 0
            do
               s1 = (ichar(b1(1:1))-ord0-1)/2.
               s2 = (ichar(a2(1:1))-ord0-1)/2.
               s3 = abs(s1-s2)
               s4 = abs(s1+s2)
               l1 = f_lval(b1(2:2))
               l2 = f_lval(a2(2:2))
               l3 = abs(l1-l2)
               l4 = abs(l1+l2)

               mbeta = ( int(s4-s3,i4)+1 )*( l4-l3+1 )

               ! internal check, make sure s4-s3 actually is an integer
               if( abs(real(int(s4-s3,i4),r8)-(s4-s3)) > 0._r8 ) then
                  print*, ' -- WARNING! (m_generate, row 301): real8 value "(s4-s3)" should be integral!'
               end if

               ! ... Generate BETA from each ALPHA

               ! There are four scratch files for storing the information about
               ! Beta(i), 1<i<6, shown as follows :
               ! --------------------------------------------
               ! |    Parent    |    Value    |    Child    |
               ! --------------------------------------------
               ! Define   Beta(i) is child of Beta(i-1) and parent of beta(i+1) ;
               ! Parent is a pointer to the parent of Beta(i), that is Beta(i-1) ;
               ! Value is one of the possible value for BETA(i), and Child is the
               ! number of children for Beta(i) .
               !
               ! MBETA  =  Number of couplings generated from ALFA(i)
               ! PARENT  =  Current pointer to the parent of Beta(i)

               ! The following loop was over non-integers prior to the f90 translation, which is illegal in modern fortran.
               ! Now instead we loop over floored values and increment the real valued kr from s3 with 1 inside the loop
               kr  = s3
               do k = floor(s3,i4), floor(s4,i4)
                  ks = int(2._r8*kr, i4)
                  ch1 = char(ord0+1+ks)
                  do l = l3, l4
                     write (fbeta(j-1,child),10000) parent, ch1//f_symb(l), 1
                     child = child + 1
                  end do
                  kr  = s3 + 1._r8
               end do

               beta(j) = beta(j) + mbeta
               if ( j==2 ) exit

               ! ... Correct the number of its children for Beta(j-1)
               locb = parent
               do k = j - 1, 2, -1
                  read (fbeta(k-1,locb),10000) ptr, ch2, n
                  write (fbeta(k-1,locb),10000) ptr, ch2, n + mbeta - 1
                  locb = ptr
               end do

               ! ... If pointer to the parent is the end of file for Beta(j-1),
               !     prepare to generate Beta(i+1) ; otherwise, generate next Beta(j)
               !     according Alfa(j-1) and Beta(j-1) .
               parent = parent + 1
               if ( parent>beta(j-1) ) exit
               read (fbeta(j-2,parent),10000) ptr, b1, n
            end do
            read (fbeta(j-1,1),10000) ptr, b1, n
         end do

         ! ... Assign values to the couplings forward ,
         !     COUPLE(I) = Alpha(i) for coupling(j) if i <= n_shells ;
         nbeta = beta(n_shells)
         if ( nbeta>nscoup ) then
            write (0,*) 'Array COUPLE in routine COUPL exceeded'
            stop
         end if
         do j = 1, nbeta
            do i = 1, n_shells
               couple(i,j) = alfa(i,nb)
            end do
            do k = 2*n_shells - 1, 9
               couple(k,j) = '   '
            end do
         end do

         ! ... Assign values to the couplings backward ,
         !     COUPLE(i) = Beta(i-n_shells) for coupling(j) if i > n_shells
         do i = n_shells, 2, -1
            n = 1
            do j = 1, beta(i)
               read (fbeta(i-1,j),10000) ptr, ch2, nt
               do k = 1, nt
                  couple(n_shells+i-1,n) = ch2//'0'
                  n = n + 1
                  if ( n>nscoup ) then
                     write (0,*) 'Array COUPLE exceeded'
                     stop
                  end if
               end do
            end do
         end do

         ! ... Selection from generated couplings according the following rules
 50      do i = 1, nbeta
            n = 2*n_shells - 1

            ! ... If the first time to call COUPLD, not compute MAX and MIN
            if ( maxfinal/=-5 ) then

               ! ... Compute MAX and MIN value for each final term, keep it if
               !     intersection is non-empty
               ch2  = couple(n,i)(1:2)
               ch3  = ch2(1:1)
               ch1  = ch2(2:2)
               s    = (f_ictoi(ch3)-1)/2._r8
               l    = f_lval(ch1)
               lmin = int( 2._r8 * abs(s-real(l,r8)), i4)
               lmax = int( 2._r8 * abs(s+real(l,r8)), i4)

               if ( lmin>maxfinal .or. lmax<minfinal ) cycle
            end if

            ! ... If Final Terms are given, do selection
            if ( nftm/=0 ) then
               ch2 = couple(n,i)(:2)
               do k = 1, nftm
                  if ( ch2==ftm(k) ) go to 60
               end do
               cycle
            end if

            ! ... Waining if the number of couplings > NSCOUP
            if ( nc==nscoup ) then
               write (0,*) '          WARNING !'
               write (0,*) '          The number of couplings',                 &
                          &' is greater than 2500 .',                           &
                          &'          Please select the Final Term .'
               return 1
            end if

            ! ... Write configurations and couplings to CI.LST
 60         nc = nc + 1
            write (file3(nc),'(9(a3))') (couple(j,i),j=1,n)
            write (7,10001) (e_lbl_o(j),q_numb(j),j=1,n_shells)
            write (7,trim(coupling_fmt)) (couple(j,i),j=1,n)

         end do
      end do

10000 format (i3,a2,i3)

      ! subshell and occupation number
10001 format (5(1x,a3,'(',i2,')'))


   end subroutine s_coupld
   
   subroutine s_replace(e_lbl,q_numb,n_shells,*)

      ! Routine for replacing electrons
      ! ===============================
      
      ! input/output
      !       e_lbl  =  number of shells
      !      q_numb  =  electon label array
      !    n_shells  =  occupation number array

      use m_globals
      use m_tools, only: f_lval
      implicit none
      
      ! dummy arguments
      integer(i4),                   intent (inout)  :: n_shells 
      character(3), dimension(nels), intent (inout)  :: e_lbl    
      integer(i4),  dimension(nels), intent (inout)  :: q_numb   
     
      ! local variables
      character(3)                  :: ch3
      character(3), dimension(nels) :: elc
      integer(i4)                   :: i, j, k, l, mark, mc, n, np, qq
      integer(i4),  dimension(ncfg) :: pp
 
      do i = 1, nels
         elc(i) = el(i)
         qc(i) = q(i)
      end do
 
      ! ... Correct Q(i) by subtraction
      do i = 1, ml
         ch3 = ell(i)
         mark = 1
         do j = 1, m
            if ( ch3==elc(j) ) then
               qc(j) = qc(j) - ql(i)
               if ( qc(j)>=0 ) mark = 0
            end if
         end do
      end do
      if ( mark/=0 ) return 1
 
      ! ... Correct QC(i) by adding
      mc = m
      do i = 1, mr
         ch3 = elr(i)
         mark = 0
         do j = 1, mc
            if ( ch3==elc(j) ) then
               qc(j) = qc(j) + qr(i)
               mark = 1
            end if
         end do
         if ( mark==0 ) then
            mc = mc + 1
            elc(mc) = ch3
            qc(mc) = qr(i)
         end if
      end do
 
      ! ... Delete EL(i) if Q(i) = 0
      n_shells = 0
      do i = 1, mc
         if ( qc(i)/=0 .or. mc==2 ) then
            n_shells = n_shells + 1
            e_lbl(n_shells) = elc(i)
            q_numb(n_shells) = qc(i)
         end if
      end do
 
      ! ... Check the input error after replacement
      j = 0
      k = 0
      do i = 1, n_shells
         ch3 = e_lbl(i)
         l = f_lval(ch3(2:2))
         if ( q_numb(i)>l*4+2 ) return 1
         j = j + q_numb(i)
         k = k + q_numb(i)*l
      end do
 
      if ( j/=const ) return 1
 
      if ( mod(k,2)/=parityval ) return 1
 
      ! ... If the replacement duplicates a configuration in the
      !     active set, it should not be sent to CI.LST
      do i = 1, nq
         read (file1(i),10001) (qa(j),j=1,ma)
10001    format (15(i2))
         n = 0
         do j = 1, n_shells
            ch3 = e_lbl(j)
            do k = 1, ma
               if ( ch3==ela(k) .and. q_numb(k)==qa(k) ) n = n + 1
            end do
         end do
         if ( n==n_shells ) return 1
      end do
 
      ! ... If the replacement duplicates a configuration in the
      !     Reference Set, it should not be sent to CI.LST
      do i = 1, nref
         l = ms(i)
         n = 0
         do j = 1, l
            ch3 = els(j)
            qq = qs(j,i)
            do k = 1, n_shells
               if ( ch3==e_lbl(k) .and. qq==q_numb(k) ) n = n + 1
            end do
         end do
         if ( n==n_shells ) return 1
      end do
 
      ! ... If the replacement duplicates a configuration in the previous
      !     replacement, it should not sent to CI.LST
      print*, np
      do i = 1, np
         l = pp(i)
         read (file2(i),10002) (elc(j),qc(j),j=1,l)
10002    format (8(a3,i2))
         n = 0
         do j = 1, l
            ch3 = elc(j)
            do k = 1, n_shells
               if ( ch3==e_lbl(k) .and. qc(k)==q_numb(k) ) n = n + 1
            end do
         end do
         if ( n==n_shells ) return 1
      end do
 
      np     = np + 1
      pp(np) = n_shells
      
      write (file2(np),10003) (e_lbl(j),q_numb(j),j=1,n_shells)

10003 format (8(a3,i2))

   end subroutine s_replace

   subroutine s_vpair(elv,mv,pl,llmin,llmax,str,*)

      ! Generate occupied or virtual pairs for D-Replacement
 
      ! input
      !
      !   elv  =  eli for virtual set
      !    mv  =  number of eli for virtual set
      !    pl  =  parity of qi for reference set
      ! llmin  =  minimum angular coupling of the pair
      ! llmax  =  maximum angular coupling of the pair
      !   str  =  string to be packed as output for replacement
 
      use m_globals
      use m_tools, only: f_lval, s_del
      use m_io, only: s_print
      
      implicit none
      
      ! dummy arguments
      integer(kind=i4), intent(in)  :: llmax, llmin, mv, pl
      character(72),    intent(in)  :: str
      character(3), dimension(nels), intent(in) :: elv
     
      ! local variables
      character(3), dimension(nels)     :: elb, elc
      character(72)                     :: sstr
      integer(kind=i4)                  :: i, j, k, lr1, lr2, lrmax, lrmin, n, nc, pr
      integer(kind=i4), dimension(nels) :: q_numb
 
      n  = index(str,'=')
      nr = nels
 
      ! ... D-Replacement for the pair of two single ELi
      mr    = 2_i4
      qr(1) = 1_i4
      qr(2) = 1_i4

      do i = 1_i4, mv - 1_i4
         elr(1) = elv(i)
         lr1    = f_lval(elr(1)(2:2))

         do j = i + 1, mv
            elr(2) = elv(j)
            lr2    = f_lval(elr(2)(2:2))
            lrmin  = iabs(lr1-lr2)
            lrmax  = lr1 + lr2
            pr     = mod(lrmax,2)
            
            ! ... If the pair has the same parity with the left side, and the
            !     angular coupling of the two pairs have values in common, replace
            !     them, then generate couplings for the new configuration
            
            if ( pl==pr .and. lrmin<=llmax .and. lrmax>=llmin ) then
               
               call s_replace(elb,q_numb,k,*50)
               call s_coupld(elb,q_numb,k,elc,nc,*200)

               if ( nc>0 ) then
               
                  write (6,*) '   FOR VIRTUAL SET, GENERATE CONFIGURATION &
                            & AND COUPLINGS FOR D-REPLACEMENT'
                  
                  sstr     = str(:n)//elr(1)//'.'//elr(2)
                  call s_del(sstr)
                  repl(nr) = sstr
                  call s_print(elc,q_numb,k,nc,4)
               
               end if
            end if
 50      end do
      end do
 
      ! ... D-Replacement for the pairs which has the value Qi=2
      mr    = 1_i4
      qr(1) = 2_i4

      do i = 1, mv
         elr(1)   = elv(i)
         lr1      = f_lval(elr(1)(2:2))
         lrmax    = lr1 + lr1
         pr       = mod(lrmax,2)

         ! ... If it has the same parity with the left side, replace them,
         !     then generate couplings for the new configuration
         if ( pl==pr .and. lrmax>=llmin ) then
            call s_replace(elb,q_numb,k,*100)
            call s_coupld(elb,q_numb,k,elc,nc,*200)

            if ( nc>0 ) then
               
               write (6,*) '   FOR VIRTUAL SET, GENERATE CONFIGURATION &
                             & AND COUPLINGS FOR D-REPLACEMENT'
               
               sstr = str(:n)//elr(1)//'(2)'
               call s_del(sstr)
               repl(nr) = sstr
               call s_print(elc,q_numb,k,nc,4)

            end if
         end if
 100  end do
 
      return
 200  return 1

   end subroutine s_vpair

   subroutine s_config
      
      ! Examine if the new configuration has the same electron number and parity
      
      ! input
      !  - q, qi =  occupation number
      !  - l     =  L-value corresponding eli
      ! output
      !  - nq    =  number of configurations
      
      use m_globals
      
      implicit none
      
      ! local variables
      character(3) :: ch3
      integer(i4)  :: i, j, k, n, iq, newp
 
      ! return if the new configuration has different pairty
      newp = 0
      do i = 1, ma
         newp = newp + qa(i)*rl(i)
      end do
 
      newp = mod(newp,2)
      if ( newp/=parityval ) return
 
      ! return if the new condiguration is the same as the Reference Set
      do i = 1, nref
         m = ms(i)
         n = 0
         do j = 1, m
            ch3 = els(j)
            iq  = qs(j,i)
            do k = 1, ma
               if ( ch3==ela(k) .and. iq==qa(k) ) n = n + 1
            end do
         end do
         if ( n==m ) return
      end do
 
      ! Otherwise, write them into the configuration file
      nq = nq + 1
      write (file1(nq),'(15(i2))') (qa(j),j=1,ma)

   end subroutine s_config

end module m_generate
