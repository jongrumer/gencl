module m_parameters

   ! module containing parameter declarations and definitions
   ! ========================================================
   
   ! note: the intrinsic module iso_fortran_env is a Fortran 2008 feature

   use, intrinsic :: iso_fortran_env

   implicit none

   ! ... integer kind parameters
   integer(integer_kinds(3)), &
                     parameter :: i4     = integer_kinds(3) ! 4  byte integer kind
   integer(kind=i4), parameter :: i8     = integer_kinds(4) ! 8  byte integer kind
   integer(kind=i4), parameter :: r4     = real_kinds(1)    ! 4  byte real kind
   integer(kind=i4), parameter :: r8     = real_kinds(2)    ! 8  byte real kind
   
   ! ... output file name
   character(len=24),parameter :: fname  = 'cfg.inp'

   ! ... other parameters
   integer(kind=i4), parameter :: nels   = 21_i4
   integer(kind=i4), parameter :: nshel  = 8_i4
   integer(kind=i4), parameter :: ncoupl = 2_i4*nshel - 1_i4
   integer(kind=i4), parameter :: ncfg   = 100000_i4
   integer(kind=i4), parameter :: nscoup = 10000_i4
   
   integer(kind=i4), parameter :: ord0  = ichar('0') ! ascii value for lower digit boundary [ 0
   integer(kind=i4), parameter :: ord9  = ichar('9') ! ascii value for upper digit boundary   9 ]
   integer(kind=i4), parameter :: ordla = ichar('a') ! ascii value for lower small   case character boundary [ a
   integer(kind=i4), parameter :: ordlz = ichar('z') ! ascii value for upper small   case character boundary  z  ]
   integer(kind=i4), parameter :: ordua = ichar('A') ! ascii value for lower capital case character boundary [ A
   integer(kind=i4), parameter :: orduz = ichar('Z') ! ascii value for upper capital case character boundary   Z ]

end module m_parameters
