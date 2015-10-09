module m_globals

   ! module containing global variable declarations
   ! ==============================================
      
   use m_parameters
   implicit none

   character(len=1)                       :: answ_coupl_fmt_new      ! [input data] (y/n) coupling format (new atsp2k style is denser)
   character(len=72)                      :: header                  ! [input data] header of output
   character(len=72)                      :: shells                  ! [input data] closed shells
   character(len=72)                      :: virtul                  ! [input data] virtual set
   character(len=60)                      :: act                     ! [input data] active set     (org len=60)
   character(len=60), dimension(nels)     :: ref                     ! [input data] reference set  (org len=60)
   character(len=72), dimension(nels)     :: repl                    ! [input data] replacement    (org len=60)
   character(len=2),  dimension(nels)     :: ftm                     ! [input data] final terms

   integer(kind=i4), dimension(nels)      :: q, qa, qb, qc, ql, qr   ! occupation number arrays
   integer(kind=i4), dimension(nels,nels) :: qs                      ! 
   integer(kind=i4)                       :: nq                      ! number of configurations
   integer(kind=i4)                       :: parityval               ! parity of the first reference set
   integer(kind=i4)                       :: const                   ! total number of q(i) for the first configuration
   integer(kind=i4)                       :: maxfinal                ! max final term
   integer(kind=i4)                       :: minfinal                ! min final term
   integer(kind=i4)                       :: nftm                    ! number of final terms
   integer(kind=i4)                       :: nref                    ! ?
   integer(kind=i4)                       :: nf                      ! ?
   integer(kind=i4)                       :: nr                      ! ?
   integer(kind=i4)                       :: m                       ! 
   integer(kind=i4)                       :: ma                      ! ?
   integer(kind=i4)                       :: ml                      ! ?
   integer(kind=i4)                       :: mr                      ! number of replacements

   integer(kind=i4),  dimension(nels)     :: ms                      ! ?
   integer(kind=i4),  dimension(nels)     :: rl                      ! electron L-values
  
   character(len=3),  dimension(nels)     :: el                      ! electron labels
   character(len=3),  dimension(nels)     :: ela                     ! 
   character(len=3),  dimension(nels)     :: ell                     ! 
   character(len=3),  dimension(nels)     :: elr                     ! 
   character(len=3),  dimension(nels)     :: els                     ! 

   character(len=8),  dimension(7,nscoup) :: fbeta                   ! ?
   character(len=32), dimension(ncfg)     :: file1                   ! ?
   character(len=40), dimension(ncfg)     :: file2                   ! ?
   character(len=32), dimension(nscoup)   :: file3                   ! ?

end module m_globals
