!! Convenient definitions of numerical constants.  This should not
!! be considered an exhaustive list of such constants.
module bl_constants_module

  use bl_types
  use bl_fort_module
  implicit none

  real(kind = bl_real), parameter :: ZERO    =  0.0_bl_real
  real(kind = bl_real), parameter :: ONE     =  1.0_bl_real
  real(kind = bl_real), parameter :: TWO     =  2.0_bl_real
  real(kind = bl_real), parameter :: THREE   =  3.0_bl_real
  real(kind = bl_real), parameter :: FOUR    =  4.0_bl_real
  real(kind = bl_real), parameter :: FIVE    =  5.0_bl_real
  real(kind = bl_real), parameter :: SIX     =  6.0_bl_real
  real(kind = bl_real), parameter :: SEVEN   =  7.0_bl_real
  real(kind = bl_real), parameter :: EIGHT   =  8.0_bl_real
  real(kind = bl_real), parameter :: NINE    =  9.0_bl_real
  real(kind = bl_real), parameter :: TEN     = 10.0_bl_real

  real(kind = bl_real), parameter :: ELEVEN  = 11.0_bl_real
  real(kind = bl_real), parameter :: TWELVE  = 12.0_bl_real
  real(kind = bl_real), parameter :: FIFTEEN = 15.0_bl_real
  real(kind = bl_real), parameter :: SIXTEEN = 16.0_bl_real

  real(kind = bl_real), parameter :: HALF    = 0.5_bl_real
  real(kind = bl_real), parameter :: THIRD   = ONE/THREE
  real(kind = bl_real), parameter :: FOURTH  = 0.25_bl_real
  real(kind = bl_real), parameter :: FIFTH   = ONE/FIVE
  real(kind = bl_real), parameter :: SIXTH   = ONE/SIX
  real(kind = bl_real), parameter :: SEVENTH = ONE/SEVEN
  real(kind = bl_real), parameter :: EIGHTH  = 0.125_bl_real
  real(kind = bl_real), parameter :: NINETH  = ONE/NINE
  real(kind = bl_real), parameter :: TENTH   = 0.10_bl_real
  real(kind = bl_real), parameter :: TWELFTH = ONE/TWELVE

  real(kind = bl_real), parameter :: TWO3RD    = TWO/THREE
  real(kind = bl_real), parameter :: FOUR3RD   = FOUR/THREE
  real(kind = bl_real), parameter :: FIVE3RD   = FIVE/THREE
  real(kind = bl_real), parameter :: FIVE6TH   = FIVE/SIX

  real(kind = bl_real), parameter :: THREE4TH  = 0.75_bl_real

  real(kind = bl_real), parameter :: FIVE12TH = FIVE/TWELVE
  real(kind = bl_real), parameter :: SEVEN12TH = SEVEN/TWELVE

  real(kind = bl_real), parameter :: FIVE32ND = FIVE/32.0_bl_real

  !! Pi
  real(kind = bl_real), parameter :: M_PI    = &
       3.141592653589793238462643383279502884197_bl_real
  real(kind = bl_real), parameter :: M_SQRT_PI  = &
       1.772453850905516027298167483341145182798_bl_real

  !! Roots
  real(kind = bl_real), parameter :: M_SQRT_2  = &
       1.414213562373095048801688724209698078570_bl_real

end module bl_constants_module


module amrex_constants_module
  use bl_constants_module
  implicit none
end module amrex_constants_module
