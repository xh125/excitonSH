module constants
  use kinds , only : dp,dpc
  implicit none
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% SET PHYSICAL CONSTANTS                              %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/ATOMIC_UNIT       %!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/ATOMIC_MASS_UNIT  %!
!% REF: HTTP://EN.WIKIPEDIA.ORG/WIKI/PHYSICAL_CONSTANT %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!  
  integer,parameter,public        :: maxlen = 256
  real(kind=dp),parameter,public  ::             &
    au2cm = 2.194887656d5                       ,&
    au2ev = 2.7211d1                            ,&
    au2mev= 2.7211d4                            ,&
    au2j  = 4.35974417d-18                      
  real(kind=dp),parameter,public  ::             &
    au2s  = 2.418884326505d-17                  ,&
    au2fs = 2.418884326505d-2                   ,&
    au2ps = 2.418884326505d-5                  
  real(kind=dp),parameter,public  ::  au2amu= 5.485798701848d-4
  real(kind=dp),parameter,public  ::  au2k  = 3.1577464d5 
  real(kind=dp),parameter,public  ::  au2ang= 5.291772108d-1
  real(kind=dp),parameter,public  ::             &  
    sqrt3 = dsqrt(3.0d0)                        ,&
    sqrt5 = dsqrt(5.0d0)                        ,&
    sqrt7 = dsqrt(7.0d0)                        ,&
    pi    = 3.141592653589793238462643383279_dp ,&  
    hbar  = 1.05457180013d-34   
    !au2amu 原子单位质量(me)转换为相对原子质量(C12/12)
    !the atomic mass of a carbon-12 atom is about 1.998467052 × 10−26 kg
    !one twelfth of the mass of an unbound neutral atom of carbon-12 
    !in its nuclear and electronic ground state and at rest has a value 
    !of 1.660539040(20)×10−27 kg,
    !electron mass 9.10938291(40)×10−31 kg
  real(kind=dp),public :: kb
  real(kind=dp),parameter,public  ::  twopi = 2*pi
  
  complex(kind=dpc),parameter,public   :: eye=(0.0d0,1.0d0)
  complex(kind=dpc), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
  !! i as a complex variable
  complex(kind=dpc), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  !! 0 as a complex variable
  complex(kind=dpc), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)
  !! 1 as a complex variable
  
  !~~ NUMERICAL CONVERGENCE CONSTANTS ~~!
  real(kind=dp), parameter, public    :: eps2  = 1.0e-2_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps5  = 1.0e-5_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps6  = 1.0e-6_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps7  = 1.0e-7_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps8  = 1.0e-8_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp
   !! numerical convergence constant

  !~~ PHYSICAL CONSTANTS ~~!
  !
  ! Values of the fundamental constants taken from 
  ! http://physics.nist.gov/cuu/Constants/index.html
  ! ##### CODATA 2010 ##### !
  ! #warning "SCSH INFO: Using CODATA 2010 constant values"
  real(kind=dp), parameter, public :: Avo_con_SI=6.022140857e23_dp        !mol-1
  real(kind=dp), parameter, public :: elem_charge_SI=1.602176565e-19_dp   ! C
  real(kind=dp), parameter, public :: sqrt_elem_charge_SI=elem_charge_SI**2
  !! elemental charge   ->  e
  real(kind=dp), parameter, public :: elec_mass_SI=9.10938291e-31_dp      ! kg
  real(kind=dp), parameter, public :: amu_mass_SI =1.660539040e-27_dp     ! kg
  !! electron mass      ->  $$m_e$$
  real(kind=dp), parameter, public :: hbar_SI=1.054571726e-34_dp          ! J*s
  !! hbar               ->  $$\hbar$$
  real(kind=dp), parameter, public :: k_B_SI=1.3806488e-23_dp             ! J/K
  !! Boltzman Constant  ->  $$k_B$$ 
  real(kind=dp), parameter, public :: bohr_magn_SI=927.400968e-26_dp      ! J/T
  !! Bohr magneton      ->  $$\mu_B$$
  real(kind=dp), parameter, public :: eps0_SI=8.854187817e-12_dp          !F / m
  !! Vacuum Dielectric Constant ->  $$\epsilon_0$$
  real(kind=dp), parameter, public :: speedlight_SI=299792458.0_dp        ! m / s
  !! Speed of light     ->  $$c$$
  real(kind=dp), parameter, public :: eV_au=3.674932379e-2_dp             ! (see table of Conv. Factors)  
  !! Electron Volt in atomic units
  real(kind=dp), parameter, public :: bohr_angstrom_internal=0.52917721092_dp
  !! Bohr to Anstrom Conversion factor
  real(kind=dp), parameter, public :: bohr = bohr_angstrom_internal
  !! 4*pi*eps0
  real(kind=dp), parameter, public :: fopieps0 = 4*pi*eps0_SI
  !real(kind=dp), parameter, public :: THz2womiga = 1.0e12_dp                     !f=  rad/s
  
  ! Leave the length to this value, and don't exceed in length (needed for output formatting)
  character(len=75), parameter, public :: constants_version_str1 = "-> Using CODATA 2010 constant values"
  character(len=75), parameter, public :: constants_version_str2 = "   (http://physics.nist.gov/cuu/Constants/index.html)"

end module constants
  