module kinddef
#ifdef MPAS
  use mpas_kind_types, only: R4KIND, R8KIND, RKIND, I8KIND, StrKIND
#endif
  implicit none

  private

  public :: kind_phys
  public :: kind_dbl_prec, kind_qdt_prec
  public :: kind_io8
#ifdef MPAS
  public :: RKIND, StrKIND
#endif

  ! FV3
#ifdef FV3
  integer, parameter :: kind_dbl_prec = 8
#ifdef CCPP_32BIT
  integer, parameter :: kind_phys     = 4
#else
  integer, parameter :: kind_phys     = 8
#endif
  integer, parameter :: kind_io8      = kind_dbl_prec
#endif

  ! MPAS
#ifdef MPAS
  integer, parameter :: kind_dbl_prec = R8KIND
#ifdef CCPP_32BIT
  integer, parameter :: kind_phys     = R4KIND
#else
  integer, parameter :: kind_phys     = R8KIND
#endif
  integer, parameter :: kind_io8      = I8KIND
#endif

#ifndef MPAS
#ifndef FV3
  integer, parameter :: kind_phys     = 8
  integer, parameter :: kind_dbl_prec = 8
  integer, parameter :: kind_io8      = kind_dbl_prec
#endif
#endif
  
#ifdef NO_QUAD_PRECISION
  integer, parameter :: kind_qdt_prec = 8
#else
  integer, parameter :: kind_qdt_prec = 16
#endif

end module kinddef
