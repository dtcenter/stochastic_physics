module kinddef
#ifdef MPAS
  use mpas_kind_types
#endif
  
      implicit none

      private

      public :: kind_phys
      public :: kind_dbl_prec, kind_qdt_prec
      public :: kind_io8

      ! kind_phys must match CCPP Physics kind_phys
#ifdef CCPP_32BIT
      integer, parameter :: kind_phys     = 4
#else
      integer, parameter :: kind_phys     = 8
#endif

#ifdef FV3
      integer, parameter :: kind_dbl_prec = 8
#endif
#ifdef MPAS
      integer, parameter :: kind_dbl_prec => RKIND
#endif
      integer, parameter :: kind_io8      = kind_dbl_prec

#ifdef NO_QUAD_PRECISION
      integer, parameter :: kind_qdt_prec = 8
#else
      integer, parameter :: kind_qdt_prec = 16
#endif

end module kinddef
