   !***********************************************************************
   !
   !  Module stochastic namelist records
   !
   !> \brief   module accesses the MPAS stochastic physics namelist 
   !> \author  Ning Wang
   !> \date    Oct 2024
   !
   !-----------------------------------------------------------------------
module stochy_nml_rec
   use kinddef, only: kind_phys, RKIND, StrKIND
   implicit none

   contains

   !***********************************************************************
   !
   !  routine get_nml_rec
   !
   !> \brief   routine to retrieve stochastic physics namelist records 
   !> \author  Ning Wang
   !> \date    Oct 2024
   !> \details 
   !>  This routine retrieves stochastic physics namelist records which
   !>  are defined in the registry files. 
   !
   !-----------------------------------------------------------------------
   subroutine get_nml_rec (domain,me,deltim,iret)
      use stochy_namelist_def
      use mpas_pool_routines

      implicit none

      type(domain_type),    intent(inout):: domain
      integer,              intent(out) :: iret
      integer,              intent(in)  :: me
      real,                 intent(in)  :: deltim

      real l_min
      real :: r_earth,circ,tmp_lat,tol
      integer k,ios
      integer,parameter :: four=4

      type(mpas_pool_type) :: configPool

! Local variables in this subroutine that will be filled in with values from
! the MPAS namelist using the mpas_pool_get_config() subroutine (which is a
! generic function) must have the same types and kinds or lengths as the "value"
! arguments for the type-specific versions of that function.  For example, the
! subroutine mpas_pool_get_config_real() that reads in a real value, assumes
! its "value" argument is a "real (kind=RKIND)" (see file mpas_pool_routines.F
! in the MPAS-Model code base).  Thus, here, any local real variables passed
! to mpas_pool_get_config() must also be delcared as "real (kind=RKIND)" because
! otherwise, the build of MPAS-Model with stochastic_physics will fail.  Similarly,
! any local character variables must be declared as "character (len=StrKIND)".
! Logical and integer variables can be declared simply as "logical" and "integer"
! because that's how the subroutines mpas_pool_get_config_logical() and
! mpas_pool_get_config_int() in mpas_pool_routines.F declare their "value"
! arguments.
      real (kind=RKIND), pointer :: config_sppt_1
      real (kind=RKIND), pointer :: config_sppt_2
      real (kind=RKIND), pointer :: config_sppt_3
      real (kind=RKIND), pointer :: config_sppt_tau_1
      real (kind=RKIND), pointer :: config_sppt_tau_2
      real (kind=RKIND), pointer :: config_sppt_tau_3
      real (kind=RKIND), pointer :: config_sppt_lscale_1
      real (kind=RKIND), pointer :: config_sppt_lscale_2
      real (kind=RKIND), pointer :: config_sppt_lscale_3
      real (kind=RKIND), pointer :: config_sppt_hgt_top1
      real (kind=RKIND), pointer :: config_sppt_hgt_top2
      logical, pointer :: config_do_sppt
      logical, pointer :: config_sppt_logit
      logical, pointer :: config_sppt_sfclimit
!      integer, pointer :: config_iseed_sppt1, config_iseed_sppt2, config_iseed_sppt3
      character (len=StrKIND), pointer :: config_iseed_sppt1, config_iseed_sppt2, config_iseed_sppt3
      integer io_stat
      integer, pointer :: config_spptint

      logical, pointer :: config_do_skeb
      logical, pointer :: config_stochini

!     spectral resolution defintion
      ntrunc=-999
      lon_s=-999
      lat_s=-999
      sppt             = -999.  ! stochastic physics tendency amplitude
      iseed_sppt       = 0      ! random seeds (if 0 use system clock)
! logicals
      do_sppt = .false.
      use_zmtnblck = .false.
      new_lscale = .true.
! parameters to control vertical tapering of stochastic physics with
! height
      spptint      = 0
      sppt_hgt_top1 = 15000.0
      sppt_hgt_top2 = 27000.0
! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_sfclimit = .false.
      pbl_taper = (/0.0,0.5,1.0,1.0,1.0,1.0,1.0/)

      sppt_logit        = .false. ! logit transform for sppt to bounded interval [-1,+1]
      stochini          = .false. ! true= read in pattern, false=initialize from seed

! retrieve namelist rec
      configPool = domain % blocklist % configs
      call mpas_pool_get_config(configPool, 'do_sppt', config_do_sppt)
      call mpas_pool_get_config(configPool, 'config_spptint', config_spptint)
      call mpas_pool_get_config(configPool, 'config_sppt_1', config_sppt_1)
      call mpas_pool_get_config(configPool, 'config_sppt_2', config_sppt_2)
      call mpas_pool_get_config(configPool, 'config_sppt_3', config_sppt_3)
      call mpas_pool_get_config(configPool, 'config_sppt_tau_1', config_sppt_tau_1)
      call mpas_pool_get_config(configPool, 'config_sppt_tau_2', config_sppt_tau_2)
      call mpas_pool_get_config(configPool, 'config_sppt_tau_3', config_sppt_tau_3)
      call mpas_pool_get_config(configPool, 'config_sppt_lscale_1', config_sppt_lscale_1)
      call mpas_pool_get_config(configPool, 'config_sppt_lscale_2', config_sppt_lscale_2)
      call mpas_pool_get_config(configPool, 'config_sppt_lscale_3', config_sppt_lscale_3)
      call mpas_pool_get_config(configPool, 'config_sppt_logit', config_sppt_logit)
      call mpas_pool_get_config(configPool, 'config_sppt_sfclimit', config_sppt_sfclimit)
      call mpas_pool_get_config(configPool, 'config_iseed_sppt1', config_iseed_sppt1)
      call mpas_pool_get_config(configPool, 'config_iseed_sppt2', config_iseed_sppt2)
      call mpas_pool_get_config(configPool, 'config_iseed_sppt3', config_iseed_sppt3)
      call mpas_pool_get_config(configPool, 'config_sppt_hgt_top1', config_sppt_hgt_top1)
      call mpas_pool_get_config(configPool, 'config_sppt_hgt_top2', config_sppt_hgt_top2)

      call mpas_pool_get_config(configPool, 'do_skeb', config_do_skeb)
      call mpas_pool_get_config(configPool, 'config_stochini', config_stochini)

! Assign values read in from the MPAS-Model namelist to variables native to
! the stochastic_physics code.
      do_sppt = config_do_sppt

! spptint is declared as real/kind_dbl_prec, so convert to that type/kind.
! Note that config_spptint is declared as an integer.
      spptint = real(config_spptint, kind=kind_dbl_prec)

! sppt, sppt_tau, and sppt_lscale are declared as real/kind_dbl_prec, so
! convert to that type/kind.  It is not clear why they aren't defined as
! kind=kind_phys.
      sppt(1) = real(config_sppt_1, kind=kind_dbl_prec)
      sppt(2) = real(config_sppt_2, kind=kind_dbl_prec)
      sppt(3) = real(config_sppt_3, kind=kind_dbl_prec)
      sppt_tau(1) = real(config_sppt_tau_1, kind=kind_dbl_prec)
      sppt_tau(2) = real(config_sppt_tau_2, kind=kind_dbl_prec)
      sppt_tau(3) = real(config_sppt_tau_3, kind=kind_dbl_prec)
      sppt_lscale(1) = real(config_sppt_lscale_1, kind=kind_dbl_prec)
      sppt_lscale(2) = real(config_sppt_lscale_2, kind=kind_dbl_prec)
      sppt_lscale(3) = real(config_sppt_lscale_3, kind=kind_dbl_prec)

      sppt_logit = config_sppt_logit
      sppt_sfclimit = config_sppt_sfclimit

      read(config_iseed_sppt1, *, iostat=io_stat) iseed_sppt(1)
      if (io_stat /= 0) print*, 'String to integer failed, sppt seed 1 IOSTAT value:', io_stat
      read(config_iseed_sppt2, *, iostat=io_stat) iseed_sppt(2)
      if (io_stat /= 0) print*, 'String to integer failed, sppt seed 2 IOSTAT value:', io_stat
      read(config_iseed_sppt3, *, iostat=io_stat) iseed_sppt(3)
      if (io_stat /= 0) print*, 'String to integer failed, sppt seed 3 IOSTAT value:', io_stat

!      iseed_sppt(1) = config_iseed_sppt1
!      iseed_sppt(2) = config_iseed_sppt2
!      iseed_sppt(3) = config_iseed_sppt3

! sppt_hgt_top1 and sppt_hgt_top2 are declared as real/kind_dbl_prec, so
! convert to that type/kind.  It is not clear why they aren't defined as
! kind=kind_phys.
      sppt_hgt_top1 = real(config_sppt_hgt_top1, kind=kind_dbl_prec)
      sppt_hgt_top2 = real(config_sppt_hgt_top2, kind=kind_dbl_prec)

      do_skeb = config_do_skeb
      stochini = config_stochini

      r_earth  =6.3712e+6      ! radius of earth (m)
      tol=0.01  ! tolerance for calculations
      if (sppt(1) > 0 ) then
        do_sppt=.true.
      endif

      if (do_sppt) then
          if (spptint == 0.) spptint=deltim
          nssppt=nint(spptint/deltim)                              ! spptint in seconds
          if(nssppt<=0 .or. abs(nssppt-spptint/deltim)>tol) then
#ifdef STOCH_PHYS_DIAG
             write(0,*) "SPPT interval is invalid",spptint
#endif
            iret=9
            return
          endif
      endif 

!calculate ntrunc if not supplied
      if (ntrunc .LT. 1) then  
#ifdef STOCH_PHYS_DIAG
        if (me==0) print*,'ntrunc not supplied, calculating'
#endif
        circ=2*3.1415928*r_earth ! start with lengthscale that is circumference of the earth
        l_min=circ
        do k=1,5
           if (sppt(k).GT.0) l_min=min(sppt_lscale(k),l_min)
        enddo
        ntrunc=circ/l_min
#ifdef STOCH_PHYS_DIAG
        if (me==0) print*,'ntrunc calculated from l_min',l_min,ntrunc
#endif
      endif
     ! ensure lat_s is a mutiple of 4 with a reminader of two
      ntrunc=INT((ntrunc+1)/four)*four+2
#ifdef STOCH_PHYS_DIAG
      if (me==0) print*,'NOTE ntrunc adjusted for even nlats',ntrunc
#endif

! set up gaussian grid for ntrunc if not already defined. 
      if (lon_s.LT.1 .OR. lat_s.LT.1) then
        lat_s=ntrunc*1.5+1
        lon_s=lat_s*2+4
! Grid needs to be larger since interpolation is bi-linear
        lat_s=lat_s*2
        lon_s=lon_s*2
#ifdef STOCH_PHYS_DIAG
        if (me==0) print*,'gaussian grid not set, defining here',lon_s,lat_s
#endif
      endif

      iret = 0

      return
    end subroutine get_nml_rec

end module stochy_nml_rec
