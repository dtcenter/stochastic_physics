   !***********************************************************************
   !
   !  Module stochastic namelist records
   !
   !> \brief   module accesses the MPAS stochastic physics namelist 
   !> \author  Ning Wang
   !> \date    Oct 2024
   !
   !-----------------------------------------------------------------------
module stoch_nml_rec

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
      integer, pointer :: config_iseed_sppt1, config_iseed_sppt2, config_iseed_sppt3
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

      do_sppt = config_do_sppt
      spptint = config_spptint
      sppt(1) = config_sppt_1
      sppt(2) = config_sppt_2
      sppt(3) = config_sppt_3

      sppt_tau(1) = config_sppt_tau_1
      sppt_tau(2) = config_sppt_tau_2
      sppt_tau(3) = config_sppt_tau_3

      sppt_lscale(1) = config_sppt_lscale_1
      sppt_lscale(2) = config_sppt_lscale_2
      sppt_lscale(3) = config_sppt_lscale_3

      sppt_logit = config_sppt_logit
      sppt_sfclimit = config_sppt_sfclimit
      iseed_sppt(1) = config_iseed_sppt1
      iseed_sppt(2) = config_iseed_sppt2
      iseed_sppt(3) = config_iseed_sppt3

      sppt_hgt_top1 = config_sppt_hgt_top1
      sppt_hgt_top2 = config_sppt_hgt_top2

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
             write(0,*) "SPPT interval is invalid",spptint
            iret=9
            return
          endif
      endif 

!calculate ntrunc if not supplied
      if (ntrunc .LT. 1) then  
        if (me==0) print*,'ntrunc not supplied, calculating'
        circ=2*3.1415928*r_earth ! start with lengthscale that is circumference of the earth
        l_min=circ
        do k=1,5
           if (sppt(k).GT.0) l_min=min(sppt_lscale(k),l_min)
        enddo
        ntrunc=circ/l_min
        if (me==0) print*,'ntrunc calculated from l_min',l_min,ntrunc
      endif
     ! ensure lat_s is a mutiple of 4 with a reminader of two
      ntrunc=INT((ntrunc+1)/four)*four+2
      if (me==0) print*,'NOTE ntrunc adjusted for even nlats',ntrunc

! set up gaussian grid for ntrunc if not already defined. 
      if (lon_s.LT.1 .OR. lat_s.LT.1) then
        lat_s=ntrunc*1.5+1
        lon_s=lat_s*2+4
! Grid needs to be larger since interpolation is bi-linear
        lat_s=lat_s*2
        lon_s=lon_s*2
        if (me==0) print*,'gaussian grid not set, defining here',lon_s,lat_s
      endif

      iret = 0

      return
    end subroutine get_nml_rec

end module stoch_nml_rec
