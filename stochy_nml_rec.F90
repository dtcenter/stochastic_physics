module stoch_nml_rec

   implicit none

   contains

   subroutine get_nml_rec (domain, me,deltim,iret)
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
      logical, pointer :: config_sppt_logit
      logical, pointer :: config_sppt_sfclimit
      integer, pointer :: config_iseed_sppt
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
      new_lscale = .false.
! parameters to control vertical tapering of stochastic physics with
! height
      sppt_sigtop1 = 0.1
      sppt_sigtop2 = 0.025
! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_sfclimit = .false.
      pbl_taper = (/0.0,0.5,1.0,1.0,1.0,1.0,1.0/)

      sppt_logit        = .false. ! logit transform for sppt to bounded interval [-1,+1]
      stochini          = .false. ! true= read in pattern, false=initialize from seed

! retrieve namelist rec
      configPool = domain % blocklist % configs
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
      call mpas_pool_get_config(configPool, 'config_iseed_sppt', config_iseed_sppt)
      call mpas_pool_get_config(configPool, 'config_stochini', config_stochini)

      sppt(1) = config_sppt_1
      sppt(2) = config_sppt_2
      sppt(3) = config_sppt_3
      iseed_sppt = config_iseed_sppt

      r_earth  =6.3712e+6      ! radius of earth (m)
      tol=0.01  ! tolerance for calculations
      IF (sppt(1) > 0 ) THEN
        do_sppt=.true.
      ENDIF

      IF (do_sppt) THEN
          IF (spptint == 0.) spptint=deltim
          nssppt=nint(spptint/deltim)                              ! spptint in seconds
          IF(nssppt<=0 .or. abs(nssppt-spptint/deltim)>tol) THEN
             WRITE(0,*) "SPPT interval is invalid",spptint
            iret=9
            return
          ENDIF
      ENDIF 

!calculate ntrunc if not supplied
     if (ntrunc .LT. 1) then  
        if (me==0) print*,'ntrunc not supplied, calculating'
        circ=2*3.1415928*r_earth ! start with lengthscale that is circumference of the earth
        l_min=circ
        do k=1,5
           if (sppt(k).GT.0) l_min=min(sppt_lscale(k),l_min)
           if (shum(k).GT.0) l_min=min(shum_lscale(k),l_min)
           if (skeb(k).GT.0) l_min=min(skeb_lscale(k),l_min)
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
!
      if (me == 0) then
         print *, 'stochastic physics'
         print *, ' do_sppt : ', do_sppt
      endif
      iret = 0
!
      return
      end subroutine get_nml_rec

end module stoch_nml_rec
