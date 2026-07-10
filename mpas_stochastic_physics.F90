module mpas_stochastic_physics
  use kinddef, only: kind_phys, RKIND
  use mpi_f08
  use mpas_pool_routines
  use mpas_log, only : mpas_log_write, mpas_log_info
  implicit none

!local pointers:
  type(mpas_pool_type),pointer::  configs,       &
                                  mesh,         &
                                  state,        &
                                  tend_physics, &
                                  atm_input,    &
                                  sfc_input

  ! For stochastic physics pattern generation
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sppt_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: shum_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebu_wts
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: skebv_wts
!  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sfc_wts
  real(kind=kind_phys), dimension(:,:,:,:), allocatable, save :: spp_wts

  logical, save :: is_initialized = .false.
  integer, save :: lsoil = -999

  !roughness length for land
!  real(kind=kind_phys), dimension(:,:),   allocatable, save :: zorll
!  integer, dimension(:,:),   allocatable, save :: stype

  logical, pointer :: do_sppt_config, do_skeb_config, do_shum_config, do_spp_config 
  logical :: do_sppt = .true.
  logical :: do_skeb, do_shum, do_spp 

  real(kind=RKIND), dimension(:),pointer:: lonCell, latCell
  real(kind=RKIND), dimension(:,:),pointer:: stoch_pat_sppt  ! stoch(nVertLevels, nCells)
  real(kind=RKIND), dimension(:,:),pointer:: zgrid  ! zgrid(nVertLevelsP1, nCells)
  real(kind=RKIND), dimension(:,:),pointer:: stoch_pat_gg  ! stoch(lon_for, lat_leg)
  integer, dimension(:), pointer :: blksz
  integer, pointer :: nCells, nEdges, nVertLevels, nVertLevelsP1
  integer, pointer :: lon_f, lat_g
  real(kind=RKIND), pointer :: dt, sppt_1, sppt_2, sppt_3
  real(kind=RKIND), pointer :: sppt_tau_1, sppt_tau_2, sppt_tau_3
  real(kind=RKIND), pointer :: sppt_lscale_1, sppt_lscale_2, sppt_lscale_3
  integer, pointer :: spptint

  integer, parameter :: tend_name_len = 32
  character(len=tend_name_len), allocatable, dimension(:) :: &
    tend_names_phys, tend_names_prog, tend_names_sppt
  character(len=tend_name_len), allocatable, target, dimension(:) :: &
    tends_to_perturb_sppt_phys, tends_to_perturb_sppt_prog
  logical :: sppt_setup_complete = .false.

  integer, save :: nblks

!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_pattern_init
  public stochastic_physics_pattern_adv
  public stochastic_physics_pattern_apply
  public dosppt

  contains

   !***********************************************************************
   !
   !  routine stochastic_physics_pattern_init
   !
   !> \brief   stochastic physics setup routine
   !> \author  Ning Wang
   !> \date    Oct 2024
   !> \details 
   !>  This routine is intended to setup the initial state for perturbation
   !>  pattern generation.  
   !
   !-----------------------------------------------------------------------
  subroutine stochastic_physics_pattern_init (domain)

    use stochastic_physics_m,  only: init_stochastic_physics
    implicit none

    type(domain_type),intent(inout):: domain


!local variables:
    type(mpas_pool_type) :: configPool
    type(block_type),pointer:: block
    integer :: maxblk, blk_sz(16)
    type(MPI_comm) ::  mpi_comm
    integer :: mpi_root, pid
    integer :: i, iret
    real(kind=kind_phys):: sppt_amp, dtp
    real(kind=kind_phys), dimension(:), pointer:: zk
    real(kind=kind_phys), dimension(:,:), pointer:: xlon, xlat
    character(len=80) :: msg
    character(len=:), allocatable :: tend_names_str

    pid = domain%dminfo%my_proc_id  
    configPool = domain%blocklist%configs

    call mpas_log_write(' ')
    call mpas_log_write('--- enter stochastic_physics_pattern_init( )')

    call mpas_pool_get_config(configPool, 'do_sppt', do_sppt_config)

    if (do_sppt_config .eqv. .false.) then
      call mpas_log_write('    do_sppt is false, do not perturb physics tendencies.')
      return
    endif

    call mpas_pool_get_config(configPool, 'config_dt', dt)
    call mpas_pool_get_config(configPool, 'config_spptint', spptint)

    call mpas_log_write('    config_dt = $r', realArgs=(/dt/))
    call mpas_log_write('    config_spptint = $i', intArgs=(/spptint/))
    call mpas_log_write('    number of processors = $i', intArgs=(/domain%dminfo%nprocs/))

    call mpas_pool_get_config(configPool, 'config_sppt_1', sppt_1)
    call mpas_pool_get_config(configPool, 'config_sppt_2', sppt_2)
    call mpas_pool_get_config(configPool, 'config_sppt_3', sppt_3)
    call mpas_pool_get_config(configPool, 'config_sppt_tau_1', sppt_tau_1)
    call mpas_pool_get_config(configPool, 'config_sppt_tau_2', sppt_tau_2)
    call mpas_pool_get_config(configPool, 'config_sppt_tau_3', sppt_tau_3)
    call mpas_pool_get_config(configPool, 'config_sppt_lscale_1', sppt_lscale_1)
    call mpas_pool_get_config(configPool, 'config_sppt_lscale_2', sppt_lscale_2)
    call mpas_pool_get_config(configPool, 'config_sppt_lscale_3', sppt_lscale_3)

    call mpas_log_write('    config_sppt_1 = $r', realArgs=(/sppt_1/))
    call mpas_log_write('    config_sppt_2 = $r', realArgs=(/sppt_2/))
    call mpas_log_write('    config_sppt_3 = $r', realArgs=(/sppt_3/))
    call mpas_log_write('    config_sppt_tau_1 = $r', realArgs=(/sppt_tau_1/))
    call mpas_log_write('    config_sppt_tau_2 = $r', realArgs=(/sppt_tau_2/))
    call mpas_log_write('    config_sppt_tau_3 = $r', realArgs=(/sppt_tau_3/))
    call mpas_log_write('    config_sppt_lscale_1 = $r', realArgs=(/sppt_lscale_1/))
    call mpas_log_write('    config_sppt_lscale_2 = $r', realArgs=(/sppt_lscale_2/))
    call mpas_log_write('    config_sppt_lscale_3 = $r', realArgs=(/sppt_lscale_3/))
    call mpas_log_write(' ')

    nblks = 0; maxblk = 1
    blk_sz = 0

    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      call mpas_pool_get_subpool(block%structs,'mesh' ,mesh)
      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      call mpas_pool_get_dimension(mesh,'nEdges',nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(mesh, 'lon_for', lon_f)
      call mpas_pool_get_dimension(mesh, 'lat_leg', lat_g)

      call mpas_pool_get_array(mesh,'zgrid',zgrid)

      blk_sz(nblks) = nCells
      if (nCells > maxblk) then
        maxblk = nCells
      endif 
      call mpas_pool_get_subpool(block%structs,'tend_physics' ,tend_physics)
      call mpas_pool_get_array(tend_physics,'stoch_pattern_sppt',stoch_pat_sppt)
      call mpas_pool_get_array(tend_physics,'stoch_pattern_gg',stoch_pat_gg)
      stoch_pat_sppt(:,:) = real(domain%dminfo%my_proc_id, kind=RKIND)
      block => block%next
    enddo

    allocate(blksz(nblks))
    blksz(1:nblks) = blk_sz(1:nblks)

    allocate(lonCell(maxblk))
    allocate(latCell(maxblk))
    allocate(xlon(nblks,maxblk))
    allocate(xlat(nblks,maxblk))
    allocate(zk(nVertLevelsP1))
!    allocate(ak(nVertLevels), bk(nVertLevels))

    call get_zk(real(zgrid,kind=kind_phys), zk)

    nblks = 0
    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      call mpas_pool_get_subpool(block%structs,'mesh' ,mesh)
      call mpas_pool_get_array(mesh,'lonCell',lonCell)
      call mpas_pool_get_array(mesh,'latCell',latCell)
      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      xlon(nblks,1:nCells) = real(lonCell(1:nCells), kind=kind_phys)
      xlat(nblks,1:nCells) = real(latCell(1:nCells), kind=kind_phys)
      block => block%next
    enddo

!    if (pid == 0) print*, 'Call init_stochastic_physics (domain, ...).'
    mpi_root = 0; mpi_comm = domain%dminfo%comm
    dtp = real(dt, kind=kind_phys)
    call init_stochastic_physics(domain, nVertLevels, blksz, dtp, sppt_amp, &
         xlon, xlat, zk, mpi_comm, mpi_root, iret) 

!    if (pid == 0) print*, 'Done init_stochastic_physics (domain, ...).', iret

    ! If the initialization fails, print out an error message and exit.
    if (iret /= 0) then
!      do_sppt = .false.
!      do_skeb = .false.
!      do_shum = .false.
!      do_spp = .false.
      write(msg, '("Pattern initialization failed:  iret = ", I0, A, "Stopping.")') &
           iret, new_line('a')
      call mpas_log_write(msg)
      stop trim(msg)
    endif
!
! The subroutine stochastic_physics_pattern_apply() can apply perturbations
! (from any given stochastic physics scheme, e.g. SPPT, SKEB, etc) to one
! of two categories of tendencies:
!
! 1) Process level physics tendencies, e.g. wind component tendencies from
!    convection.
! 2) Accumulated physics tendencies that appear in (on the right-hand side
!    of) the prognostic equation for a state variable, e.g. the sum of the
!    physics tendencies in the thermodynamic evolution equation that MPAS
!    solves (named "tend_rtheta_physics").
!
! Specify the names of these two categories of tendencies.
!
    allocate(tend_names_phys(4))
    tend_names_phys(1) = "rucuten"
    tend_names_phys(2) = "rvcuten"
    tend_names_phys(3) = "rublten"
    tend_names_phys(4) = "rvblten"
    allocate(tend_names_prog(2))
    tend_names_prog(1) = "tend_rtheta_physics"
    !tend_names_prog(2) = "tend_rho_physics"
    tend_names_prog(2) = "tend_ru_physics"
!
! Print out the list of physics tendencies that stochastic schemes may perturb.
!
    call char_array_to_str(tend_names_phys, tend_names_str, separator=', ')
    call mpas_log_write( &
         'The tendencies in the "phys" (physics) category that may be perturbed by ' // &
         'stochastic physics schemes are:')
    call mpas_log_write('    ' // trim(tend_names_str))
!
! Print out the list of prognostic tendencies that stochastic schemes may perturb.
!
    call char_array_to_str(tend_names_prog, tend_names_str, separator=', ')
    call mpas_log_write( &
         'The tendencies in the "prog" (prognostic) category that may be perturbed by ' // &
         'stochastic physics schemes are:')
    call mpas_log_write('    ' // trim(tend_names_str))

!    print*, 'Exit stochastic_physics_pattern_init', pid

  end subroutine stochastic_physics_pattern_init


   !***********************************************************************
   !
   !  routine stochastic_physics_pattern_adv
   !
   !> \brief   stochastic physics pattern advance routine 
   !> \author  Ning Wang
   !> \date    Oct 2024
   !> \details 
   !>  This routine advances the perturbation pattern.  
   !
   !-----------------------------------------------------------------------
  subroutine stochastic_physics_pattern_adv (domain, ts_ct)

    use stochastic_physics_m,  only: run_stochastic_physics
    implicit none

    type(domain_type),intent(inout):: domain
    integer, intent(in) :: ts_ct

    real(kind=kind_phys), dimension(:,:,:),pointer:: st_pat 
    integer :: iret

    type(block_type),pointer:: block
    integer :: k, kdt, i, j
    integer :: pid
    logical do_pattern
    character(len=80) :: msg

    pid = domain%dminfo%my_proc_id  

!debug code ...
!    do_pattern = .false.
!    if (do_pattern .eqv. .false.) then
!      stoch_pat_sppt = 1.0
!      return
!    endif
!
!    if (do_sppt .eqv. .false.) return

    call mpas_log_write('--- enter stochastic_physics_pattern_adv( )')

    kdt = int(dt)
    nblks = 0
    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      block => block%next
    enddo
    
    if (nblks > 1)  then 
      pid = domain%dminfo%my_proc_id  
      call mpas_log_write('--- There are more than one block in this MPI rank! pid = $i', intArgs=(/pid/), messageType=MPAS_LOG_WARN)
    endif

    allocate(st_pat(1, blksz(1), nVertLevels)) 
    call run_stochastic_physics(nVertLevels, ts_ct, blksz, st_pat, iret)
!    call run_stochastic_physics(nVertLevels, ts_ct, blksz, st_pat, stoch_pat_gg, iret)
    if (iret == 0) then ! if pattern is advanced (updated)
      do k = 1, nVertLevels
        stoch_pat_sppt(k,1:nCells) = real(st_pat(1,1:nCells,k), kind=RKIND)
      enddo
    else
      write(msg, '("Pattern advancement failed:  iret = ", I0, A, "Stopping.")') &
           iret, new_line('a')
      call mpas_log_write(msg)
      stop trim(msg)
    endif
    deallocate(st_pat)

  end subroutine stochastic_physics_pattern_adv

   !***********************************************************************
   !
   !  routine stochastic_physics_pattern_apply
   !
   !> \brief   stochastic physics application routine
   !> \author  Ning Wang
   !> \date    Oct 2024
   !> \details 
   !>  This routine applies a perturbation pattern to the specified fields.  
   !>  It is called to perturb the cell-centered tendencies before they are
   !>  used in the shallow-water equations.
   !
   !>  For velocity components, the pattern is applied to the tendencies 
   !>  from both PBL and convection schemes, at the cell centers.
   !
   !-----------------------------------------------------------------------
  subroutine stochastic_physics_pattern_apply(domain, tend_category)
    type(domain_type),intent(in):: domain
    character(len=4), intent(in) :: tend_category

    real(kind=RKIND), dimension(:,:), pointer:: tend_array 
    character(len=tend_name_len), dimension(:), pointer :: tends_to_perturb
    character(len=:), allocatable :: tends_to_perturb_str
    type(block_type), pointer:: block
    integer :: i, nblks, pid
    character(len=tend_name_len) :: tend_name
    character(len=:), allocatable :: stoch_scheme_name, tend_category_long

    pid = domain%dminfo%my_proc_id

#ifdef STOCH_PHYS_DIAG
    if (pid == 0) print*, 'Enter stochastic_physics_pattern_apply, pid = ', pid
#endif

    nblks = 0
    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      block => block%next
    enddo
    
    if (nblks > 1)  then 
      call mpas_log_write('--- There are more than one block in this MPI rank! pid = $i', intArgs=(/pid/), messageType=MPAS_LOG_WARN)
    endif
!
! If SPPT is enabled, retrieve the specified tendencies that need to be
! perturbed and apply the patterns (at all levels) to the tendency.
!
    if (do_sppt) then

      stoch_scheme_name = "SPPT"
!
! If SPPT has not yet been set up, do so.  This should happen only once
! per forecast.
!
      if (.not. sppt_setup_complete) then
        call setup_sppt_mpas()
        sppt_setup_complete = .true.
      end if

      if (trim(tend_category) == 'phys') then
        tend_category_long = "physics"
        tends_to_perturb => tends_to_perturb_sppt_phys
      else if (trim(tend_category) == 'prog') then
        tend_category_long = "prognostic"
        tends_to_perturb => tends_to_perturb_sppt_prog
      else
         stop "Invalid tend_category: tend_category = " // trim(tend_category)
      end if

      call char_array_to_str(tends_to_perturb, tends_to_perturb_str, separator=', ')
      call mpas_log_write( &
           'Applying ' // trim(stoch_scheme_name) // ' perturbations to the following ' // &
           trim(tend_category_long) // ' tendencies:  ' // tends_to_perturb_str)
      do i = 1, size(tends_to_perturb)
        tend_name = tends_to_perturb(i)
        call mpas_pool_get_array(tend_physics, trim(tend_name), tend_array)
        call apply_pattern(tend_array, stoch_pat_sppt)
      enddo

    endif

  end subroutine stochastic_physics_pattern_apply

  subroutine apply_pattern(tendency, stoch_pattern)
    real(kind=RKIND), dimension(nVertLevels,nCells) :: tendency
    real(kind=RKIND), dimension(nVertLevels,nCells) :: stoch_pattern

    tendency = tendency*stoch_pattern

  end subroutine apply_pattern

  !***********************************************************************
  !
  !  routine setup_sppt_mpas
  !
  !> \brief   Subroutine to set up SPPT to run with MPAS.
  !> \author  Gerard Ketefian
  !> \date    July 2026
  !> \details 
  !>  This routine is intended to setup variables needed by SPPT to run
  !>  with the MPAS dycore.
  !
  !-----------------------------------------------------------------------
  subroutine setup_sppt_mpas()

    character(len=tend_name_len), allocatable, dimension(:) :: tend_names_all_categ
    character(len=tend_name_len) :: tend_name
    character(len=:), allocatable :: stoch_scheme_name, tend_names_str
    integer :: i, indx_phys, indx_prog
    real(kind=RKIND), dimension(:,:), pointer:: tend_array 
    logical :: sppt_tend_is_valid
    logical, allocatable, dimension(:) :: keep_mask_phys, keep_mask_prog
    character(len=:), allocatable :: &
      tend_names_sppt_str, tends_to_perturb_sppt_phys_str, tends_to_perturb_sppt_prog_str

    stoch_scheme_name = "SPPT"
    call mpas_log_write('')
    call mpas_log_write( &
         'Setting up ' // trim(stoch_scheme_name) // '-associated variables...')
!
! Set all tendencies (both process-level and accumulated state tendencies) that
! SPPT may (but not necessarily will) perturb.
!
    allocate(tend_names_sppt(5))
    tend_names_sppt(1) = "rucuten"
    tend_names_sppt(2) = "rvcuten"
    tend_names_sppt(3) = "rublten"
    tend_names_sppt(4) = "rvblten"
    tend_names_sppt(5) = "tend_rtheta_physics"
!
! Print out the list of tendencies that SPPT may perturb.
!
    call char_array_to_str(tend_names_sppt, tend_names_sppt_str, separator=', ')
    call mpas_log_write( &
         'The tendencies that ' // trim(stoch_scheme_name) // &
         ' may (but not necessarily will) perturb are:')
    call mpas_log_write('    ' // trim(tend_names_sppt_str))
!
! Form list of all tendencies (i.e. either in the physics or the prognostic 
! categories) that may be perturbed by stochastic schemes.  This is useful
! in the various consistency checks below.
! 
! Note that this is a whole-array assignment, i.e. tend_names_all_categ is
! auto-(re)allocated to size(tend_names_phys) + size(tend_names_prog)
!
    tend_names_all_categ = [tend_names_phys, tend_names_prog]
!
! Form lists of the tendencies in each category (physics or progonostic)
! that will (not may) be perturbed by SPPT.  Use masking arrays to select
! these from the general list of tendencies above that may be perturbed.
!
    allocate(keep_mask_phys(size(tend_names_sppt)))
    keep_mask_phys(:) = .false.

    allocate(keep_mask_prog(size(tend_names_sppt)))
    keep_mask_prog(:) = .false.

    call mpas_log_write( &
         'Selecting from the physics and prognostic tendency lists those tendencies that ' // &
         'will be perturbed by ' // trim(stoch_scheme_name) // '...')

    do i=1, size(tend_names_sppt)

      tend_name = tend_names_sppt(i)
      call mpas_log_write( &
         '  Considering tendency: ' // trim(tend_name))
!
! First, check that the current tendency is valid, i.e. that it exists in
! either the list of physics tendencies or prognostic tendencies that may
! be perturbed by the stochastic schemes.  If not, print out a message and
! stop.
!
      sppt_tend_is_valid = any(tend_name == tend_names_all_categ)
      if (.not. sppt_tend_is_valid) then
        call mpas_log_write( &
             'The current SPPT tendency is not valid because it is neither ' // &
             'in the list of physics nor prognostic tendencies that may be perturbed ' // &
             'by the stochastic physics schemes:')
        call mpas_log_write( &
             '  tend_name = "' // trim(tend_name) // '"')
        call char_array_to_str(tend_names_phys, tend_names_str, &
                               separator=", ", quote_elems=.true.)
        call mpas_log_write( &
             '  tend_names_phys = (/ ' // trim(tend_names_str) // ' /)')
        call char_array_to_str(tend_names_prog, tend_names_str, &
                               separator=", ", quote_elems=.true.)
        call mpas_log_write( &
             '  tend_names_prog = (/ ' // trim(tend_names_str) // ' /)')
        stop
      end if
!
! Now check whether the pointer to the tendency array is associated (i.e.
! doesn't point to a null address).  If so, set its mask value in the list
! of physics or prognostic tendency names to true.
!
      call mpas_pool_get_array(tend_physics, trim(tend_name), tend_array)
      if (associated(tend_array)) then

        call mpas_log_write( &
             '    Will apply ' // trim(stoch_scheme_name) // ' perturbations to ' // &
             trim(tend_name))

        if (any(tend_name == tend_names_phys)) then
          indx_phys = findloc(tend_names_phys, tend_name, dim=1)
          keep_mask_phys(indx_phys) = .true.
        else if (any(tend_name == tend_names_prog)) then
          indx_prog = findloc(tend_names_prog, tend_name, dim=1)
          keep_mask_prog(indx_prog) = .true.
        end if

      else
        call mpas_log_write( &
             '  Will NOT apply ' // trim(stoch_scheme_name) // ' perturbations to ' // &
             trim(tend_name) // ' because the latter is not associated.')
      endif

    end do
!
! Set and then print out the arrays containing the lists of physics and
! prognostic tendencies, respectively, that SPPT will perturb.
!
    tends_to_perturb_sppt_phys = pack(tend_names_phys, keep_mask_phys)
    call char_array_to_str( &
         tends_to_perturb_sppt_phys, tends_to_perturb_sppt_phys_str, separator=', ')
    call mpas_log_write( &
         'Tendencies in the physics category to be perturbed by ' // &
         trim(stoch_scheme_name) // ' are: ')
    call mpas_log_write('    ' // trim(tends_to_perturb_sppt_phys_str))

    tends_to_perturb_sppt_prog = pack(tend_names_prog, keep_mask_prog)
    call char_array_to_str( &
         tends_to_perturb_sppt_prog, tends_to_perturb_sppt_prog_str, separator=', ')
    call mpas_log_write( &
         'Tendencies in the prognostic category to be perturbed by ' // &
         trim(stoch_scheme_name) // ' are: ')
    call mpas_log_write('    ' // trim(tends_to_perturb_sppt_prog_str))

    call mpas_log_write( &
         'Done setting up ' // trim(stoch_scheme_name) // '-associated variables.')

  end subroutine setup_sppt_mpas

   subroutine get_zk(zgrid, zk)
    real(kind=kind_phys), dimension(nVertLevelsP1,nCells) :: zgrid
    real(kind=kind_phys), dimension(nVertlevelsP1) :: zk
    
    real(kind=kind_phys) :: z_btm
    integer :: i
    
    z_btm = 10000.0 
    do i = 1, nCells
      if (zgrid(1,i) == 0) then 
        zk(:) = zgrid(:,i)
        exit
      else
        if (z_btm > zgrid(1,i)) then
          zk(:) = zgrid(:,i)
          z_btm = zgrid(1,i)
        endif
      endif
    enddo

   end subroutine get_zk

   logical function dosppt(domain)
     type(domain_type),intent(in):: domain

     dosppt = .false.
     call mpas_pool_get_config(domain % blocklist % configs, 'do_sppt', do_sppt_config)
     if (do_sppt_config .eqv. .true. .and. do_sppt .eqv. .true.) then
       dosppt = .true. 
     endif 

   end function dosppt

   subroutine char_array_to_str(char_array, str, separator, quote_elems)
     character(len=*), dimension(:), intent(in) :: char_array
     character(len=:), allocatable, intent(out) :: str
     character(len=*), intent(in), optional :: separator
     logical, intent(in), optional :: quote_elems

     character(len=:), allocatable :: quote, sep
     integer :: i

     if (present(separator)) then
       sep = separator
     else
       sep = ' '
     end if

     quote = ''
     if (present(quote_elems)) then
       if (quote_elems) quote = '"'
     end if

     str = ''
     do i=1, size(char_array)
       if (i == 1) then
         str = quote // trim(char_array(i)) // quote
       else
         str = str // sep // quote // trim(char_array(i)) // quote
       end if
     end do

   end subroutine char_array_to_str

end module mpas_stochastic_physics
