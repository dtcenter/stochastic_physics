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
  character(len=tend_name_len), pointer, dimension(:) :: tend_names_phys
  character(len=tend_name_len), pointer, dimension(:) :: tend_names_prog
  character(len=tend_name_len), pointer, dimension(:) :: tend_names_sppt

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
    integer :: iret
    real(kind=kind_phys):: sppt_amp, dtp
    real(kind=kind_phys), dimension(:), pointer:: zk
    real(kind=kind_phys), dimension(:,:), pointer:: xlon, xlat
    character(len=80) :: msg

    pid = domain%dminfo%my_proc_id  
    configPool = domain%blocklist%configs

    call mpas_log_write(' ')
    call mpas_log_write('--- enter stochastic_physics_pattern_init( )')

    call mpas_pool_get_config(configPool, 'do_sppt', do_sppt_config)

    if (do_sppt_config .eqv. .false.) then
      call mpas_log_write('    do_sppt is false, do not perturb physics tendencies.')
      return
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
    ! Set all tendencies (process-level and accumulated state tendencies) that
    ! SPPT may perturb.
    !
    allocate(tend_names_sppt(5))
    tend_names_sppt(1) = "rucuten"
    tend_names_sppt(2) = "rvcuten"
    tend_names_sppt(3) = "rublten"
    tend_names_sppt(4) = "rvblten"
    tend_names_sppt(5) = "tend_rtheta_physics"

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

    real(kind=RKIND), dimension(:,:),pointer:: tend_array 
    character(len=tend_name_len), dimension(:), pointer :: tend_names
    type(block_type),pointer:: block
    integer :: i, j, k, nblks, pid, num_tends_perturbed
    character(len=tend_name_len) :: tend_name
    character(len=:), allocatable :: tend_category_desc, stoch_scheme_name, tend_list_str

    pid = domain%dminfo%my_proc_id

    if (trim(tend_category) == 'phys') then
       tend_names => tend_names_phys
       tend_category_desc = 'tendencies from individual physics schemes'
    else if (trim(tend_category) == 'prog') then
       tend_names => tend_names_prog
       tend_category_desc = 'accumulated physics tendencies in the MPAS prognostic equations for state variables'
    else
       stop "Invalid tend_category: tend_category ="//trim(tend_category)
    end if

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

! If SPPT is enabled, retrieve the specified tendencies that need to be perturbed
! and apply the patterns (at all levels) to the tendency.
   if (do_sppt) then

     stoch_scheme_name = "SPPT"

     tend_list_str = ''
     do i = 1, size(tend_names_sppt)
       if (i == 1) then
         tend_list_str = '"' // trim(tend_names_sppt(i)) // '"'
       else
         tend_list_str = tend_list_str // ', "' // trim(tend_names_sppt(i)) // '"'
       end if
     end do

     call mpas_log_write('') 
     call mpas_log_write( &
          'Attempting to apply ' // trim(stoch_scheme_name) // ' perturbations to ' // &
          'tendencies in the "' // trim(tend_category) // '" category, i.e.')
     call mpas_log_write( &
          'to ' // trim(tend_category_desc) // '...')
     call mpas_log_write( &
          'The tendencies that ' // trim(stoch_scheme_name) // ' may perturb are: ' // &
          tend_list_str)

     num_tends_perturbed = 0
     do i = 1, size(tend_names)
       tend_name = tend_names(i)
       call mpas_log_write('  tend_name = "' // trim(tend_name) // '":')
       ! If the current tendency is in the list of tendencies that SPPT may perturb
       ! (tend_names_sppt), try to perturb it.
       if (any(tend_name == tend_names_sppt)) then
         call mpas_pool_get_array(tend_physics, trim(tend_name), tend_array)
         if (associated(tend_array)) then
           call mpas_log_write( &
                '    Applying ' // trim(stoch_scheme_name) // ' perturbation to "' // &
                trim(tend_name) // '"...')
           call apply_pattern(tend_array, stoch_pat_sppt)
           num_tends_perturbed = num_tends_perturbed + 1
           call mpas_log_write( '    Done.')
         else
           call mpas_log_write( &
                '    Not applying ' // trim(stoch_scheme_name) // ' perturbation to "' // &
                trim(tend_name) // '" because it is not allocated.')
         endif
       else
         call mpas_log_write( &
              '    Not applying ' // trim(stoch_scheme_name) // ' perturbation to "' // &
              trim(tend_name) // '" because it is not in the list of tendencies ' // &
              'that ' // trim(stoch_scheme_name) // ' may perturb.')
       end if
     enddo

     call mpas_log_write( &
          'Done applying ' // trim(stoch_scheme_name) // ' perturbations to ' // &
          'tendencies in the "' // trim(tend_category) // '" category.  A total ' // &
          'of $i tendencies were perturbed.', &
          intArgs=(/num_tends_perturbed/))

   endif

   nullify(tend_names)

  end subroutine stochastic_physics_pattern_apply

  subroutine apply_pattern(tendency, stoch_pattern)
    real(kind=RKIND), dimension(nVertLevels,nCells) :: tendency
    real(kind=RKIND), dimension(nVertLevels,nCells) :: stoch_pattern

    tendency = tendency*stoch_pattern

   end subroutine apply_pattern

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

end module mpas_stochastic_physics
