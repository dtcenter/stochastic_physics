module stochastic_physics_mpas
  use kinddef, only: kind_phys
  use mpi_f08
  use mpas_pool_routines
  use mpi_wrapper, only: mype
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
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: sfc_wts
  real(kind=kind_phys), dimension(:,:,:,:), allocatable, save :: spp_wts

  logical, save :: is_initialized = .false.
  integer, save :: lsoil = -999
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: smc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: stc
  real(kind=kind_phys), dimension(:,:,:), allocatable, save :: slc
  !
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: vfrac
  !albedo
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: snoalb
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnsf
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: alnwf
  !emissivity
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: semis
  !roughness length for land
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: zorll

  !real(kind=kind_phys), dimension(:,:),   allocatable, save :: stype
  integer, dimension(:,:),   allocatable, save :: stype

  ! For cellular automata
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: sst
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lmsk
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: lake
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: condition
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca_deep_cpl, ca_turb_cpl, ca_shal_cpl
  real(kind=kind_phys), dimension(:,:),   allocatable, save :: ca1_cpl, ca2_cpl, ca3_cpl
  logical, pointer :: do_sppt 
  real(kind=RKIND), dimension(:),pointer:: lonCell, latCell
  real(kind=RKIND), dimension(:,:),pointer:: stoch_pat  ! stoch(nVertLevels, nCells)
  real(kind=RKIND), dimension(:,:),pointer:: stoch_pat_gg  ! stoch(log_for, lat_leg)
  integer, dimension(:), pointer :: blksz
  integer, pointer :: nCells, nVertLevels
  integer, pointer :: lon_f, lat_g
  real(kind=RKIND), pointer :: dt, sppt_1, sppt_2, sppt_3
  integer, save :: nblks

!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_pattern_init
  public stochastic_physics_pattern_adv

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
  subroutine stochastic_physics_pattern_init (domain, ierr)

    use stochastic_physics,  only: init_stochastic_physics
    implicit none

    type(domain_type),intent(inout):: domain
    integer, intent(out) :: ierr

!local variables:
    type(block_type),pointer:: block
    integer :: maxblk, blk_sz(16)
    type(MPI_comm) ::  mpi_comm
    integer :: mpi_root
    integer :: iret
    real(kind=kind_phys):: sppt_amp, dtp
    real(kind=kind_phys), dimension(:), pointer:: ak, bk
    real(kind=kind_phys), dimension(:,:), pointer:: xlon, xlat
    
    print*, 'Enter stochastic_physics_pattern_init( )...', mype, domain % dminfo % nprocs

    call mpas_pool_get_config(domain % blocklist % configs, 'do_sppt', do_sppt)

    if (do_sppt .eqv. .false.) then
      print*, 'do_sppt is false, do not perturb physics tendencies.'
      return
    endif

    call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)
    print*, 'config - dt', dt
    print*,'domain%dm_info%comm', domain%dminfo%comm 
    print*,'domain%dm_info%nprocs', domain%dminfo%nprocs  
    print*,'domain%dm_info%my_proc_id', domain%dminfo%my_proc_id  
    call mpas_pool_get_config(domain % blocklist % configs, 'config_sppt_1', sppt_1)
    call mpas_pool_get_config(domain % blocklist % configs, 'config_sppt_2', sppt_2)
    call mpas_pool_get_config(domain % blocklist % configs, 'config_sppt_3', sppt_3)
    print*, 'config_sppt_1:', sppt_1
    print*, 'config_sppt_2:', sppt_2
    print*, 'config_sppt_3:', sppt_3

    nblks = 0; maxblk = 1
    blk_sz = 0

    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      call mpas_pool_get_subpool(block%structs,'mesh' ,mesh)
      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'lon_for', lon_f)
      call mpas_pool_get_dimension(mesh, 'lat_leg', lat_g)
      blk_sz(nblks) = nCells
      if (nCells > maxblk) then
        maxblk = nCells
      endif 
      call mpas_pool_get_subpool(block%structs,'tend_physics' ,tend_physics)
      call mpas_pool_get_array(tend_physics,'stoch_pattern',stoch_pat)
      call mpas_pool_get_array(tend_physics,'stoch_pattern_gg',stoch_pat_gg)
      stoch_pat(:,:) = real(domain%dminfo%my_proc_id)
!      print*, 'size(stoch_pat,:) - dimensions:', size(stoch_pat,1), size(stoch_pat,2)
!      print*, 'In do while loop, nblks and nCells : ', nblks, nCells
!      print*, 'In do while loop, lonf and latg : ', lon_f, lat_g
      block => block%next
    enddo

    allocate(blksz(nblks))
    blksz(1:nblks) = blk_sz(1:nblks)

    allocate(lonCell(maxblk))
    allocate(latCell(maxblk))
    allocate(xlon(nblks,maxblk))
    allocate(xlat(nblks,maxblk))
    allocate(ak(nVertLevels), bk(nVertLevels))

    nblks = 0
    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      call mpas_pool_get_subpool(block%structs,'mesh' ,mesh)
      call mpas_pool_get_array(mesh,'lonCell',lonCell)
      call mpas_pool_get_array(mesh,'latCell',latCell)
      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      xlon(nblks,1:nCells) = lonCell(1:nCells) 
      xlat(nblks,1:nCells) = latCell(1:nCells) 
      block => block%next
    enddo

    print*, 'Call init_stochastic_physics (domain, ...).'
    mpi_root = 0; mpi_comm = domain%dminfo%comm
    dtp = dt
    ierr = 0
    call init_stochastic_physics(domain, nVertLevels, blksz, dtp, sppt_amp, &
         xlon, xlat, ak, bk, mpi_comm, mpi_root, iret) 

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
  subroutine stochastic_physics_pattern_adv (domain, ierr)

    use stochastic_physics,  only: run_stochastic_physics
    implicit none

    type(domain_type),intent(inout):: domain
    integer, intent(out) :: ierr

    real(kind=kind_phys), dimension(:,:,:),pointer:: st_pat 

    type(block_type),pointer:: block
    integer :: k, kdt, i, j, pid

    if (do_sppt .eqv. .false.) return

    pid = domain%dminfo%my_proc_id  
    if (pid == 0) print*, 'Enter stochastic_physics_pattern_adv...'
!    print*, 'Enter stochastic_physics_pattern_adv, proc_id = ', pid

    kdt = int(dt)
    nblks = 0
    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      block => block%next
    enddo
    
    if (nblks > 1)  then 
      print*,'There are more than one block, in this MPI rank!',  domain%dminfo%my_proc_id  
    endif

    allocate(st_pat(1, blksz(1), nVertLevels)) 
    call run_stochastic_physics(nVertLevels, kdt, blksz, st_pat, stoch_pat_gg)
    do k = 1, nVertLevels
      stoch_pat(k,1:nCells) = st_pat(1,1:nCells,k)
    enddo

!    pid = domain%dminfo%my_proc_id  
!    do i = 1, 20
!      do j = 1, 20
!        stoch_pat_gg(pid*20+i,pid*20+j) = real(pid)
!      enddo
!    enddo  

    deallocate(st_pat)
    ierr = 0

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
  subroutine stochastic_physics_pattern_apply(domain, num_tends, tend_names, ierr)
    type(domain_type),intent(in):: domain
    integer, intent(in) :: num_tends
    character(len=32), intent(in) :: tend_names(num_tends)
    integer, intent(out) :: ierr

    real(kind=RKIND), dimension(:,:),pointer:: tend_array 

    type(block_type),pointer:: block
    integer :: i, j, k, nblks, pid
    character(len=32) :: tend

    if (do_sppt .eqv. .false.) return

    pid = domain%dminfo%my_proc_id  
    if (pid == 0) print*, 'Enter stochastic_physics_pattern_apply...'
!    print*, 'Enter stochastic_physics_pattern_apply, proc_id = ', pid

    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      block => block%next
    enddo
    
    if (nblks > 1)  then 
      print*,'There are more than one block, in this MPI rank',  domain%dminfo%my_proc_id  
    endif

 ! retrieve the specified tendencies that needs to be perturbed and
 ! apply the patterns (at all levels) to the tendency      
   do i = 1, num_tends
     tend = tend_names(i)
!     print*, "tendency from physics parameterization to be perturbed:", trim(tend)
     select case (trim(tend))
       case ("rucuten")
         call mpas_pool_get_array(tend_physics,'rucuten',tend_array)
         call apply_pattern(tend_array, stoch_pat)
         if (pid == 0) print*, 'apply perturbation pattern to rucuten'
       case ("rvcuten")
         call mpas_pool_get_array(tend_physics,'rvcuten',tend_array)
         call apply_pattern(tend_array, stoch_pat)
         if (pid == 0) print*, 'apply perturbation pattern to rvcuten'
       case ("rublten")
         call mpas_pool_get_array(tend_physics,'rublten',tend_array)
         call apply_pattern(tend_array, stoch_pat)
         if (pid == 0) print*, 'apply perturbation pattern to rublten'
       case ("rvblten")
         call mpas_pool_get_array(tend_physics,'rvblten',tend_array)
         call apply_pattern(tend_array, stoch_pat)
         if (pid == 0) print*, 'apply perturbation pattern to rvblten'
       case ("tend_rtheta_physics")
         call mpas_pool_get_array(tend_physics,'tend_rtheta_physics',tend_array)
         call apply_pattern(tend_array, stoch_pat)
         if (pid == 0) print*, 'apply perturbation pattern to rtheta'
       case ("tend_rho_physics")
         call mpas_pool_get_array(tend_physics,'tend_rho_physics',tend_array)
         call apply_pattern(tend_array, stoch_pat)
         if (pid == 0) print*, 'apply perturbation pattern to rho'
       case default 
         print*, "The specified tendency is unknown and not perturbed." 
     end select
   enddo

   ierr = 0

  end subroutine stochastic_physics_pattern_apply

  subroutine apply_pattern(tendency, stoch_pattern)
    real(kind=RKIND), dimension(nVertLevels,nCells) :: tendency
    real(kind=RKIND), dimension(nVertLevels,nCells) :: stoch_pattern

    tendency = tendency*stoch_pattern

   end subroutine apply_pattern

end module stochastic_physics_mpas
