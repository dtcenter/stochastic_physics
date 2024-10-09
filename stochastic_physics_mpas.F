module stochastic_physics_mpas
  use kinddef, only: kind_phys
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
  real(kind=RKIND), dimension(:),pointer:: latCell,lonCell
  real(kind=RKIND), dimension(:,:),pointer:: st_pat 


!----------------
! Public Entities
!----------------
! functions
  public stochastic_physics_patinit
  public stochastic_physics_patadv

  contains

subroutine stochastic_physics_patinit (domain, ierr)

    use stochastic_physics,  only: init_stochastic_physics
    implicit none

    type(domain_type),intent(inout):: domain
    integer, intent(out) :: ierr

!local variables:
    type(block_type),pointer:: block
    integer :: maxblk, nblks, n, v
    logical :: param_update_flag
    integer, pointer :: nCells, nVertLevels, nThreads
    integer :: blksz(16)
    real(kind=RKIND), pointer :: dt
    real(kind=kind_phys):: sppt_amp
    real(kind=kind_phys), allocatable :: xlon(:,:)
    real(kind=kind_phys), allocatable :: xlat(:,:)
    
    print*, 'Enter stochastic_physics_patinit...', mype, domain % dminfo % nprocs

    call mpas_pool_get_config(domain % blocklist % configs, 'config_dt', dt)
    print*, 'config - dt', dt
    print*,'domain%dm_info%comm', domain%dminfo%comm 
    print*,'domain%dm_info%nprocs', domain%dminfo%nprocs  
    print*,'domain%dm_info%my_proc_id', domain%dminfo%my_proc_id  
    nblks = 0; maxblk = 1
    blksz = 0
    block => domain % blocklist
    do while(associated(block))
      nblks = nblks + 1
      call mpas_pool_get_subpool(block%structs,'mesh' ,mesh)
      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      blksz(nblks) = nCells
      if (nCells > maxblk) then
        maxblk = nCells
      endif 
      call mpas_pool_get_subpool(block%structs,'tend_physics' ,tend_physics)
      call mpas_pool_get_array(tend_physics,'stoch_pattern',st_pat)
      st_pat(:,:) = real(domain%dminfo%my_proc_id)
      print*, 'size(st_pat,:) - dimensions:', size(st_pat,1), size(st_pat,2)
      print*, 'st_pat(:,:), assignment', st_pat(1:10,1)
      print*, 'In do while loop, nblks and nCells : ', nblks, nCells
      block => block%next
    enddo

!    allocate(xlon(nblks,maxblk))
!    allocate(xlat(nblks,maxblk))
    ierr = 0
!    call init_stochastic_physics(nVertLevels, & ! number of vertical levels
!         blksz, dt, &  ! physics timestep in seconds   
!         sppt_amp, & !pjp cloud perturbations (output)
!         input_nml_file, fn_nml, nlunit, &   ! namelist file
!         xlon, xlat, do_sppt, do_shum, do_skeb, &
!         lndp_type, n_var_lndp, & !land perturbation scheme and # of variables
!         use_zmtnblck, skeb_npass, & ! stochastic physics control variables
!         lndp_var_list, lndp_prt_list, & ! land perturbation var and print list
!         n_var_spp, spp_var_list, spp_prt_list, & ! spp var infomation
!         spp_stddev_cutoff, do_spp,  & ! spp flags                    
!         ak, bk, nthreads, & ! openMP threads available for ohysics
!         master, communicator, & ! master atmosphere processor MPI-rank, MPI communicator 
!         ierr)
!     if (ierr/=0)  then
!         write(6,*) 'call to init_stochastic_physics failed'
!         return
!     endif

  end subroutine stochastic_physics_patinit

  subroutine stochastic_physics_patadv (domain, levs, ierr)

    use stochastic_physics,  only: run_stochastic_physics
    implicit none

    type(domain_type),intent(inout):: domain
    integer, intent(in)  :: levs                 
    integer, intent(out) :: ierr

    integer :: n_threads, nb, maxblk, nblks, n, v
    logical :: param_update_flag

    ierr = 0

  end subroutine stochastic_physics_patadv

end module stochastic_physics_mpas
