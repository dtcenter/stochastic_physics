!>@brief Modifications were made to the routines in the module to work
!with the MPAS software framework.  

!>@brief The module 'stochastic_physics' is for initialization and running of
!! the stochastic physics random pattern generators
module stochastic_physics
use mpi
use kinddef, only : kind_phys, kind_dbl_prec
use mpas_pool_routines

implicit none

private

public :: init_stochastic_physics
public :: run_stochastic_physics

contains

!>@brief The subroutine 'init_stochastic_physics' initializes the stochastic
!!pattern genertors
!>@details It reads the stochastic physics namelist (nam_stoch and nam_sfcperts)
!allocates and polulates the necessary arrays

subroutine init_stochastic_physics(domain, levs, blksz, dtp, sppt_amp, &
         xlon, xlat, zk, mpicomm, mpiroot, iret) 
!\callgraph
!use stochy_internal_state_moa
use stochy_data_mod, only : init_stochdata,gg_lats,gg_lons,nsppt, &
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy, &
                            vfact_sppt,sl
use stochy_namelist_def
use spectral_transforms,only:colrad_a,latg,lonf,skeblevs
use mpi_wrapper, only : mpi_wrapper_initialize,mype,npes,is_rootpe, mp_bcst

implicit none

! Interface variables

type(domain_type),intent(inout)     :: domain
integer,              intent(in)    :: levs, mpiroot
integer,              intent(in)    :: mpicomm
integer,              intent(in)    :: blksz(:)
real(kind=kind_phys), intent(in)    :: dtp
real(kind=kind_phys), intent(out)   :: sppt_amp
real(kind=kind_phys), intent(in)    :: xlon(:,:)
real(kind=kind_phys), intent(in)    :: xlat(:,:)
!real(kind=kind_phys), intent(in)    :: ak(:), bk(:) 
real(kind=kind_phys), intent(in)    :: zk(:) 
integer, intent(out)                :: iret

! Local variables
real(kind=kind_phys), parameter     :: con_pi =4.0d0*atan(1.0d0)
integer :: nblks,len
real*8 :: dx
real, allocatable :: skeb_vloc(:)
integer :: k,kflip,latghf,blk,k2,v,i
integer :: k_top(2)

! Initialize MPI and OpenMP
call mpi_wrapper_initialize(mpiroot,mpicomm)

gis_stochy%nodes = npes
gis_stochy%mype=mype
gis_stochy%nx=maxval(blksz)
nblks = size(blksz)
gis_stochy%ny=nblks
rad2deg=180.0/con_pi

! ------------------------------------------

allocate(gis_stochy%len(nblks))
allocate(gis_stochy%parent_lons(gis_stochy%nx,gis_stochy%ny))
allocate(gis_stochy%parent_lats(gis_stochy%nx,gis_stochy%ny))
do blk=1,nblks
   len=blksz(blk)
   gis_stochy%parent_lons(1:len,blk)=xlon(blk,1:len)*rad2deg
   gis_stochy%parent_lats(1:len,blk)=xlat(blk,1:len)*rad2deg
   gis_stochy%len(blk)=len
enddo

! replace
INTTYP=0 ! bilinear interpolation
call init_stochdata(domain,mype,levs,dtp,iret)

if (is_rootpe()) then
  print*, 'return value from init_stochdata():', iret
endif

! update remaining model configuration parameters from namelist
!use_zmtnblck_out=use_zmtnblck

if  (.NOT. do_sppt) return

!allocate(sl(levs))
!do k=1,levs
!   sl(k)= 0.5*(ak(k)/101300.+bk(k)+ak(k+1)/101300.0+bk(k+1)) ! si are now sigmas
!    sl(k)=1.0 - real(k)/real(levs) + 1.0/real(levs)
!enddo

!print*, "zk:", zk
if (is_rootpe()) then
  print*, 'sppt_hgt_top1,2:', sppt_hgt_top1, sppt_hgt_top2 
  call find_ktop(zk, levs+1, sppt_hgt_top1,sppt_hgt_top2,k_top)
endif

if (is_rootpe()) print*, 'k_top:', k_top
call mp_bcst(k_top, 2) 
if (is_rootpe()) print*, 'after broadcast k_top:', k_top

if (do_sppt) then
   allocate(vfact_sppt(levs))
   do k=1,levs
      if (k .lt. k_top(1)) then
          vfact_sppt(k) = 1.0
      elseif (k .gt. k_top(2)) then
          vfact_sppt(k) = 0.0
      elseif (k .ge. k_top(1) .and. k .le. k_top(2)) then
         vfact_sppt(k) = 1.0 - (real(k)-real(k_top(1)))/(real(k_top(2)) - real(k_top(1)))
      endif
   enddo

!if (is_rootpe()) print*,'after vfact_sppt(k) computation:', vfact_sppt

if (is_rootpe()) print*, 'sppt_sfclimit:', sppt_sfclimit

   if (sppt_sfclimit) then
       do k=1,7
       vfact_sppt(k)=pbl_taper(k)
       enddo
   endif
   if (is_rootpe()) then
      do k=1,levs
         print *,'sppt vert profile',k,vfact_sppt(k)
      enddo
   endif
   sppt_amp=sqrt(SUM(sppt(1:nsppt)**2))
endif

!if (is_rootpe()) print*, 'sppt_amp:', sppt_amp

! get interpolation weights
! define gaussian grid lats and lons
latghf=latg/2
allocate(gg_lats(latg))
allocate(gg_lons(lonf))
do k=1,latghf
   gg_lats(k)=-1.0*colrad_a(latghf-k+1)*rad2deg
   gg_lats(latg-k+1)=-1*gg_lats(k)
enddo
dx=360.0/lonf
do k=1,lonf
  gg_lons(k)=dx*(k-1)
enddo
WLON=gg_lons(1)-(gg_lons(2)-gg_lons(1))
RNLAT=gg_lats(1)*2-gg_lats(2)

if (is_rootpe()) then
  print*, 'latg,lonf =', latg, lonf
  print*, 'gis_stochy%lats_nodes_a', gis_stochy%lats_nodes_a(:)
  print*, 'gis_stochy%lats_node_a', gis_stochy%lats_node_a
endif

!if (is_rootpe()) print*, 'exit init_stochastic_physics, iret:', iret
end subroutine init_stochastic_physics

subroutine run_stochastic_physics(levs, kdt, blksz, sppt_wts, sppt_wts_gg, iret)

!\callgraph
!use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_sppt,nsppt, gis_stochy,vfact_sppt 
use get_stochy_pattern_mod,only : get_random_pattern_scalar,get_random_pattern_vector, & 
                                  get_random_pattern_sfc,get_random_pattern_spp
use stochy_namelist_def, only : do_sppt,nssppt,sppt_logit 
use mpi_wrapper, only: is_rootpe
implicit none

! Interface variables
integer,                  intent(in) :: levs, kdt
integer,                  intent(in) :: blksz(:)
integer,                  intent(out) :: iret 
real(kind=kind_phys), intent(inout) :: sppt_wts(:,:,:)
real(kind=RKIND), intent(inout) :: sppt_wts_gg(:,:)  

real(kind_dbl_prec),allocatable :: tmp_wts(:,:)
integer :: k,v
integer j,ierr,i
integer :: nblks, blk, length
logical :: do_advance_pattern

if (.NOT. do_sppt) return

! Update number of threads in shared variables in spectral_layout_mod and set block-related variables
nblks = size(blksz)
if (is_rootpe()) then
  print*, 'kdt, nssppt, nblks, blkzs(1) =', kdt, nssppt,nblks,blksz(1)
!  print*, 'nsppt =', nsppt
!  print*, 'sppt_logit =', sppt_logit
endif

iret = 1
allocate(tmp_wts(gis_stochy%nx,gis_stochy%ny))
if (do_sppt) then
   if (mod(kdt,nssppt) == 1 .or. nssppt == 1) then
     if (is_rootpe()) print*, 'advance pattern', kdt, nssppt
      call get_random_pattern_scalar(rpattern_sppt,nsppt,gis_stochy,tmp_wts,sppt_wts_gg)
      DO blk=1,nblks
         length=blksz(blk)
         DO k=1,levs
            !sppt_wts(blk,1:len,k)=tmp_wts(blk,1:len)*vfact_sppt(k)
            sppt_wts(blk,1:length,k)=tmp_wts(1:length,blk)*vfact_sppt(k)
         ENDDO
         if (sppt_logit) sppt_wts(blk,:,:) = (2./(1.+exp(sppt_wts(blk,:,:))))-1.
         sppt_wts(blk,:,:) = sppt_wts(blk,:,:)+1.0
      ENDDO
      iret = 0
   endif
endif
deallocate(tmp_wts)

end subroutine run_stochastic_physics

subroutine find_ktop(zk, nvlsp1,  hgt_top1, hgt_top2, k_top)
  real(kind=kind_phys), intent(in) :: zk(nvlsp1), hgt_top1, hgt_top2 
  integer, intent(in) :: nvlsp1
  integer, intent(out) :: k_top(2)

  integer :: k

  do k = 1, nvlsp1-1
    if (zk(k) < hgt_top1 .and. zk(k+1) > hgt_top1) then
       k_top(1) = k
    endif
    if (zk(k) < hgt_top2 .and. zk(k+1) > hgt_top2) then
       k_top(2) = k+1
    endif
  enddo   
 
end subroutine find_ktop

end module stochastic_physics
