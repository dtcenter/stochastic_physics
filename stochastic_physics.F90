!>@brief Modifications were made to the routines in the module to work
!with the MPAS software framework.  

!>@brief The module 'stochastic_physics' is for initialization and running of
!! the stochastic physics random pattern generators
module stochastic_physics
use mpi_f08
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
!use stochy_data_mod, only : init_stochdata,gg_lats,gg_lons,nsppt, &
!                            rad2deg,INTTYP,wlon,rnlat,gis_stochy, &
!                            vfact_sppt,sl
use stochy_data_mod, only : init_stochdata,gg_lats,gg_lons,nsppt, &
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy, &
                            vfact_skeb,vfact_sppt,vfact_shum,skeb_vpts,skeb_vwts,sl, &
                            nspp, vfact_spp
use stochy_namelist_def
use spectral_transforms,only:colrad_a,latg,lonf,skeblevs
use mpi_wrapper, only : mpi_wrapper_initialize,mype,npes,is_rootpe, mp_bcst

implicit none

! Interface variables

type(domain_type),intent(inout)     :: domain
integer,              intent(in)    :: levs, mpiroot
type(MPI_Comm),       intent(in)    :: mpicomm
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
integer :: sppt_k_top(2), skeb_k_top(2)

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

if (is_rootpe()) then
  print*, 'do_sppt, do_skeb:', do_sppt, do_skeb
endif

iret = 0
if (do_sppt .eqv. .false.) then
  iret = ior(iret, 1)
endif
if (do_skeb .eqv. .false.) then
  iret = ior(iret, 2)
endif
  
if  (.NOT. do_sppt .AND. .NOT. do_skeb) then
  return
endif

!allocate(sl(levs))
!do k=1,levs
!   sl(k)= 0.5*(ak(k)/101300.+bk(k)+ak(k+1)/101300.0+bk(k+1)) ! si are now sigmas
!    sl(k)=1.0 - real(k)/real(levs) + 1.0/real(levs)
!enddo

!print*, "zk:", zk
if (is_rootpe()) then
  print*, 'sppt_hgt_top1,2:', sppt_hgt_top1, sppt_hgt_top2 
  print*, 'skeb_hgt_top1,2:', skeb_hgt_top1, skeb_hgt_top2 
  call find_ktop(zk, levs+1, sppt_hgt_top1,sppt_hgt_top2,sppt_k_top)
  call find_ktop(zk, levs+1, skeb_hgt_top1,skeb_hgt_top2,skeb_k_top)
endif

if (is_rootpe()) print*, 'sppt_k_top:', sppt_k_top
if (is_rootpe()) print*, 'skeb_k_top:', skeb_k_top

call mp_bcst(sppt_k_top, 2) 
call mp_bcst(skeb_k_top, 2) 

if (is_rootpe()) print*, 'after broadcast sppt_k_top:', sppt_k_top
if (is_rootpe()) print*, 'after broadcast skeb_k_top:', skeb_k_top

! initialization for sppt 
if (do_sppt) then
   allocate(vfact_sppt(levs))
   do k=1,levs
      if (k .lt. sppt_k_top(1)) then
          vfact_sppt(k) = 1.0
      elseif (k .gt. sppt_k_top(2)) then
          vfact_sppt(k) = 0.0
      elseif (k .ge. sppt_k_top(1) .and. k .le. sppt_k_top(2)) then
         vfact_sppt(k) = 1.0 - (real(k)-real(sppt_k_top(1)))/(real(sppt_k_top(2)) - real(sppt_k_top(1)))
      endif
   enddo

   if (is_rootpe()) print*,'sppt vert profile:', vfact_sppt
   if (is_rootpe()) print*, 'sppt_sfclimit:', sppt_sfclimit

   if (sppt_sfclimit) then
       do k=1,7
       vfact_sppt(k)=pbl_taper(k)
       enddo
   endif
   if (is_rootpe()) then
      do k=1,levs
         print *,'after pbl_tapering, vfact_sppt(k)',k,vfact_sppt(k)
      enddo
   endif
   sppt_amp=sqrt(SUM(sppt(1:nsppt)**2))
endif

if (is_rootpe()) print*, 'sppt_amp:', sppt_amp

! initialization for skeb 
if (do_skeb) then
   allocate(vfact_skeb(levs))
!   allocate(skeb_vloc(skeblevs)) ! local
!   allocate(skeb_vwts(levs,2)) ! save for later
!   allocate(skeb_vpts(levs,2)) ! save for later
   do k=1,levs
      if (k .lt. skeb_k_top(1)) then
          vfact_skeb(k) = 1.0
      elseif (k .gt. skeb_k_top(2)) then
          vfact_skeb(k) = 0.0
      elseif (k .ge. skeb_k_top(1) .and. k .le. skeb_k_top(2)) then
         vfact_skeb(k) = 1.0 - (real(k)-real(skeb_k_top(1)))/(real(skeb_k_top(2)) - real(skeb_k_top(1)))
      endif
      if (is_rootpe())  print *,'skeb vert profile',k,vfact_skeb(k)
   enddo
! calculate vertical interpolation weights
!   do k=1,skeblevs
!      skeb_vloc(k)=sl(1)-real(k-1)/real(skeblevs-1.0)*(sl(1)-sl(levs))
!   enddo
! surface
!skeb_vwts(1,2)=0
!skeb_vpts(1,1)=1
! top
!skeb_vwts(levs,2)=1
!skeb_vpts(levs,1)=skeblevs-2
! internal
!DO k=2,levs-1
!   DO k2=1,skeblevs-1
!      IF (sl(k) .LE. skeb_vloc(k2) .AND. sl(k) .GT. skeb_vloc(k2+1)) THEN
!        skeb_vpts(k,1)=k2
!        skeb_vwts(k,2)=(skeb_vloc(k2)-sl(k))/(skeb_vloc(k2)-skeb_vloc(k2+1))
!      ENDIF
!   ENDDO
!ENDDO
!deallocate(skeb_vloc)
!if (is_rootpe()) then
!DO k=1,levs
!   print*,'skeb vpts ',skeb_vpts(k,1),skeb_vwts(k,2)
!ENDDO
!endif
!skeb_vwts(:,1)=1.0-skeb_vwts(:,2)
!skeb_vpts(:,2)=skeb_vpts(:,1)+1.0

endif

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

if (is_rootpe()) print*, 'exit init_stochastic_physics, iret:', iret
end subroutine init_stochastic_physics

subroutine run_stochastic_physics(levs,kdt,blksz,sppt_wts,skebu_wts,skebv_wts,skebt_wts,wts_gg,iret)

!\callgraph
!use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_sppt,nsppt,gis_stochy,vfact_sppt,rpattern_skebuv,rpattern_skebth,nskeb,vfact_skeb 
use get_stochy_pattern_mod,only : get_random_pattern_scalar,get_random_pattern_vector, & 
                                  get_random_pattern_sfc,get_random_pattern_spp, &
                                  get_one_random_pattern_vector 
use stochy_namelist_def, only : do_sppt,do_skeb,nssppt,nsskeb,sppt_logit 
use mpi_wrapper, only: is_rootpe
implicit none

! Interface variables
integer,                  intent(in) :: levs, kdt
integer,                  intent(in) :: blksz(:)
integer,                  intent(out) :: iret 
real(kind=kind_phys), intent(inout) :: sppt_wts(:,:,:)
real(kind=kind_phys), intent(inout) :: skebu_wts(:,:,:)
real(kind=kind_phys), intent(inout) :: skebv_wts(:,:,:)
real(kind=kind_phys), intent(inout) :: skebt_wts(:,:,:)

real(kind=RKIND), intent(inout) :: wts_gg(:,:)  

real(kind_dbl_prec),allocatable :: tmp_wts(:,:),tmpu_wts(:,:),tmpv_wts(:,:),tmpl_wts(:,:,:),tmp_spp_wts(:,:,:)
integer :: k,v
integer j,ierr,i
integer :: nblks, blk, length
logical :: do_advance_pattern

if (.NOT. do_sppt) return

! Update number of threads in shared variables in spectral_layout_mod and set block-related variables
nblks = size(blksz)
if (is_rootpe()) then
  print*, 'kdt, nssppt, nblks, blkzs(1) =', kdt, nssppt,nblks,blksz(1)
  print*, 'nsppt =', nsppt
  print*, 'sppt_logit =', sppt_logit
endif

iret = 1
if (do_sppt) then
   allocate(tmp_wts(gis_stochy%nx,gis_stochy%ny))
   if (mod(kdt,nssppt) == 1 .or. nssppt == 1) then
     if (is_rootpe()) print*, 'advance pattern', kdt, nssppt
      call get_random_pattern_scalar(rpattern_sppt,nsppt,gis_stochy,tmp_wts,wts_gg)
      do blk=1,nblks
         length=blksz(blk)
         do k=1,levs
            !sppt_wts(blk,1:len,k)=tmp_wts(blk,1:len)*vfact_sppt(k)
            sppt_wts(blk,1:length,k)=tmp_wts(1:length,blk)*vfact_sppt(k)
         enddo
         if (sppt_logit) sppt_wts(blk,:,:) = (2./(1.+exp(sppt_wts(blk,:,:))))-1.
         sppt_wts(blk,:,:) = sppt_wts(blk,:,:)+1.0
      enddo
      iret = 0
   endif
   deallocate(tmp_wts)
endif
if (do_skeb) then
!vector pattern
   allocate(tmpu_wts(gis_stochy%nx,gis_stochy%ny))
   allocate(tmpv_wts(gis_stochy%nx,gis_stochy%ny))
   if (mod(kdt,nsskeb) == 1 .or. nsskeb == 1) then
      call get_one_random_pattern_vector(rpattern_skebuv,nskeb,gis_stochy,tmpu_wts,tmpv_wts)
      do blk=1,nblks
         length=blksz(blk)
         do k=1,levs
            skebu_wts(blk,1:length,k)=tmpu_wts(1:length,blk)*vfact_skeb(k)
            skebv_wts(blk,1:length,k)=tmpv_wts(1:length,blk)*vfact_skeb(k)
         enddo
      enddo
   endif
   deallocate(tmpu_wts)
   deallocate(tmpv_wts)
! scalar pattern
   allocate(tmp_wts(gis_stochy%nx,gis_stochy%ny))
   if (mod(kdt,nsskeb) == 1 .or. nsskeb == 1) then
      call get_random_pattern_scalar(rpattern_skebth,nskeb,gis_stochy,tmp_wts,wts_gg)
      do blk=1,nblks
         length=blksz(blk)
         do k=1,levs
            skebt_wts(blk,1:length,k)=tmp_wts(1:length,blk)*vfact_skeb(k)
         enddo
      enddo
      iret = 0
   endif
   deallocate(tmp_wts)
endif


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
