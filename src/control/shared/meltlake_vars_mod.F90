! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
!
! Module containing melt lake variables.
!
! Code Description:
!   Language: FORTRAN 90
!
! Code Owner: Please refer to ModuleLeaders.txt
!

MODULE meltlake_vars_mod


USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Implementation for field variables:
! Each variable is declared in both the 'data' TYPE and the 'pointer' type.
! Instances of these types are declared at at high level as required
! This is to facilitate advanced memory management features, which are generally
! not visible in the science code.
! Checklist for adding a new variable:
! -add to data_type
! -add to pointer_type
! -add to the allocate routine, passing in any new dimension sizes required
!  by argument (not via USE statement)
! -add to the deallocate routine
! -add to the assoc and nullify routines
!  Sarah copy fire variables 

TYPE :: meltlake_vars_data_type    
 
  REAL(KIND=real_jlslsm), ALLOCATABLE :: sfrac_ml(:,:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lfrac_ml(:,:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: refreeze_ml(:,:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: melt_ml(:,:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lake_depth_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lake_albedo_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lake_temp_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ice_lens_depth(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ice_lens_index(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lake_inflow(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: kdtdz_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ksnow0_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lid_depth_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: vlid_depth_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lid_temp_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lake_state_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: dhdt_lake_snow_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: dhdt_lid_lake_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lid_snow_depth_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lid_snow_temp_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: cold_puddle_hrs_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: snow_on_lid_melt_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: water_on_lid_depth_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: ei_surft_ml(:,:)
  REAL(KIND=real_jlslsm), ALLOCATABLE :: lake_frac_ml(:,:) ! not used yet
  
  
  LOGICAL, ALLOCATABLE :: exposed_water(:,:)
  LOGICAL, ALLOCATABLE :: has_lid(:,:)
  LOGICAL, ALLOCATABLE :: has_vlid(:,:)
  LOGICAL, ALLOCATABLE :: has_lake(:,:)
  LOGICAL, ALLOCATABLE :: did_insert_lid(:,:)
  LOGICAL, ALLOCATABLE :: snow_on_lid(:,:)
  
  
END TYPE meltlake_vars_data_type

TYPE :: meltlake_vars_type
  
  REAL(KIND=real_jlslsm), POINTER :: sfrac_ml(:,:,:)
  REAL(KIND=real_jlslsm), POINTER :: lfrac_ml(:,:,:)
  REAL(KIND=real_jlslsm), POINTER :: refreeze_ml(:,:,:)
  REAL(KIND=real_jlslsm), POINTER :: melt_ml(:,:,:)
  REAL(KIND=real_jlslsm), POINTER :: lake_depth_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lake_albedo_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lake_temp_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: ice_lens_depth(:,:)
  REAL(KIND=real_jlslsm), POINTER :: ice_lens_index(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lake_inflow(:,:)
  REAL(KIND=real_jlslsm), POINTER :: kdtdz_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: ksnow0_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lid_depth_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: vlid_depth_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lid_temp_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lake_state_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: dhdt_lake_snow_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: dhdt_lid_lake_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lid_snow_depth_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lid_snow_temp_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: cold_puddle_hrs_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: snow_on_lid_melt_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: water_on_lid_depth_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: ei_surft_ml(:,:)
  REAL(KIND=real_jlslsm), POINTER :: lake_frac_ml(:,:)
  
  LOGICAL, POINTER :: exposed_water(:,:)
  LOGICAL, POINTER :: has_lid(:,:)
  LOGICAL, POINTER :: has_vlid(:,:)
  LOGICAL, POINTER :: has_lake(:,:)
  LOGICAL, POINTER :: did_insert_lid(:,:)
  LOGICAL, POINTER :: snow_on_lid(:,:)
  
END TYPE meltlake_vars_type

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MELTLAKE_VARS_MOD'

CONTAINS

!===============================================================================
SUBROUTINE meltlake_vars_alloc(land_pts,nsurft,nsmax_ml,meltlake_vars_data)

!No USE statements other than Dr Hook
USE parkind1,               ONLY: jprb, jpim
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts, nsurft, nsmax_ml
TYPE(meltlake_vars_data_type), INTENT(IN OUT) :: meltlake_vars_data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE_VARS_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Allocate space for meltlake diagnostic variables
!-----------------------------------------------------------------------

ALLOCATE(meltlake_vars_data%sfrac_ml(land_pts,nsurft,nsmax_ml))
ALLOCATE(meltlake_vars_data%lfrac_ml(land_pts,nsurft,nsmax_ml))
ALLOCATE(meltlake_vars_data%refreeze_ml(land_pts,nsurft,nsmax_ml))
ALLOCATE(meltlake_vars_data%melt_ml(land_pts,nsurft,nsmax_ml))
ALLOCATE(meltlake_vars_data%lake_depth_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lake_albedo_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lake_temp_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%ice_lens_depth(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%ice_lens_index(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lake_inflow(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%kdtdz_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%ksnow0_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lid_depth_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%vlid_depth_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lid_temp_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lake_state_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%dhdt_lake_snow_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%dhdt_lid_lake_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lid_snow_depth_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lid_snow_temp_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%cold_puddle_hrs_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%snow_on_lid_melt_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%water_on_lid_depth_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%ei_surft_ml(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%lake_frac_ml(land_pts,nsurft))

ALLOCATE(meltlake_vars_data%exposed_water(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%has_lid(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%has_vlid(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%has_lake(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%did_insert_lid(land_pts,nsurft))
ALLOCATE(meltlake_vars_data%snow_on_lid(land_pts,nsurft))


meltlake_vars_data%sfrac_ml(:,:,:)         = 0.0
meltlake_vars_data%lfrac_ml(:,:,:)         = 0.0
meltlake_vars_data%refreeze_ml(:,:,:)      = 0.0
meltlake_vars_data%melt_ml(:,:,:)          = 0.0
meltlake_vars_data%lake_depth_ml(:,:)      = 0.0
meltlake_vars_data%lake_albedo_ml(:,:)     = 0.85  
meltlake_vars_data%lake_temp_ml(:,:)       = 273.15  
meltlake_vars_data%ice_lens_depth(:,:)     = -1.0
meltlake_vars_data%ice_lens_index(:,:)     = 0.0
meltlake_vars_data%lake_inflow(:,:)        = 0.0
meltlake_vars_data%kdtdz_ml(:,:)           = 0.0
meltlake_vars_data%ksnow0_ml(:,:)          = 0.0
meltlake_vars_data%lid_depth_ml(:,:)       = 0.0
meltlake_vars_data%vlid_depth_ml(:,:)      = 0.0
meltlake_vars_data%lid_temp_ml(:,:)        = 273.15
meltlake_vars_data%lake_state_ml(:,:)      = 0.0
meltlake_vars_data%dhdt_lake_snow_ml(:,:)  = 0.0
meltlake_vars_data%dhdt_lid_lake_ml(:,:)   = 0.0
meltlake_vars_data%lid_snow_depth_ml(:,:)  = 0.0
meltlake_vars_data%lid_snow_temp_ml(:,:)   = 273.15
meltlake_vars_data%cold_puddle_hrs_ml(:,:) = 0.0
meltlake_vars_data%snow_on_lid_melt_ml(:,:)= 0.0
meltlake_vars_data%water_on_lid_depth_ml(:,:)= 0.0
meltlake_vars_data%lake_frac_ml(:,:)       = 0.0
meltlake_vars_data%ei_surft_ml(:,:)        = 0.0

meltlake_vars_data%exposed_water(:,:)      = .FALSE.
meltlake_vars_data%has_lid(:,:)            = .FALSE.
meltlake_vars_data%has_vlid(:,:)           = .FALSE.
meltlake_vars_data%has_lake(:,:)           = .FALSE.
meltlake_vars_data%did_insert_lid(:,:)     = .FALSE.
meltlake_vars_data%snow_on_lid(:,:)        = .FALSE.


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE meltlake_vars_alloc


!===============================================================================
SUBROUTINE meltlake_vars_dealloc(meltlake_vars_data)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(meltlake_vars_data_type), INTENT(IN OUT) :: meltlake_vars_data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE_VARS_DEALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Deallocate space for melt lake diagnostic variables
!-----------------------------------------------------------------------
DEALLOCATE(meltlake_vars_data%sfrac_ml)
DEALLOCATE(meltlake_vars_data%lfrac_ml)
DEALLOCATE(meltlake_vars_data%refreeze_ml)
DEALLOCATE(meltlake_vars_data%melt_ml)
DEALLOCATE(meltlake_vars_data%lake_depth_ml)
DEALLOCATE(meltlake_vars_data%lake_albedo_ml)
DEALLOCATE(meltlake_vars_data%lake_temp_ml)
DEALLOCATE(meltlake_vars_data%ice_lens_depth)
DEALLOCATE(meltlake_vars_data%ice_lens_index)
DEALLOCATE(meltlake_vars_data%lake_inflow)
DEALLOCATE(meltlake_vars_data%kdtdz_ml)
DEALLOCATE(meltlake_vars_data%ksnow0_ml)
DEALLOCATE(meltlake_vars_data%lid_depth_ml)
DEALLOCATE(meltlake_vars_data%vlid_depth_ml)
DEALLOCATE(meltlake_vars_data%lid_temp_ml)
DEALLOCATE(meltlake_vars_data%lake_state_ml)
DEALLOCATE(meltlake_vars_data%dhdt_lake_snow_ml)
DEALLOCATE(meltlake_vars_data%dhdt_lid_lake_ml)
DEALLOCATE(meltlake_vars_data%lid_snow_depth_ml)
DEALLOCATE(meltlake_vars_data%lid_snow_temp_ml)
DEALLOCATE(meltlake_vars_data%cold_puddle_hrs_ml)
DEALLOCATE(meltlake_vars_data%snow_on_lid_melt_ml)
DEALLOCATE(meltlake_vars_data%water_on_lid_depth_ml)
DEALLOCATE(meltlake_vars_data%ei_surft_ml)
DEALLOCATE(meltlake_vars_data%lake_frac_ml)

DEALLOCATE(meltlake_vars_data%exposed_water)
DEALLOCATE(meltlake_vars_data%has_lid)
DEALLOCATE(meltlake_vars_data%has_vlid)
DEALLOCATE(meltlake_vars_data%has_lake)
DEALLOCATE(meltlake_vars_data%did_insert_lid)
DEALLOCATE(meltlake_vars_data%snow_on_lid)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE meltlake_vars_dealloc

!===============================================================================
SUBROUTINE meltlake_vars_assoc(meltlake_vars, meltlake_vars_data)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(meltlake_vars_type), INTENT(IN OUT) :: meltlake_vars
TYPE(meltlake_vars_data_type), INTENT(IN OUT), TARGET :: meltlake_vars_data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE_VARS_ASSOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL meltlake_vars_nullify(meltlake_vars)

meltlake_vars%sfrac_ml => meltlake_vars_data%sfrac_ml
meltlake_vars%lfrac_ml => meltlake_vars_data%lfrac_ml
meltlake_vars%refreeze_ml => meltlake_vars_data%refreeze_ml
meltlake_vars%melt_ml => meltlake_vars_data%melt_ml
meltlake_vars%lake_depth_ml => meltlake_vars_data%lake_depth_ml
meltlake_vars%lake_albedo_ml => meltlake_vars_data%lake_albedo_ml
meltlake_vars%lake_temp_ml => meltlake_vars_data%lake_temp_ml
meltlake_vars%ice_lens_depth => meltlake_vars_data%ice_lens_depth
meltlake_vars%ice_lens_index => meltlake_vars_data%ice_lens_index
meltlake_vars%lake_inflow => meltlake_vars_data%lake_inflow
meltlake_vars%kdtdz_ml => meltlake_vars_data%kdtdz_ml
meltlake_vars%ksnow0_ml => meltlake_vars_data%ksnow0_ml
meltlake_vars%lid_depth_ml => meltlake_vars_data%lid_depth_ml
meltlake_vars%vlid_depth_ml => meltlake_vars_data%vlid_depth_ml
meltlake_vars%lid_temp_ml => meltlake_vars_data%lid_temp_ml
meltlake_vars%lake_state_ml => meltlake_vars_data%lake_state_ml
meltlake_vars%dhdt_lake_snow_ml => meltlake_vars_data%dhdt_lake_snow_ml
meltlake_vars%dhdt_lid_lake_ml => meltlake_vars_data%dhdt_lid_lake_ml
meltlake_vars%lid_snow_depth_ml => meltlake_vars_data%lid_snow_depth_ml
meltlake_vars%lid_snow_temp_ml => meltlake_vars_data%lid_snow_temp_ml
meltlake_vars%cold_puddle_hrs_ml => meltlake_vars_data%cold_puddle_hrs_ml
meltlake_vars%snow_on_lid_melt_ml => meltlake_vars_data%snow_on_lid_melt_ml
meltlake_vars%water_on_lid_depth_ml => meltlake_vars_data%water_on_lid_depth_ml
meltlake_vars%ei_surft_ml => meltlake_vars_data%ei_surft_ml
meltlake_vars%lake_frac_ml => meltlake_vars_data%lake_frac_ml

meltlake_vars%exposed_water => meltlake_vars_data%exposed_water
meltlake_vars%has_lid => meltlake_vars_data%has_lid
meltlake_vars%has_vlid => meltlake_vars_data%has_vlid
meltlake_vars%has_lake => meltlake_vars_data%has_lake
meltlake_vars%did_insert_lid => meltlake_vars_data%did_insert_lid
meltlake_vars%snow_on_lid => meltlake_vars_data%snow_on_lid

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE meltlake_vars_assoc

!===============================================================================
SUBROUTINE meltlake_vars_nullify(meltlake_vars)

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(meltlake_vars_type), INTENT(IN OUT) :: meltlake_vars

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE_VARS_NULLIFY'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(meltlake_vars%sfrac_ml)
NULLIFY(meltlake_vars%lfrac_ml)
NULLIFY(meltlake_vars%refreeze_ml)
NULLIFY(meltlake_vars%melt_ml)
NULLIFY(meltlake_vars%lake_depth_ml)
NULLIFY(meltlake_vars%lake_albedo_ml)
NULLIFY(meltlake_vars%lake_temp_ml)
NULLIFY(meltlake_vars%ice_lens_depth)
NULLIFY(meltlake_vars%ice_lens_index)
NULLIFY(meltlake_vars%lake_inflow)
NULLIFY(meltlake_vars%kdtdz_ml)
NULLIFY(meltlake_vars%ksnow0_ml)
NULLIFY(meltlake_vars%lid_depth_ml)
NULLIFY(meltlake_vars%vlid_depth_ml)
NULLIFY(meltlake_vars%lake_state_ml)
NULLIFY(meltlake_vars%dhdt_lake_snow_ml)
NULLIFY(meltlake_vars%dhdt_lid_lake_ml)
NULLIFY(meltlake_vars%lid_snow_depth_ml)
NULLIFY(meltlake_vars%lid_snow_temp_ml)
NULLIFY(meltlake_vars%cold_puddle_hrs_ml)
NULLIFY(meltlake_vars%snow_on_lid_melt_ml)
NULLIFY(meltlake_vars%water_on_lid_depth_ml)
NULLIFY(meltlake_vars%ei_surft_ml)
NULLIFY(meltlake_vars%lake_frac_ml)

NULLIFY(meltlake_vars%exposed_water)
NULLIFY(meltlake_vars%has_lid)
NULLIFY(meltlake_vars%has_vlid)
NULLIFY(meltlake_vars%has_lake)
NULLIFY(meltlake_vars%did_insert_lid)
NULLIFY(meltlake_vars%snow_on_lid)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE meltlake_vars_nullify

END MODULE meltlake_vars_mod
