! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE ADJUST_ICE_LENS_DEPTH-----------------------------------------------

! Description:
!  Bury the lens if snowpack grows (accumulation)
!  Make lens shallower if snowpack decreses (compaction, melting, sublimination)
!  Set -1 values for nonsense values

! Subroutine Interface:
MODULE adjust_ice_lens_depth_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ADJUST_ICE_LENS_DEPTH_MOD'

CONTAINS

  SUBROUTINE adjust_ice_lens_depth ( land_pts, surft_pts, surft_index,         &
       dz_snowdepth, snowdepth, ice_lens_depth, ice_lens_index)


USE model_time_mod, ONLY: timestep_number
    
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Number of land points.
  surft_pts
    ! Number of tile points.
  
!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_index(land_pts)   ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  dz_snowdepth(land_pts),                                                      &
    ! Snow depth change (m).
  snowdepth(land_pts)
    ! Snowdepth of current timestep (m) 

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
   ice_lens_depth(land_pts),                                                   &
      ! Ice lens depth (m).
   ice_lens_index(land_pts)
      ! -1 if lens does not exist. prob should change to integer
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  k,                                                                           &
    ! Tile point index.
  i
    ! Land point index
 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ADJUST_ICE_LENS_DEPTH'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
!$OMP DO SCHEDULE(STATIC)
DO k = 1,surft_pts
   i = surft_index(k)
  
   IF (ice_lens_depth(i) >= 0.0) THEN
      ice_lens_depth(i) = ice_lens_depth(i) + dz_snowdepth(i)
   END IF

   ! If surface lowered past the lens, lens is exposed/removed
   IF (ice_lens_depth(i) < 0.0) THEN
      print *, 'surface lowered past the lens, lens is exposed/removed'
      ice_lens_depth(i) = -1.0
      ice_lens_index(i) = -1.0
   END IF

   ! If lens lies deeper than the snowpack, clear it
   IF (ice_lens_depth(i) > snowdepth(i)) THEN
      print *, 'lens depth is deeper than the snowpack depth, resetting lens'
      ice_lens_depth(i) = -1.0
      ice_lens_index(i) = -1.0
   END IF



   ! If surface lowered past the lens, lens is exposed/removed
   !IF (ice_lens_depth(i) < 0.0) THEN
     ! print *, 'lens is exposed'
      !ice_lens_depth(i) = -1.0
    !  ice_lens_index(i) = 0.0
   !END IF
      
   !If lens lies deeper than the snowpack, clear it
   !IF (ice_lens_depth(i) > snowdepth(i)) THEN
   !   print *, 'lens depth is deeper than the snowpack depth'
   !   stop
   !   ice_lens_depth(i) = -1.0
   !   ice_lens_index(i) = -1.0
   !END IF
      
END DO
!$OMP END DO


!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE adjust_ice_lens_depth
END MODULE adjust_ice_lens_depth_mod
