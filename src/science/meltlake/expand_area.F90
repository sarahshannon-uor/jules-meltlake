! ****************************COPYRIGHT*******************************
!
! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]
!
! *****************************COPYRIGHT*******************************
! SUBROUTINE EXPAND_AREA-----------------------------------------------
!
! Update melt-lake tile fraction using the observed volume-area
! relationship:
!
!     V = a A**b
!
! Code Owner: s.r.shannon@reading.ac.uk
!---------------------------------------------------------------------

MODULE expand_area_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EXPAND_AREA_MOD'

CONTAINS

  SUBROUTINE expand_area(land_pts,                    &
                         nsurft,                      &
                         surft_pts,                   &
                         surft_index,                 &
                         exposed_water,               &
                         frac,                        &
                         lake_depth_ml,               &
                         lake_area_ml)

USE model_grid_mod, ONLY: grid_area_ij
    
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                       &
  land_pts,                                                                  &
  nsurft,                                                                    &
  surft_pts

INTEGER, INTENT(IN) ::                                                       &
  surft_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  frac(land_pts,nsurft),                                                     &
! Fractional coverage of each surface tile.
  lake_depth_ml(land_pts),                                                   &
! Mean lake depth over the lake tile (m).
  lake_area_ml(land_pts)
! lake area in (m2)

LOGICAL, INTENT(IN) ::                                                       &
  exposed_water(land_pts)
     
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------

INTEGER ::                                                                   &
  i,                                                                         &
  k

INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
  nonlake_ml = 10

REAL(KIND=real_jlslsm) ::                                                    &
  lake_volume,                                                               &
! Lake volume in gridbox (m3).
  new_lake_frac, &
! New lake tile fraction.
 lake_area
!-----------------------------------------------------------------------------
! Volume-area parameters for AWS18
!
! V = a_va * A**b_va
!
! V in m3, A in m2.
! aws18_lat = -66.411999
! aws18_lon = -63.435142
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), PARAMETER ::                                         &
  ! most obs   
  a_va         = 0.401,                                                      &
  b_va         = 1.018,                                                      &
  ! regional pooled good_fit
  !a_va         = 6.788,                                                     &
  !b_va         = 0.864,                                                     &
  ! synthetic
  ! a_va         = 0.01,                                                     &
  ! b_va         = 1.20,                                                     &
  gridbox_area = 3.0814E8
  ! ERA5 0.25 degree gridbox area at AWS18 (m2).
  


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND_AREA'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,lake_volume,new_lake_frac)

DO k = 1, surft_pts

   i = surft_index(k)

   IF (exposed_water(i) .AND. &
        lake_depth_ml(i) > 0.0) THEN
      
  
!-----------------------------------------------------------------------------
! Calculate current lake volume.
!-----------------------------------------------------------------------------

      lake_volume = lake_depth_ml(i) &
                  * frac(i,lake_ml)          &
                  * gridbox_area

      !WRITE(*,*) '------------------------------'
      !WRITE(*,*) 'gridbox_area = ', gridbox_area
      !WRITE(*,*) 'lake_depth before = ', lake_depth_ml(i)
      !WRITE(*,*) 'lake_frac before  = ', frac(i,lake_ml)
      !WRITE(*,*) 'lake_volume       = ', lake_volume

!-----------------------------------------------------------------------------
! Calculate lake area from V = a A^b.
!-----------------------------------------------------------------------------

      lake_area_ml(i) = (lake_volume / a_va)**(1.0 / b_va)

      !WRITE(*,*) 'lake_area from VA m2  = ', lake_area_ml(i)
      !WRITE(*,*) 'lake_area from VA km2 = ', lake_area_ml(i) / 1.0E6

      new_lake_frac = lake_area_ml(i) / gridbox_area

      !WRITE(*,*) 'new frac before cap = ', new_lake_frac

      new_lake_frac = MIN(1.0, new_lake_frac)

!-----------------------------------------------------------------------------
! Keep area, fraction and depth consistent.
!-----------------------------------------------------------------------------

      lake_area_ml(i) = &
         new_lake_frac * gridbox_area

      lake_depth_ml(i) = &
         lake_volume / lake_area_ml(i)

      frac(i,lake_ml)    = new_lake_frac
      frac(i,nonlake_ml) = 1.0 - new_lake_frac

      !WRITE(*,*) 'lake_area final m2  = ', lake_area_ml(i)
      !WRITE(*,*) 'lake_area final km2 = ', &
      !   lake_area_ml(i) / 1.0E6
      !WRITE(*,*) 'lake_frac after     = ', frac(i,lake_ml)
      !WRITE(*,*) 'lake_depth after    = ', lake_depth_ml(i)
      !WRITE(*,*) 'volume after        = ', &
      !   lake_depth_ml(i) * frac(i,lake_ml) * gridbox_area

   END IF

END DO


!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE expand_area

END MODULE expand_area_mod

