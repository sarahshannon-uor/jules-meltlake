
! *****************************COPYRIGHT*******************************
!
! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]
!
! *****************************COPYRIGHT*******************************
! SUBROUTINE ROUTE_LAKE_INFLOW-----------------------------------------
!
! Description:
!   Route lateral meltwater from the nonlake snowpack to the lake
!   and surface runoff.
!
! Method:
!   When lake water is exposed, transfer 50% of nonlake outflow
!   into the lake and route the remaining 50% to surface runoff.
!
!   When lake water is not exposed, route all nonlake outflow
!   to surface runoff.
!
!   Surface runoff is returned as a gridbox-average rate
!   (kg m-2 s-1).
!
! Code Owner: s.r.shannon@reading.ac.uk
!
!---------------------------------------------------------------------

MODULE route_lake_inflow_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                         &
  ModuleName='ROUTE_LAKE_INFLOW_MOD'

CONTAINS

SUBROUTINE route_lake_inflow(land_pts,       &
                             timestep,      &
                             nsurft,        &
                             surft_pts,     &
                             surft_index,   &
                             frac,          &
                             exposed_water, &
                             lake_inflow,   &
                             lake_depth_ml, &
                             lateral_runoff_gb)

USE water_constants_mod, ONLY: rho_water

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!---------------------------------------------------------------------
! Scalar arguments with intent(in)
!---------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                           &
  land_pts,                                                      &
    ! Number of land points.
  nsurft
    ! Number of surface tiles.

REAL(KIND=real_jlslsm), INTENT(IN) ::                           &
  timestep
    ! Timestep length (s).

!---------------------------------------------------------------------
! Array arguments with intent(in)
!---------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                           &
  surft_pts(nsurft),                                             &
    ! Number of points on each tile.
  surft_index(land_pts,nsurft)
    ! Land-point indices for each tile.

LOGICAL, INTENT(IN) ::                                           &
  exposed_water(land_pts,nsurft)
    ! True where lake water is exposed.

REAL(KIND=real_jlslsm), INTENT(IN) ::                            &
  frac(land_pts,nsurft),                                         &
    ! Fractional coverage of surface tiles.
  lake_inflow(land_pts,nsurft)
    ! Water leaving snowpack (kg m-2 per timestep).

!---------------------------------------------------------------------
! Array arguments with intent(inout)
!---------------------------------------------------------------------

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                        &
  lake_depth_ml(land_pts,nsurft)
    ! Lake depth (m).

!---------------------------------------------------------------------
! Array arguments with intent(out)
!---------------------------------------------------------------------

REAL(KIND=real_jlslsm), INTENT(OUT) ::                          &
  lateral_runoff_gb(land_pts)
    ! Gridbox surface runoff rate (kg m-2 s-1).

!---------------------------------------------------------------------
! Local parameters
!---------------------------------------------------------------------

INTEGER, PARAMETER ::                                           &
  lake_ml    = 9,                                               &
    ! Melt-lake tile index.
  nonlake_ml = 10
    ! Nonlake elevated-ice tile index.

!---------------------------------------------------------------------
! Local scalars
!---------------------------------------------------------------------

INTEGER ::                                                      &
  i,                                                            &
    ! Land-point index.
  k
    ! Tile-point index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ROUTE_LAKE_INFLOW'

!---------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,            &
                        zhook_in,zhook_handle)

lateral_runoff_gb(:) = 0.0

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i)

!-----------------------------------------------------------------------------
! Add lateral runoff from nonlake tile to lake before lid evolution. 
! Laternal water transport is only allowed if there is exposed water
!-----------------------------------------------------------------------------
DO k = 1, surft_pts(nonlake_ml)

   i = surft_index(k,nonlake_ml)

   IF (exposed_water(i,lake_ml) .AND. &
        frac(i,lake_ml) > 0.0) THEN

      ! Transfer 50% of nonlake outflow to exposed lake
      lake_depth_ml(i,lake_ml) = lake_depth_ml(i,lake_ml) &
           + 0.025 * frac(i,nonlake_ml) / frac(i,lake_ml) &
           * lake_inflow(i,nonlake_ml) / rho_water

      ! Route remaining 50% to surface runoff
      lateral_runoff_gb(i) = 0.5 * frac(i,nonlake_ml) &
           * lake_inflow(i,nonlake_ml) / timestep

   ELSE

      ! No exposed lake: all nonlake outflow becomes runoff
      lateral_runoff_gb(i) = frac(i,nonlake_ml) &
           * lake_inflow(i,nonlake_ml) / timestep

   END IF

END DO

!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,            &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE route_lake_inflow

END MODULE route_lake_inflow_mod
