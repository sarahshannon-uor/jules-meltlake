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
! Description:
!
! Grau et al. 2025
! https://www.nature.com/articles/s41467-025-61798-8
!
! Uses cumulative water supply from the lake and non-lake elevated
! ice tiles to estimate lake area fraction and mean lake depth.
! 
! Code Owner: s.r.shannon@reading.ac.uk
! Subroutine Interface:
MODULE expand_area_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EXPAND_AREA_MOD'

CONTAINS

  SUBROUTINE expand_area(land_pts,                    &
                         nsurft,                      &
                         surft_pts,                   &
                         surft_index,                 &
                         frac,                        &
                         lake_inflow,                 &
                         dhdt_lake_snow_ml,           &
                         lid_rain_water_ml,           &
                         lid_snowmelt_water_ml,       &
                         lake_depth_ml,               &
                         grau_water_supply_ml,        &
                         grau_mean_lake_depth_ml,     &
                         exposed_water,               &
                         has_lid,                     &
                         has_vlid,                    &
                         has_lake,                    &
                         did_insert_lid)

USE model_time_mod, ONLY: timestep_number

USE water_constants_mod, ONLY:                                               &
  rho_ice,                                                                   &
! Density of solid ice (kg m-3).
  rho_water
! Density of pure water (kg m-3).

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                       &
  land_pts,                                                                  &
! Number of land points.
  nsurft,                                                                    &
! Number of surface tiles.
  surft_pts
! Number of lake tile points.

INTEGER, INTENT(IN) ::                                                       &
  surft_index(land_pts)
! Index of lake tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
  lake_inflow(land_pts,nsurft),                                              &
! Water supplied from snow/firn hydrology on each elevated-ice tile
! (kg m-2).
  dhdt_lake_snow_ml(land_pts,nsurft),                                        &
! Stefan lake-snow boundary movement (m ice equivalent per timestep).
  lid_rain_water_ml(land_pts,nsurft),                                        &
! Rainwater supplied on the lake/lid tile this timestep (m).
  lid_snowmelt_water_ml(land_pts,nsurft)
! Snowmelt water supplied on the lake/lid tile this timestep (m).

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  frac(land_pts,nsurft),                                                     &
! Fractional coverage of each surface tile.
  grau_water_supply_ml(land_pts),                                            &
! Cumulative gridbox-equivalent water supply used by Grau (m).
  grau_mean_lake_depth_ml(land_pts), &
! Grau-predicted mean water depth over lake-covered area (m).
  lake_depth_ml(land_pts,nsurft)

LOGICAL, INTENT(IN OUT) ::                                                    &
 exposed_water(land_pts,nsurft),                                              &
 has_lid(land_pts,nsurft),                                                    &
 has_vlid(land_pts,nsurft),                                                   &
 has_lake(land_pts,nsurft),                                                   &
 did_insert_lid(land_pts,nsurft)
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------

INTEGER ::                                                                   &
  i,                                                                         &
! Land-point index.
  k
! Lake tile point index.

INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
! Melt-lake tile index.
  nonlake_ml = 10
! Non-lake elevated-ice tile index.

REAL(KIND=real_jlslsm) ::                                                    &
  old_lake_frac,                                                             &
! Lake frac before Grau update .
  delta_water_supply,                                                        &
! Gridbox-equivalent water supplied this timestep (m).
  wd_max,                                                                    &
! Grau characteristic depression-storage depth (m).
  supply_ratio,                                                              &
! Grau dimensionless water-supply ratio.
  grau_frac, &
! Grau-predicted lake area fraction.
 lake_tile_supply
!-----------------------------------------------------------------------------
! Grau parameters
! Nearest satellite measurement to AWS 18
! Distance: 3.188 km
! Latitude: -66.411999
! Longitude: -63.435142
! H: 0.409942
! Sigma: 3.2476 m
!-----------------------------------------------------------------------------

REAL(KIND=real_jlslsm), PARAMETER ::                                         &
  hurst_ml = 0.41,                                                           &
! Hurst exponent of elevated-ice surface topography.
  sigma_ml = 3.25
! Standard deviation of elevated-ice surface elevation (m).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND_AREA'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,lake_tile_supply,old_lake_frac, &
!$OMP         delta_water_supply,wd_max,supply_ratio,grau_frac)
DO k = 1, surft_pts

   i = surft_index(k)

!-----------------------------------------------------------------------------
! Step 1: Transfer lateral runoff from non-lake tile to exposed lake.
! Use the existing tile fractions before any area changes. 
!-----------------------------------------------------------------------------

   IF (exposed_water(i,lake_ml) .AND. frac(i,lake_ml) > 0.0) THEN

      lake_depth_ml(i,lake_ml) = lake_depth_ml(i,lake_ml) &
           + frac(i,nonlake_ml) / frac(i,lake_ml) &
           * lake_inflow(i,nonlake_ml) / rho_water

   END IF

!-----------------------------------------------------------------------------
! Step 2: Calculate Grau water supply.
! Keep non-lake runoff in Grau's diagnostic supply even though
! it has also been transferred physically in step 1.
!-----------------------------------------------------------------------------

   IF (.NOT. has_lake(i,lake_ml)) THEN

      lake_tile_supply = lake_inflow(i,lake_ml) / rho_water

   ELSE IF (exposed_water(i,lake_ml)) THEN

      lake_tile_supply = dhdt_lake_snow_ml(i,lake_ml) &
           * rho_ice / rho_water

   ELSE IF (has_lid(i,lake_ml) .OR. has_vlid(i,lake_ml)) THEN

      lake_tile_supply = lid_rain_water_ml(i,lake_ml) &
           + lid_snowmelt_water_ml(i,lake_ml)

   ELSE

      lake_tile_supply = 0.0

   END IF

   delta_water_supply = &
        frac(i,lake_ml) * lake_tile_supply &
        + frac(i,nonlake_ml) &
        * lake_inflow(i,nonlake_ml) / rho_water

!-----------------------------------------------------------------------------
! Accumulate total water supply used by Grau.
!-----------------------------------------------------------------------------

   grau_water_supply_ml(i) = grau_water_supply_ml(i) &
                           + delta_water_supply

!-----------------------------------------------------------------------------
! Max depression storage, Eqn 1.
!-----------------------------------------------------------------------------

   wd_max = sigma_ml * (0.2 - 0.12 * hurst_ml**0.6)

!-----------------------------------------------------------------------------
! Supply ratio, Eqn 2.
!-----------------------------------------------------------------------------

   supply_ratio = grau_water_supply_ml(i) / wd_max

!-----------------------------------------------------------------------------
! Grau mean lake depth, Eqn 3. Diagnostic only
!-----------------------------------------------------------------------------

   grau_mean_lake_depth_ml(i) = &
        0.6 * sigma_ml * ERF(67.0 * supply_ratio) &
        * (1.0 - 0.41 * hurst_ml**0.6)

!-----------------------------------------------------------------------------
! Grau mean lake area fraction, Eqn 4.
!-----------------------------------------------------------------------------

   grau_frac = 0.13 * ERF(55.0 * supply_ratio) &
        * (1.0 - 0.13 * hurst_ml**1.3 &
        + 2.0 * ERF(0.08 * supply_ratio))

   grau_frac = MAX(0.0, MIN(1.0, grau_frac))

!-----------------------------------------------------------------------------
! Tile fraction changes 
!-----------------------------------------------------------------------------

   old_lake_frac = frac(i,lake_ml)
   
   frac(i,lake_ml) = grau_frac
   frac(i,nonlake_ml) = 1.0 - frac(i,lake_ml)

  IF (frac(i,lake_ml) > 0.0 .AND. old_lake_frac > 0.0) THEN

     lake_depth_ml(i,lake_ml) = lake_depth_ml(i,lake_ml) &
          * old_lake_frac / frac(i,lake_ml)

  END IF

!-----------------------------------------------------------------------------
! Reset Grau accumulation following complete lake refreezing.
!-----------------------------------------------------------------------------

   IF (did_insert_lid(i,lake_ml)) THEN

      grau_water_supply_ml(i) = 0.0
      grau_mean_lake_depth_ml(i) = 0.0

   END IF


!-----------------------------------------------------------------------------
! Diagnostics
!-----------------------------------------------------------------------------

   IF (delta_water_supply > 0.0) THEN

      PRINT *, '===================================================='


      PRINT *, 'EXPAND AREA'
      PRINT *, 't                     = ', timestep_number
      PRINT *, 'land point            = ', i
      PRINT *, 'lake inflow           = ',                         &
           lake_inflow(i,lake_ml) / rho_water
      PRINT *, 'nonlake inflow        = ',                         &
           lake_inflow(i,nonlake_ml) / rho_water
      PRINT *, 'lake bottom melt      = ',                         &
           dhdt_lake_snow_ml(i,lake_ml) * rho_ice / rho_water
      PRINT *, 'lid rain water        = ',                         &
           lid_rain_water_ml(i,lake_ml)
      PRINT *, 'lid snowmelt water    = ',                         &
           lid_snowmelt_water_ml(i,lake_ml)
      PRINT *, 'delta water supply    = ', delta_water_supply
      PRINT *, 'Grau water supply     = ', grau_water_supply_ml(i)
      PRINT *, 'supply ratio          = ', supply_ratio
      PRINT *, 'Grau mean lake depth  = ', grau_mean_lake_depth_ml(i)
      PRINT *, 'Grau fraction         = ', grau_frac
   
      PRINT *, 'lake frac             = ', frac(i,lake_ml)
      PRINT *, 'nonlake frac          = ', frac(i,nonlake_ml)
      !PRINT *, 'water before          = ', old_lake_frac * old_lake_depth
      !PRINT *, 'water after           = ', frac(i,lake_ml) * lake_depth_ml(i,lake_ml)


      PRINT *, '===================================================='
      !stop
   END IF

END DO

!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE expand_area

END MODULE expand_area_mod
