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
! V in m3, A in m2.
! aws18_lat = -66.411999
! aws18_lon = -63.435142
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
                         nsnow,                       &
                         sice_ml,                     &
                         sliq_ml,                     &
                         tsnow_ml,                    &
                         ds_ml,                       &
                         snowdepth,                   &
                         snow_surft,                  &
                         rho_snow_ml,                 &
                         rho_snow_grnd,               &
                         exposed_water,               &
                         frac,                        &
                         lake_depth_ml,               &
                         lake_area_ml,                &
                         sice_imbal_ml,               &
                         sliq_imbal_ml,               &
                         en_imbal_ml,                 &
                         immediate_correction,         &
                         merge_snowpacks,              &
                         rgrain,                       &
                         rgrainl_ml)
         
USE model_time_mod,     ONLY: timestep_number
USE jules_meltlake_mod, ONLY: nsmax_ml

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
  ! Index of lake tile points.

INTEGER, INTENT(IN OUT) :: nsnow(land_pts,nsurft)
  ! Number of snow layers.
  
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  sice_ml(land_pts,nsurft,nsmax_ml),                                         &
  sliq_ml(land_pts,nsurft,nsmax_ml),                                         &
  tsnow_ml(land_pts,nsurft,nsmax_ml),                                        &
  ds_ml(land_pts,nsurft,nsmax_ml),                                           &
  snowdepth(land_pts,nsurft),                                                &
  snow_surft(land_pts,nsurft),                                               &
  rho_snow_ml(land_pts,nsurft,nsmax_ml),                                      &
  rho_snow_grnd(land_pts,nsurft),                                            &
  rgrain(land_pts,nsurft),                                                   &
    ! Snow surface grain size (microns).
  rgrainl_ml(land_pts,nsurft,nsmax_ml),                                      &
    ! Snow layer grain size (microns).
  frac(land_pts,nsurft),                                                     &
   ! Fractional coverage of each surface tile.
  lake_depth_ml(land_pts),                                                   &
   ! Mean lake depth over the lake tile (m).
  lake_area_ml(land_pts),                                                    &
   ! lake area in (m2)
  sice_imbal_ml(land_pts),                                                    &
   ! Per-call compensating snow ice mass adjustment (kg m-2 gridbox).
  sliq_imbal_ml(land_pts),                                                   &
   ! Per-call compensating snow liquid water mass adjustment (kg m-2 gridbox).
  en_imbal_ml(land_pts)
   ! Per-call compensating snowpack enthalpy adjustment (J m-2 gridbox).

LOGICAL, INTENT(IN) :: merge_snowpacks
  ! TRUE: area-weight the lake and nonlake snowpacks at each fraction change.

LOGICAL, INTENT(IN) :: immediate_correction
  ! TRUE: scale the snowpack on the growing tile.
  ! FALSE: no correction; diagnostics remain zero.

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
  delta_frac


!Larger va_b: slower relative increase in lake area with volume.
!Smaller va_b: faster relative increase in lake area with volume.
!Larger va_a: Deeper lakes for a given area, or smaller areas for a given volume.
!Smaller va_a: Shallower lakes for a given area, or larger areas for a given volume.
! Reference curve: 10 m depth at approximately 20% coverage

REAL(KIND=real_jlslsm), PARAMETER ::                                         &
 a_va = 0.045,                                                           &
 b_va = 1.2,                                                                 &
 gridbox_area = 3.0814E8
  ! ERA5 0.25 degree gridbox area at AWS18 (m2).
  
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND_AREA'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Per-call adjustments; zero when correction is not applied.
sice_imbal_ml(:) = 0.0
sliq_imbal_ml(:) = 0.0
en_imbal_ml(:)   = 0.0

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,lake_volume,new_lake_frac,delta_frac)

DO k = 1, surft_pts

   i = surft_index(k)

!-----------------------------------------------------------------------------
! Only evolve frac if there is exposed water. Don't change if there is a lid
!-----------------------------------------------------------------------------
   
   IF (exposed_water(i) .AND. &
        lake_depth_ml(i) > 0.0) THEN

!-----------------------------------------------------------------------------
! Lake exists so get volume and area
!-----------------------------------------------------------------------------
      IF (frac(i,lake_ml) > 0.0) THEN

         lake_volume = lake_depth_ml(i) * frac(i,lake_ml) * gridbox_area

         lake_area_ml(i) = (lake_volume / a_va)**(1.0 / b_va)

!-----------------------------------------------------------------------------
! Lake area frac no < 0 yet
! Get area from lake_depth = a_va * lake_area**(b_va - 1)
!-----------------------------------------------------------------------------         
      ELSE

         lake_area_ml(i) = &
              (lake_depth_ml(i) / a_va)**(1.0 / (b_va - 1.0))

         lake_volume = lake_depth_ml(i) * lake_area_ml(i)

      END IF

!-----------------------------------------------------------------------------
! Update the lake frac using area 
!-----------------------------------------------------------------------------
      new_lake_frac = lake_area_ml(i) / gridbox_area

      new_lake_frac = MIN(1.0, new_lake_frac)

!-----------------------------------------------------------------------------
! Overwrite lake depth with depth consistant with volume and area
!-----------------------------------------------------------------------------

      lake_area_ml(i) = &
         new_lake_frac * gridbox_area

      IF (lake_area_ml(i) <= 0.0) THEN
         PRINT *, 'EXPAND_AREA: zero lake area'
         PRINT *, 'Timestep = ', timestep_number
         PRINT *, 'Lake depth = ', lake_depth_ml(i)
         PRINT *, 'Lake fraction = ', frac(i,lake_ml)
         PRINT *, 'Lake volume = ', lake_volume
         PRINT *, 'Lake area = ', lake_area_ml(i)
         PRINT *, 'Exposed water = ', exposed_water(i)
         STOP 'EXPAND_AREA: zero lake area'
      END IF

      lake_depth_ml(i) = &
         lake_volume / lake_area_ml(i)

!-----------------------------------------------------------------------------
! The compensating adjustment to masses from changing frac
!-----------------------------------------------------------------------------

      delta_frac = new_lake_frac - frac(i,lake_ml)

      IF (delta_frac /= 0.0) THEN

         IF (merge_snowpacks) THEN

            CALL merge_snowpacks_frac(i, land_pts, nsurft, nsnow,       &
                 sice_ml, sliq_ml, tsnow_ml, ds_ml, snowdepth,         &
                 snow_surft, rho_snow_ml, rho_snow_grnd,              &
                 rgrain, rgrainl_ml, frac(i,lake_ml))

         ELSE IF (immediate_correction) THEN

            ! Diagnose the first lake expansion after complete freeze-over
            IF (frac(i,lake_ml) == 0.0) THEN
               PRINT *, '-----------------------------------------'
               PRINT *, 'LAKE REACTIVATION, grid point:', i
               PRINT *, 'Old lake fraction:', frac(i,lake_ml)
               PRINT *, 'New lake fraction:', new_lake_frac
               PRINT *, 'BEFORE correction'
               PRINT *, 'Lake snow depth:', snowdepth(i,lake_ml)
               PRINT *, 'Nonlake snow depth:', snowdepth(i,nonlake_ml)
               PRINT *, 'Lake snow ice:', SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml)))
               PRINT *, 'Nonlake snow ice:', SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))
            END IF

            ! Diagnose adjustments and apply the correction.
            CALL correct_snowpack_frac(i, land_pts, nsurft, nsnow,           &
                 sice_ml, sliq_ml, tsnow_ml, ds_ml, snowdepth, snow_surft,   &
                 rho_snow_ml, rho_snow_grnd, frac(i,lake_ml), new_lake_frac, &
                 sice_imbal_ml(i), sliq_imbal_ml(i), en_imbal_ml(i))

            IF (frac(i,lake_ml) == 0.0) THEN
               PRINT *, 'AFTER correction'
               PRINT *, 'Lake snow depth:', snowdepth(i,lake_ml)
               PRINT *, 'Nonlake snow depth:', snowdepth(i,nonlake_ml)
               PRINT *, 'Lake snow ice:', SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml)))
               PRINT *, 'Nonlake snow ice:', SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))
               PRINT *, 'Snow ice correction:', sice_imbal_ml(i)
               PRINT *, '-----------------------------------------'
            END IF


         END IF

      END IF

!-----------------------------------------------------------------------------
! Update the fractions
!-----------------------------------------------------------------------------      
      frac(i,lake_ml)    = new_lake_frac
      frac(i,nonlake_ml) = 1.0 - new_lake_frac

     
   END IF ! exposed_water

END DO !surft_pts


!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE expand_area


!-----------------------------------------------------------------------------
! Conserve gridbox snow ice, liquid water and sensible energy when frac changes.
! Scale the snowpack on the growing tile, preserving the ice mass per unit
! snow depth.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Area-weight the two snowpacks using the OLD tile fractions. Assign the
! same snowpack to both tiles, so new fractions conserve gridbox snow ice,
! liquid water and sensible energy. Grain size is weighted by snow ice mass.
! Corresponding layer indices are merged, as in the earlier merge method.
!-----------------------------------------------------------------------------
SUBROUTINE merge_snowpacks_frac(i, land_pts, nsurft, nsnow,             &
     sice_ml, sliq_ml, tsnow_ml, ds_ml, snowdepth, snow_surft,          &
     rho_snow_ml, rho_snow_grnd, rgrain, rgrainl_ml, old_frac)

USE jules_meltlake_mod, ONLY: nsmax_ml
USE jules_snow_mod, ONLY: r0
USE water_constants_mod, ONLY: hcapi, hcapw, tm
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                 &
  i,                                                                    &
    ! Land point index.
  land_pts,                                                             &
    ! Number of land points.
  nsurft
    ! Number of surface tiles.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                   &
  old_frac
    ! Lake fraction before the area change.

INTEGER, INTENT(IN OUT) ::                                              &
  nsnow(land_pts,nsurft)
    ! Number of active snow layers on each tile.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                               &
  sice_ml(land_pts,nsurft,nsmax_ml),                                    &
    ! Snow ice mass (kg m-2 tile).
  sliq_ml(land_pts,nsurft,nsmax_ml),                                    &
    ! Snow liquid mass (kg m-2 tile).
  tsnow_ml(land_pts,nsurft,nsmax_ml),                                   &
    ! Layer temperature (K).
  ds_ml(land_pts,nsurft,nsmax_ml),                                      &
    ! Layer thickness (m).
  snowdepth(land_pts,nsurft),                                           &
    ! Total snow depth (m).
  snow_surft(land_pts,nsurft),                                          &
    ! Snow mass (kg m-2 tile).
  rho_snow_ml(land_pts,nsurft,nsmax_ml),                                &
    ! Layer bulk density (kg m-3).
  rho_snow_grnd(land_pts,nsurft),                                       &
    ! Snowpack bulk density (kg m-3).
  rgrain(land_pts,nsurft),                                              &
    ! Surface snow grain size (microns).
  rgrainl_ml(land_pts,nsurft,nsmax_ml)
    ! Layer snow grain size (microns).

INTEGER, PARAMETER ::                                                  &
  lake_ml = 9,                                                          &
    ! Lake tile index.
  nonlake_ml = 10
    ! Nonlake tile index.

INTEGER ::                                                             &
  n,                                                                    &
    ! Source layer index.
  nl,                                                                   &
    ! Number of lake snow layers.
  nn,                                                                   &
    ! Number of nonlake snow layers.
  nmerge,                                                               &
    ! Number of merged snow layers.
  surft
    ! Tile receiving the merged profile.

REAL(KIND=real_jlslsm) ::                                              &
  lake_frac,                                                            &
    ! Old lake tile fraction.
  nonlake_frac,                                                         &
    ! Old nonlake tile fraction.
  sice_l, sice_n,                                                       &
    ! Ice masses of source layers (kg m-2 tile).
  sliq_l, sliq_n,                                                       &
    ! Liquid masses of source layers (kg m-2 tile).
  ds_l, ds_n,                                                           &
    ! Source layer thicknesses (m).
  energy_l, energy_n,                                                   &
    ! Source layer sensible energies (J m-2 tile).
  sice_new, sliq_new, ds_new, energy_new, cap_new,                      &
    ! Weighted layer mass, depth, energy and heat capacity.
  grain_num
    ! Ice-mass-weighted grain size numerator (microns kg m-2).

REAL(KIND=real_jlslsm) ::                                              &
  sice_merge(nsmax_ml),                                                &
    ! Merged snow ice mass (kg m-2 tile).
  sliq_merge(nsmax_ml),                                                &
    ! Merged liquid water mass (kg m-2 tile).
  ds_merge(nsmax_ml),                                                  &
    ! Merged layer thickness (m).
  tsnow_merge(nsmax_ml),                                               &
    ! Merged layer temperature (K).
  rgrain_merge(nsmax_ml)
    ! Merged layer grain size (microns).

lake_frac = old_frac
nonlake_frac = 1.0 - old_frac

nl = nsnow(i,lake_ml)
nn = nsnow(i,nonlake_ml)

! A bulk snow mass with nsnow=0 has no resolved layer energy to merge.
IF ((nl == 0 .AND. snow_surft(i,lake_ml) > 0.0) .OR. &
    (nn == 0 .AND. snow_surft(i,nonlake_ml) > 0.0)) THEN
   STOP 'MERGE_SNOWPACKS_FRAC: unresolved bulk snow mass'
END IF

nmerge = 0

sice_merge(:) = 0.0
sliq_merge(:) = 0.0
ds_merge(:) = 0.0
tsnow_merge(:) = tm
rgrain_merge(:) = r0

DO n = 1, MAX(nl,nn)

   sice_l = 0.0
   sliq_l = 0.0
   ds_l = 0.0
   energy_l = 0.0
   sice_n = 0.0
   sliq_n = 0.0
   ds_n = 0.0
   energy_n = 0.0
   grain_num = 0.0

   IF (n <= nl) THEN
      sice_l = sice_ml(i,lake_ml,n)
      sliq_l = sliq_ml(i,lake_ml,n)
      ds_l = ds_ml(i,lake_ml,n)
      energy_l = (hcapi*sice_l + hcapw*sliq_l) * &
                 (tsnow_ml(i,lake_ml,n)-tm)
      IF (sice_l > 0.0) grain_num = &
           lake_frac*sice_l*rgrainl_ml(i,lake_ml,n)
   END IF

   IF (n <= nn) THEN
      sice_n = sice_ml(i,nonlake_ml,n)
      sliq_n = sliq_ml(i,nonlake_ml,n)
      ds_n = ds_ml(i,nonlake_ml,n)
      energy_n = (hcapi*sice_n + hcapw*sliq_n) * &
                 (tsnow_ml(i,nonlake_ml,n)-tm)
      IF (sice_n > 0.0) grain_num = grain_num + &
           nonlake_frac*sice_n*rgrainl_ml(i,nonlake_ml,n)
   END IF

   sice_new = lake_frac*sice_l + nonlake_frac*sice_n
   sliq_new = lake_frac*sliq_l + nonlake_frac*sliq_n
   ds_new = lake_frac*ds_l + nonlake_frac*ds_n
   energy_new = lake_frac*energy_l + nonlake_frac*energy_n

   ! Do not introduce empty trailing layers when one tile has no snow.
   IF (ds_new <= 0.0) THEN
      IF (sice_new+sliq_new > 0.0) THEN
         STOP 'MERGE_SNOWPACKS_FRAC: positive mass with zero thickness'
      END IF
      CYCLE
   END IF

   nmerge = nmerge + 1
   sice_merge(nmerge) = sice_new
   sliq_merge(nmerge) = sliq_new
   ds_merge(nmerge) = ds_new

   cap_new = hcapi*sice_new + hcapw*sliq_new
   IF (cap_new > 0.0) THEN
      tsnow_merge(nmerge) = tm + energy_new/cap_new
   ELSE IF (energy_new /= 0.0) THEN
      STOP 'MERGE_SNOWPACKS_FRAC: zero heat capacity'
   END IF

   IF (sice_new > 0.0) THEN
      rgrain_merge(nmerge) = grain_num/sice_new
   END IF

END DO

DO surft = lake_ml, nonlake_ml

   nsnow(i,surft) = nmerge
   sice_ml(i,surft,:) = 0.0
   sliq_ml(i,surft,:) = 0.0
   ds_ml(i,surft,:) = 0.0
   rho_snow_ml(i,surft,:) = 0.0
   tsnow_ml(i,surft,:) = tm
   rgrainl_ml(i,surft,:) = r0

   IF (nmerge > 0) THEN
      sice_ml(i,surft,1:nmerge) = sice_merge(1:nmerge)
      sliq_ml(i,surft,1:nmerge) = sliq_merge(1:nmerge)
      ds_ml(i,surft,1:nmerge) = ds_merge(1:nmerge)
      tsnow_ml(i,surft,1:nmerge) = tsnow_merge(1:nmerge)
      rgrainl_ml(i,surft,1:nmerge) = rgrain_merge(1:nmerge)

      rho_snow_ml(i,surft,1:nmerge) = &
           (sice_merge(1:nmerge)+sliq_merge(1:nmerge)) / &
           ds_merge(1:nmerge)
      rgrain(i,surft) = rgrain_merge(1)
   ELSE
      rgrain(i,surft) = r0
   END IF

   snowdepth(i,surft) = SUM(ds_merge(1:nmerge))
   snow_surft(i,surft) = SUM(sice_merge(1:nmerge)) + &
                         SUM(sliq_merge(1:nmerge))
   rho_snow_grnd(i,surft) = 0.0
   IF (snowdepth(i,surft) > 0.0) THEN
      rho_snow_grnd(i,surft) = snow_surft(i,surft) / snowdepth(i,surft)
   END IF

END DO

END SUBROUTINE merge_snowpacks_frac

SUBROUTINE correct_snowpack_frac(i, land_pts, nsurft, nsnow,             &
     sice_ml, sliq_ml, tsnow_ml, ds_ml, snowdepth, snow_surft,           &
     rho_snow_ml, rho_snow_grnd, old_frac, new_frac,                     &
     sice_imbal, sliq_imbal, en_imbal)

USE jules_meltlake_mod, ONLY: nsmax_ml
USE water_constants_mod, ONLY: hcapi, hcapw, tm
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                  &
  i,                                                                    &
    ! Land point index
  land_pts,                                                             &
    ! Number of land points
  nsurft
    ! Number of surface tiles

INTEGER, INTENT(IN) ::                                                  &
  nsnow(land_pts,nsurft)
    ! Number of active snow layers on each tile

REAL(KIND=real_jlslsm), INTENT(IN) ::                                   &
  old_frac,                                                             &
    ! Lake tile fraction before the area change
  new_frac
    ! Lake tile fraction after the area change

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                  &
  sice_imbal,                                                           &
    ! Snow ice mass correction (kg m-2 gridbox)
  sliq_imbal,                                                           &
    ! Snow liquid water mass correction (kg m-2 gridbox)
  en_imbal
    ! Snow sensible energy correction (J m-2 gridbox)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                               &
  sice_ml(land_pts,nsurft,nsmax_ml),                                    &
    ! Snow ice mass in each layer (kg m-2 tile)
  sliq_ml(land_pts,nsurft,nsmax_ml),                                    &
    ! Snow liquid water mass in each layer (kg m-2 tile)
  tsnow_ml(land_pts,nsurft,nsmax_ml),                                   &
    ! Snow layer temperature (K)
  ds_ml(land_pts,nsurft,nsmax_ml),                                      &
    ! Snow layer thickness (m)
  snowdepth(land_pts,nsurft),                                           &
    ! Total snow depth (m)
  snow_surft(land_pts,nsurft),                                          &
    ! Total snow mass (kg m-2 tile)
  rho_snow_ml(land_pts,nsurft,nsmax_ml),                                &
    ! Bulk snow layer density (kg m-3)
  rho_snow_grnd(land_pts,nsurft)
    ! Bulk snowpack density (kg m-3)

INTEGER, PARAMETER ::                                                  &
  lake_ml = 9,                                                          &
    ! Lake tile index
  nonlake_ml = 10
    ! Nonlake tile index

INTEGER ::                                                             &
  n,                                                                    &
    ! Snow layer index
  nl,                                                                   &
    ! Number of active lake snow layers
  nn,                                                                   &
    ! Number of active nonlake snow layers
  surft
    ! Index of surface tile receiving the correction

REAL(KIND=real_jlslsm) ::                                              &
  dfrac,                                                                &
    ! Change in lake tile fraction
  target_frac,                                                          &
    ! Updated fraction of the tile receiving the correction
  dsice,                                                                &
    ! Snow ice mass correction (kg m-2 gridbox)
  dsliq,                                                                &
    ! Snow liquid water mass correction (kg m-2 gridbox)
  denergy,                                                              &
    ! Sensible energy correction (J m-2 gridbox)
  energy_lake,                                                          &
    ! Lake snowpack sensible energy (J m-2 tile)
  energy_nonlake
    ! Nonlake snowpack sensible energy (J m-2 tile)

REAL(KIND=real_jlslsm) ::                                              &
  sice_total,                                                           &
    ! Original total snow ice mass on target tile (kg m-2 tile)
  sliq_total,                                                           &
    ! Original total liquid water mass on target tile (kg m-2 tile)
  sice_total_new,                                                       &
    ! Corrected total snow ice mass on target tile (kg m-2 tile)
  sliq_total_new,                                                       &
    ! Corrected total liquid water mass on target tile (kg m-2 tile)
  sice_scale,                                                           &
    ! Multiplicative snow ice mass scaling factor
  sliq_scale
    ! Multiplicative snow liquid water mass scaling factor

REAL(KIND=real_jlslsm) ::                                              &
  old_mass,                                                             &
    ! Original total snow layer mass (kg m-2 tile)
  new_mass,                                                             &
    ! Corrected total snow layer mass (kg m-2 tile)
  cap,                                                                  &
    ! Corrected snow layer heat capacity (J m-2 K-1)
  total_cap,                                                            &
    ! Corrected target snowpack heat capacity (J m-2 K-1)
  dtemp
    ! Uniform snow layer temperature correction (K)

REAL(KIND=real_jlslsm) ::                                              &
  layer_energy_old(nsmax_ml)
    ! Original sensible energy of each target snow layer (J m-2 tile)


dfrac = new_frac - old_frac

!-----------------------------------------------------------------------------
! No frac change, exit
!-----------------------------------------------------------------------------
IF (dfrac == 0.0) RETURN

!-----------------------------------------------------------------------------
! Number of active snow layers in lake and nonlake snowpacks
!-----------------------------------------------------------------------------
nl = nsnow(i,lake_ml)
nn = nsnow(i,nonlake_ml)

!-----------------------------------------------------------------------------
! Calculate corrections required to conserve gridbox snow ice, liquid water
! and sensible energy
!-----------------------------------------------------------------------------
dsice = dfrac * (SUM(sice_ml(i,nonlake_ml,1:nn)) - &
                 SUM(sice_ml(i,lake_ml,1:nl)))

dsliq = dfrac * (SUM(sliq_ml(i,nonlake_ml,1:nn)) - &
                 SUM(sliq_ml(i,lake_ml,1:nl)))

energy_lake = SUM((hcapi*sice_ml(i,lake_ml,1:nl) + &
                   hcapw*sliq_ml(i,lake_ml,1:nl)) * &
                  (tsnow_ml(i,lake_ml,1:nl) - tm))

energy_nonlake = SUM((hcapi*sice_ml(i,nonlake_ml,1:nn) + &
                      hcapw*sliq_ml(i,nonlake_ml,1:nn)) * &
                     (tsnow_ml(i,nonlake_ml,1:nn) - tm))

denergy = dfrac * (energy_nonlake - energy_lake)

sice_imbal = dsice
sliq_imbal = dsliq
en_imbal   = denergy

!-----------------------------------------------------------------------------
! Identify the snowpack receiving the correction
! Lake growth: correct lake tile
! Lake contraction or freeze-over: correct nonlake tile
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Identify the snowpack receiving the correction
! Preserve lake snowpack memory when lake area increases from zero
!-----------------------------------------------------------------------------
!IF (dfrac > 0.0 .AND. old_frac == 0.0) THEN
!   surft = nonlake_ml
!   target_frac = 1.0 - new_frac

!ELSE IF (dfrac > 0.0) THEN
!   surft = lake_ml
!   target_frac = new_frac

!ELSE
!   surft = nonlake_ml
!   target_frac = 1.0 - new_frac
!END IF



IF (dfrac > 0.0) THEN
   surft = lake_ml
   target_frac = new_frac
ELSE
   surft = nonlake_ml
   target_frac = 1.0 - new_frac
END IF

IF (old_frac == 0.0 .AND. new_frac > 0.0) THEN
   PRINT *, 'CORRECT_SNOWPACK_FRAC'
   PRINT *, 'Tile receiving correction:', surft
   PRINT *, 'Lake tile:', lake_ml
   PRINT *, 'Nonlake tile:', nonlake_ml
   PRINT *, 'Target fraction:', target_frac
END IF

!-----------------------------------------------------------------------------
! Stop if the target tile has no active snow layers
!-----------------------------------------------------------------------------
IF (target_frac <= 0.0 .OR. nsnow(i,surft) <= 0) THEN
   PRINT *, 'CORRECT_SNOWPACK_FRAC: invalid target tile'
   PRINT *, 'Grid point = ', i
   PRINT *, 'Target tile = ', surft
   PRINT *, 'Old lake fraction = ', old_frac
   PRINT *, 'New lake fraction = ', new_frac
   PRINT *, 'Target fraction = ', target_frac
   PRINT *, 'Target nsnow = ', nsnow(i,surft)
   STOP 'CORRECT_SNOWPACK_FRAC: invalid target tile'
END IF

!-----------------------------------------------------------------------------
! Add gridbox mass corrections to the target snowpack
!-----------------------------------------------------------------------------
sice_total = SUM(sice_ml(i,surft,1:nsnow(i,surft)))
sliq_total = SUM(sliq_ml(i,surft,1:nsnow(i,surft)))

sice_total_new = sice_total + dsice/target_frac
sliq_total_new = sliq_total + dsliq/target_frac

IF (sice_total_new < 0.0 .OR. sliq_total_new < 0.0) THEN
   STOP 'CORRECT_SNOWPACK_FRAC: negative snow mass'
END IF

!-----------------------------------------------------------------------------
! Calculate multiplicative scaling factors
!-----------------------------------------------------------------------------
sice_scale = 0.0
sliq_scale = 0.0

IF (sice_total > 0.0) sice_scale = sice_total_new/sice_total
IF (sliq_total > 0.0) sliq_scale = sliq_total_new/sliq_total

!-----------------------------------------------------------------------------
! Save original layer energies, then scale layer masses and thicknesses
!-----------------------------------------------------------------------------
DO n = 1, nsnow(i,surft)

   old_mass = sice_ml(i,surft,n) + sliq_ml(i,surft,n)

   layer_energy_old(n) = (hcapi*sice_ml(i,surft,n) + &
                          hcapw*sliq_ml(i,surft,n)) * &
                         (tsnow_ml(i,surft,n) - tm)

   sice_ml(i,surft,n) = sice_ml(i,surft,n)*sice_scale
   sliq_ml(i,surft,n) = sliq_ml(i,surft,n)*sliq_scale

  ! IF (n == 1) THEN
  !    IF (sice_total == 0.0) sice_ml(i,surft,n) = sice_total_new
  !    IF (sliq_total == 0.0) sliq_ml(i,surft,n) = sliq_total_new
  ! END IF

   new_mass = sice_ml(i,surft,n) + sliq_ml(i,surft,n)

!-----------------------------------------------------------------------------
! Scale layer thickness to preserve ice mass per unit snow depth
!-----------------------------------------------------------------------------
   IF (old_mass > 0.0) THEN
      ds_ml(i,surft,n) = ds_ml(i,surft,n)*sice_scale
   ELSE IF (new_mass > 0.0) THEN
      STOP 'CORRECT_SNOWPACK_FRAC: mass added to empty layer'
   END IF

!-----------------------------------------------------------------------------
! Recalculate bulk density using corrected mass and layer thickness
!-----------------------------------------------------------------------------
   IF (ds_ml(i,surft,n) > 0.0) THEN
      rho_snow_ml(i,surft,n) = new_mass/ds_ml(i,surft,n)
   END IF

END DO

!-----------------------------------------------------------------------------
! Restore original layer sensible energies and apply energy correction
!-----------------------------------------------------------------------------
total_cap = SUM(hcapi*sice_ml(i,surft,1:nsnow(i,surft)) + &
                hcapw*sliq_ml(i,surft,1:nsnow(i,surft)))

dtemp = (denergy/target_frac)/total_cap

DO n = 1, nsnow(i,surft)

   cap = hcapi*sice_ml(i,surft,n) + hcapw*sliq_ml(i,surft,n)

   IF (cap > 0.0) THEN
      tsnow_ml(i,surft,n) = tm + layer_energy_old(n)/cap + dtemp
   ELSE IF (layer_energy_old(n) /= 0.0) THEN
      STOP 'CORRECT_SNOWPACK_FRAC: lost layer energy'
   END IF

END DO

!-----------------------------------------------------------------------------
! Update target snowpack depth, total mass and bulk density
!-----------------------------------------------------------------------------
snowdepth(i,surft) = SUM(ds_ml(i,surft,1:nsnow(i,surft)))

snow_surft(i,surft) = SUM(sice_ml(i,surft,1:nsnow(i,surft))) + &
                      SUM(sliq_ml(i,surft,1:nsnow(i,surft)))

IF (snowdepth(i,surft) > 0.0) &
   rho_snow_grnd(i,surft) = snow_surft(i,surft)/snowdepth(i,surft)


END SUBROUTINE correct_snowpack_frac

END MODULE expand_area_mod
