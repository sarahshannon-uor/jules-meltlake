! ****************************rCOPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE -----------------------------------------------
! Description:
!     Calling routine for calculate_lake_bottom_melt
! Method:
!     Increment lake depth using excess water from snowpack
!     Get lake albedo for nc output
!     Get lake temp from tile surface temp using 4/3 law
!     Find boundary change between lake bottom and snowpack top using
!     Stefan condition.
!     Calculate the requested Stefan melt at the lake-snow boundary.
!     Snow mass and geometry are adjusted later in the snow routine,
!     after layersnow has created the JULES snow-layer geometry.
!     flux_lower heat from lake to snowpack interface
!     kdtdz heat conducted from snow interface to top layer of snowpack
! Code Owner: s.r.shannon@reading.ac.uk
! Subroutine Interface:
MODULE calculate_lake_bottom_melt_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALCULATE_LAKE_BOTTOM_MELT_MOD'

CONTAINS

  SUBROUTINE calculate_lake_bottom_melt(        &
       land_pts, & !IN
       timestep,                 & !IN
       nsurft,                   & !IN
       surft_pts,                & !IN
       surft_index,              & !IN
       nsnow,                    & !IN 
       lake_inflow,              & !IN
       ksnow0_ml,                & !IN
       tsnow0_ml,                & !IN
       lw_down_surft,            & !IN  
       tstar_surft,              & !IN
       sw_surft,                 & !IN  
       lid_depth_ml,             & !IN
       vlid_depth_ml,            & !IN
       snow_surft,               & !IN/OUT
       snowdepth,                & !IN/OUT
       ls_snow,                  & !IN/OUT
       con_snow,                 & !IN/OUT
       ls_rain,                  & !IN/OUT
       con_rain,                 & !IN/OUT
       lake_depth_ml,            & !IN/OUT
       lake_albedo_ml,           & !IN/OUT
       lake_temp_ml,             & !IN/OUT
       has_lake,                 & !IN/OUT
       exposed_water,            & !IN/OUT
       has_lid,                  & !IN/OUT
       has_vlid,                 & !IN/OUT
       sice_ml,                  & !IN/OUT
       sliq_ml,                  & !IN/OUT
       ds_ml,                    & !IN/OUT
       cold_puddle_hrs_ml,       & !IN/OUT 
       kdtdz_ml,                 & !OUT
       dhdt_lake_snow_ml,        &  !OUT
       !Ancil info (IN)
       l_lice_point,             & !IN (land_pts)
       l_lice_surft)               !IN (ntype)  
                     
USE model_time_mod, ONLY: timestep_number

USE ereport_mod, ONLY: ereport

USE water_constants_mod,     ONLY:                                           &
 rho_water,                                                                  &
  ! Density of pure water (kg/m3).
 hcapw,                                                                      &
  ! Specific heat capacity of water (J/kg/K).
 lf,                                                                         &
  ! Latent heat of fusion at 0degC (J kg-1).
 rho_ice,                                                                    &
  ! Density of solid ice (kg/m3).
 tm
 ! Temperature at which fresh water freezes and ice melts (K).


!USE jules_snow_mod, ONLY:                                                      &
!  snow_hcon
  ! Thermal conductivity of lying snow (Watts per m per K) snow=0.265, ice=2.2 
  
USE jules_meltlake_mod, ONLY: l_meltlake, nsmax_ml
USE jules_surface_mod,  ONLY: l_elev_land_ice
USE jules_surface_types_mod, ONLY: ntype
!USE jules_water_tracers_mod, ONLY: l_wtrac_jls

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
   land_pts,                                                                   &
     ! Total number of land points.
   nsurft
   ! Number of land tiles.
   
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
   timestep
    ! Timestep length (s).


REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  sw_surft(land_pts,nsurft),                                                   &                                                   
    ! Net shortwave radiation on tile (W/m2). Lake albedo used
  lw_down_surft(land_pts,nsurft),                                              &
    ! Surface downward LW radiation on tiles (W/m2)
  tstar_surft(land_pts,nsurft),                                                &
    ! Tile surface temperature (K)
  lake_inflow(land_pts,nsurft),                                                &
    ! melt water from snowpack (kgm-2)
  tsnow0_ml(land_pts,nsurft),                                                  &
    ! snowpack top level temp (K)  
  ksnow0_ml(land_pts,nsurft),                                                  &
    ! thermal conductivity of top snow layer (W m⁻¹ K⁻¹)
  lid_depth_ml(land_pts,nsurft),                                               &
    ! depth of frozen lid (m)
  vlid_depth_ml(land_pts,nsurft) 
    ! depth of frozen virtual lid (m)
   
!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_pts(nsurft),                                                           &
    ! Number of tile points.
  surft_index(land_pts,nsurft),                                                &
    ! Index of tile points.
  nsnow(land_pts,nsurft)                                                        
    ! Number of snow layers.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
   lake_depth_ml(land_pts,nsurft),                                             &
    ! Convective rainfall rate (kg/m2/s).
   lake_albedo_ml(land_pts,nsurft),                                            &
    ! Albedo of lake
   lake_temp_ml(land_pts,nsurft),                                              &
    ! Temperature of melt lake (K)
   snow_surft(land_pts,nsurft),                                                &
    ! Snow mass on tiles (kgm-2)
   con_rain(land_pts),                                                         &
    ! Convective rainfall rate (kg/m2/s).
   ls_rain(land_pts),                                                          &
    ! Large-scale rainfall fall rate (kg/m2/s).
   con_snow(land_pts),                                                         &
    ! Convective frozen rainfall rate (kg/m2/s).
   ls_snow(land_pts),                                                          &
    ! Large-scale frozen precip fall rate (kg/m2/s).
   sice_ml(land_pts,nsurft,nsmax_ml),                                          &
    ! Ice content of snow layers (kg/m2)
   sliq_ml(land_pts,nsurft,nsmax_ml),                                          &
    ! Liquid content of snow layers (kg/m2)
   ds_ml(land_pts, nsurft, nsmax_ml),                                          &
    ! snowpack top level depth (m)
   cold_puddle_hrs_ml(land_pts,nsurft),                                        & 
    ! Number of accum hours with lake_depth < 0.1m a cold orphan puddle
   snowdepth(land_pts,nsurft)
    ! Snowdepth (m)
   
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
   kdtdz_ml(land_pts,nsurft),                                                  &
     ! Conductive heat flux from lake into snowpack. pass this to snow module
     ! so the heat flux will be adjusted if there is exposed water 
   dhdt_lake_snow_ml(land_pts,nsurft)
     ! Stefan boundary movement (ms-1 ice equivalent) 

LOGICAL, INTENT(IN OUT) ::                                                     &
   has_lake(land_pts,nsurft),                                                  &
   exposed_water(land_pts,nsurft),                                             &
   has_lid(land_pts,nsurft),                                                   &
   has_vlid(land_pts,nsurft)

!ancil_info (IN)
LOGICAL, INTENT(IN) ::                                                         &
     l_lice_point(land_pts),                                                   &
     l_lice_surft(ntype)

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  i,                                                                           &
    ! Land point index and loop counter.
  j,                                                                           &
    ! Tile pts loop counter.
  k,                                                                           &
    ! Tile number.
  n,                                                                           &                                                                           
  ! Tile loop counter.
  ns
   ! snow layer loop counter
!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  expon_term,                                                                  &
    ! exponential term in albedo eqn (m)
  ex,                                                                          &
    ! exponential term in albedo eqn (m)
      ! Conductive heat flux into snowpack from lake
  flux_lower, &
    ! flux lake and ice
  flux_upper, &
    ! flux between atmosphere and lake
  delta_t, &
  dt, &
  dTdt, &
  sw_absorb, &
  dhdt, &
   ! Lake bottom-snowpack surface boundary change (m per sec ice equivalent)
  dh_ice, &
   ! Boundary change from Stefan condition (m of ice per timestep)
  dh_water, &
   ! Water equiv of boundary change from Stefan condition (m of water per timestep)
  dm, &
   ! mass to be removed from snowpack from Stefan boundary retreat (kgm-2)
  flux_upper_diag, &
  dh_remain, &
  frac_melt, &
  dsice, &
  dsliq, &
  ds_old, &
  rain_add, &
  snow_add, &
  dsice_tot, &
  dsliq_tot
  

LOGICAL :: cold_lingering_shallow_puddle
REAL(KIND=real_jlslsm) ::                                                 &
  cold_shallow_puddle_hours(land_pts,nsurft)

! --- diagnostic to check how much lake depth is from inflow versus stefan
!REAL, SAVE       :: cum_inflow_m = 0.0
!REAL, SAVE       :: cum_stefan_m = 0.0

REAL(KIND=real_jlslsm), PARAMETER ::                                      &
     vlid_seed_depth = 0.0001,                                             &
     tau_water       = 0.36,                                              &
                     ! shortwave attenuation coeff (m-1)
                     ! Table 1, panchromatic absorption coefficient,
                     ! Pope et al. (2016), The Cryosphere,
                     ! doi:10.5194/tc-10-15-2016
                     ! 0.5 light penatrates deeper,  2-3 shallow heating
     tau_ice         = 1.5,                                               &
                     ! black ice extinction coefficient
                     ! rough value taking into account:
                     ! https://tc.copernicus.org/articles/15/1931/2021/#section3
                     ! Cooper, M.G., et al., 2021. Spectral attenuation coefficients from
                     ! measurements of light transmission in bare ice on the Greenland Ice Sheet.
                     ! The Cryosphere, 15(4), pp.1931-1953.
     emis_water      = 0.98,                                             &
     Jturb           = 1.907e-5,                                         &
                 ! Turbulent heat flux factor (ms⁻¹ K⁻¹/3) Eqn 16 Buzzard 
     lid_min_depth  = 0.1,                                               &
     lake_min_depth = 0.1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALCULATE_LAKE_BOTTOM_MELT'

INTEGER              :: errcode            ! Error code
CHARACTER(LEN=80)    :: ERRMSG             ! Error message

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(SHARED)                                             &
!$OMP PRIVATE(k,i,n,expon_term,ex,flux_lower,flux_upper,delta_t,dt,dTdt,      &
!$OMP         sw_absorb,dhdt,dh_ice,dh_water,dm,flux_upper_diag,dh_remain,    &
!$OMP         frac_melt,dsice,dsliq,ds_old,rain_add,snow_add,dsice_tot,       &
!$OMP         dsliq_tot,cold_lingering_shallow_puddle)
DO n = 1,nsurft
   
   IF (l_elev_land_ice .AND. l_lice_surft(n)) THEN
      
      DO k = 1,surft_pts(n)
         i = surft_index(k,n)


!-----------------------------------------------------------------------------
! Get lake states. exposed_water flag used in snowpack to modify heat flux into
! snowpack top. 
!-----------------------------------------------------------------------------
   CALL update_lake_states(i,n)

!-----------------------------------------------------------------------------
! No lake yet grow one using meltwater 
!-----------------------------------------------------------------------------
   IF (.NOT. has_lake(i,n)) THEN
      IF (lake_inflow(i,n) > 0.0) THEN
                     
         lake_depth_ml(i,n) = lake_depth_ml(i,n) +                          &
              lake_inflow(i,n) / rho_water
            
      END IF
   END IF

   !cum_inflow_m = cum_inflow_m + lake_inflow(i) / rho_water

!-----------------------------------------------------------------------------
! Re-calculate lake state after lake growth
!-----------------------------------------------------------------------------
   CALL update_lake_states(i,n)

!---------------------------------------------------------------------
! Mass-conserving fate for shallow uncovered puddle after inflow stops
! Edge case: small puddle lingers if inflow stops 
! Also no rain/snow is added to lake depth if < 0.1m leaving the linerging
! puddle. 
! Refreeze back into the snowpack if the puddle lingers for 48 hours 
! Note: I should add latent heat and change tsnow(i,1) - come back to this 
!---------------------------------------------------------------------

   cold_lingering_shallow_puddle = (lake_depth_ml(i,n) > 0.0) .AND.        &
        (lake_depth_ml(i,n) < lake_min_depth) .AND.                        &
        (.NOT. has_lid(i,n)) .AND. (.NOT. has_vlid(i,n))                     &
        .AND. (lake_inflow(i,n) <= 1.0e-12)                                &
        .AND. (tstar_surft(i,n) < tm)

   

   IF (cold_lingering_shallow_puddle) THEN
      
      cold_puddle_hrs_ml(i,n) = cold_puddle_hrs_ml(i,n) + timestep/3600.0

      IF (cold_puddle_hrs_ml(i,n) >= 48.0) THEN
         
         WRITE(*,*) '--- REFREEZE LINGERING SHALLOW PUDDLE ---'
         WRITE(*,'(A,I8)')    'timestep                = ', timestep_number
         WRITE(*,'(A,I8)')    'i                       = ', i
         WRITE(*,'(A,F12.6)') 'cold_puddle_hours       = ', cold_puddle_hrs_ml(i,n)
         WRITE(*,'(A,F12.6)') 'tstar_surft (C)         = ', tstar_surft(i,n) - 273.15
         WRITE(*,'(A,F12.6)') 'lake_inflow (kg m-2)    = ', lake_inflow(i,n)
         WRITE(*,'(A,F12.6)') 'lake_depth before (m)   = ', lake_depth_ml(i,n)
         WRITE(*,'(A,F12.6)') 'snow_surft before       = ', snow_surft(i,n)
         WRITE(*,'(A,F12.6)') 'ds_ml(1) before         = ', ds_ml(i,n,1)
         WRITE(*,'(A,F12.6)') 'sice_ml(1) before       = ', sice_ml(i,n,1)
         WRITE(*,'(A,F12.6)') 'sliq_ml(1) before       = ', sliq_ml(i,n,1)
         WRITE(*,'(A,F12.6)') 'snowdepth before        = ', snowdepth(i,n)
         
         sice_ml(i,n,1)  = sice_ml(i,n,1)  + rho_water * lake_depth_ml(i,n)
         ds_ml(i,n,1)    = ds_ml(i,n,1)    + lake_depth_ml(i,n) * rho_water / rho_ice
         snow_surft(i,n) = snow_surft(i,n) + rho_water * lake_depth_ml(i,n)

         snowdepth(i,n) = 0.0
         DO ns = 1, nsnow(i,n)
            snowdepth(i,n) = snowdepth(i,n) + ds_ml(i,n,ns)
         END DO

         lake_depth_ml(i,n) = 0.0
         lake_temp_ml(i,n)  = tm
         cold_puddle_hrs_ml(i,n) = 0.0

       WRITE(*,'(A,F12.6)') 'lake_depth after (m)    = ', lake_depth_ml(i,n)
       WRITE(*,'(A,F12.6)') 'snow_surft after        = ', snow_surft(i,n)
       WRITE(*,'(A,F12.6)') 'ds_ml(1) after          = ', ds_ml(i,n,1)
       WRITE(*,'(A,F12.6)') 'sice_ml(1) after        = ', sice_ml(i,n,1)
       WRITE(*,'(A,F12.6)') 'sliq_ml(1) after        = ', sliq_ml(i,n,1)
       WRITE(*,'(A,F12.6)') 'snowdepth after         = ', snowdepth(i,n)
       WRITE(*,*) '------------------------------------------'
      ! stop
    END IF
    

   ELSE ! reset counter
      cold_shallow_puddle_hours(i,n) = 0.0
   END IF

      
!-----------------------------------------------------------------------------
! Add all rain to lake and reset fields 
!-----------------------------------------------------------------------------
   IF (exposed_water(i,n)) THEN

         IF (ls_rain(i) > 0.0 .OR. con_rain(i) > 0.0) THEN

         rain_add = ls_rain(i) + con_rain(i)
        
         lake_depth_ml(i,n) = lake_depth_ml(i,n)                                   &
              + (rain_add + snow_add) * timestep / rho_water
         
         !WRITE(*,*) '--- PRECIP ADDED TO LAKE/LID ---'
         !WRITE(*,*) 'i                  = ', i
         !WRITE(*,'(A,F16.8)') 'ls_rain       = ', ls_rain(i)
         !WRITE(*,'(A,F16.8)') 'con_rain      = ', con_rain(i)
         !WRITE(*,'(A,F16.8)') 'ls_snow       = ', ls_snow(i)
         !WRITE(*,'(A,F16.8)') 'con_snow      = ', con_snow(i)
         
         ls_rain(i)  = 0.0
         con_rain(i) = 0.0
        
      END IF
   END IF

!-----------------------------------------------------------------------------
! Get albedo - doing this again to output albedo to nc. this is already in
! src/science/radiation/jules_land_albedo_jls_mod.F90
! assume albedo is constant 0.05 for deep lakes (> 14m approx)
! that depth is only likely for idealised tests eqn 13
!-----------------------------------------------------------------------------
   flux_upper = 0.0
   flux_lower = 0.0
   flux_upper_diag = 0.0
   sw_absorb = 0.0
   dTdt = 0.0
   dhdt = 0.0
   dh_ice = 0.0
   dh_water = 0.0
   kdtdz_ml(i,n) = 0.0
   dhdt_lake_snow_ml(i,n) = 0.0
   
   IF (exposed_water(i,n)) THEN
     ! IF (lake_depth_ml(i)>=0.01) THEN
                     
      expon_term = 3.6 * lake_depth_ml(i,n)

      IF (expon_term < 50.0) THEN
         ex = EXP(-expon_term)
         lake_albedo_ml(i,n) = (9702.0 * ex + 1000.0)                 &
              / (-539.0 * ex + 20000.0)
      ELSE
         lake_albedo_ml(i,n) = 0.05
      END IF
               
   END IF

   
!-----------------------------------------------------------------------------
! Lower flux eqn 16
! lake water (lake_temp_ml)
!
!   ↓
!   │ flux_lower
!   │
!
! snowpack top (tm)
!-----------------------------------------------------------------------------
   IF (has_lake(i,n)) THEN

      delta_t = lake_temp_ml(i,n)  - tm
                 
      flux_lower = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
           * ABS(delta_t)**(4.0/3.0)

   END IF

            
!-----------------------------------------------------------------------------
! Upper flux, eqn 16
!
! air / tile surface (tstar_surft)  | lid bottom (tm)       
!                                   |
!   ↑                               |    !   ↑
!   │ flux_upper                    |    │ flux_upper
!   │                               |    |
!                                   |
! lake water (lake_temp_ml)         | lake water (lake_temp_ml)
!-----------------------------------------------------------------------------
   IF (exposed_water(i,n) .OR. has_vlid(i,n) .OR. has_lid(i,n)) THEN

      IF (exposed_water(i,n)) THEN
         delta_t = lake_temp_ml(i,n) - tstar_surft(i,n)
      ELSE
         delta_t = lake_temp_ml(i,n) - tm
      END IF
                 
      flux_upper = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
           * ABS(delta_t)**(4.0/3.0)

!-----------------------------------------------------------------------------
! Shortwave absorbed into lake, use a simple Beer-law form.
! Not using lake layers like in monarchs, can be simple if we only need a
! bulk lake temperature.                  
!-----------------------------------------------------------------------------

      IF (has_vlid(i,n)) THEN
         sw_absorb = sw_surft(i,n) * EXP(-tau_ice * vlid_depth_ml(i,n)) *   &
              (1.0 - EXP(-tau_water * lake_depth_ml(i,n)))
      ELSE IF (has_lid(i,n)) THEN
         sw_absorb = sw_surft(i,n) * EXP(-tau_ice * lid_depth_ml(i,n)) *    &
              (1.0 - EXP(-tau_water * lake_depth_ml(i,n)))
      ELSE IF (exposed_water(i,n)) THEN
         sw_absorb = sw_surft(i,n) * (1.0 - EXP(-tau_water * lake_depth_ml(i,n)))
      END IF

      
!-----------------------------------------------------------------------------
! Change in lake temperature. eqn 15
!-----------------------------------------------------------------------------
         
      dTdt = (- flux_upper - flux_lower + sw_absorb) /                  &
           (rho_water * hcapw * lake_depth_ml(i,n))
            
      lake_temp_ml(i,n) = lake_temp_ml(i,n) + timestep * dTdt

   
      IF (lake_temp_ml(i,n) > 273.15 + 27.0) THEN
         
         WRITE(*,*) '*** meltlake temp > 27C ***'
         WRITE(*,*) 'has_lid(i)                = ', has_lid(i,n)
         WRITE(*,*) 'has_lake(i)               = ', has_lake(i,n)
         WRITE(*,*) 'exposed_water(i)          = ', exposed_water(i,n)
         WRITE(*,'(A,F16.8)') 'sw_absorb       = ', sw_absorb
         WRITE(*,'(A,F16.8)') 'tstar_surft     = ', tstar_surft(i,n) - 273.15
         WRITE(*,'(A,F16.8)') 'lake_temp_ml    = ', lake_temp_ml(i,n) - 273.15
         WRITE(*,'(A,F16.8)') 'lake_depth_ml   = ', lake_depth_ml(i,n)
         WRITE(*,'(A,F16.8)') 'delta_t         = ', delta_t
         WRITE(*,'(A,F16.8)') 'flux_upper      = ', flux_upper
         WRITE(*,'(A,F16.8)') 'flux_lower      = ', flux_lower
         WRITE(*,'(A,F16.8)') 'sw_absorb       = ', sw_absorb
         WRITE(*,'(A,F16.8)') 'dTdt            = ', dTdt
         STOP 'meltlake temp exceeded 27C'
      END IF
      
   END IF ! exposed_water or has_lid
    
!-----------------------------------------------------------------------------
! Refreezing will trigger lid formation, see lid_evolve subroutine
!-----------------------------------------------------------------------------
   IF (lake_temp_ml(i,n) < tm) lake_temp_ml(i,n) = tm 


!-----------------------------------------------------------------------------
! Bottom Stefan condition is active when virtual or permanent lids are 
! present or there is exposed water
!-----------------------------------------------------------------------------         
   IF (has_lake(i,n)) THEN
      
      !IF (has_lake(i,n) .AND. (.NOT. has_lid(i,n)) .AND. (.NOT. has_vlid(i,n))) THEN ! stop if lid
!-----------------------------------------------------------------------------
! Re-calculate flux_upper before lake temp update for cross checking against 
! surf_ht_flux in sf_flux. flux_upper=surf_ht_flux in magnitude but have 
! opposite signs. They match all good :)
!-----------------------------------------------------------------------------
      
      IF (exposed_water(i,n)) THEN
         delta_t =  lake_temp_ml(i,n) - tstar_surft(i,n)

         flux_upper_diag = SIGN(1.0, delta_t) * rho_water * hcapw *  &
              Jturb * ABS(delta_t)**(4.0/3.0)
      ELSE
         flux_upper_diag = 0.0
      END IF

!-----------------------------------------------------------------------------
! heat conducted away from lake bottom into colder snowpack 
! eqn 14. direction format:  from the T1 to the T2
! if lake present on top of snowpack then use kdtdz as heat flux into
! snowpack top (surf_htf_surft)
!
! interface / snowpack top (tm)
!
!   ↓
!   │ kdtdz
!   │
!
! top snow layer (tsnow0_ml)
!-----------------------------------------------------------------------------                   
            
      kdtdz_ml(i,n) = ksnow0_ml(i,n) * (tm - tsnow0_ml(i,n)) / ds_ml(i,n,1)

      !IF (timestep_number >= 330 .AND. timestep_number <= 350) THEN

   !WRITE(*,*) '--- STEFAN FLUX INPUTS ---'
   !WRITE(*,'(A,I8)')     'timestep_number        = ', timestep_number
   !WRITE(*,'(A,I8)')     'i                      = ', i
   !WRITE(*,'(A,I8)')     'n                      = ', n

   !WRITE(*,'(A,F16.8)')  'ds_ml(1)               = ', ds_ml(i,n,1)
   !WRITE(*,'(A,F16.8)')  'sice_ml(1)             = ', sice_ml(i,n,1)
   !WRITE(*,'(A,F16.8)')  'sliq_ml(1)             = ', sliq_ml(i,n,1)
   !WRITE(*,'(A,F16.8)')  'tsnow0_ml (C)          = ', &
   !     tsnow0_ml(i,n) - tm
   !WRITE(*,'(A,F16.8)')  'ksnow0_ml              = ', ksnow0_ml(i,n)
   !WRITE(*,'(A,F16.8)')  'kdtdz_ml               = ', kdtdz_ml(i,n)
   !WRITE(*,'(A,F16.8)')  'flux_lower             = ', flux_lower
   !WRITE(*,'(A,F16.8)')  'flux_lower - kdtdz     = ', &
   !     flux_lower - kdtdz_ml(i,n)
   !WRITE(*,'(A,F16.8)')  'dhdt_lake_snow_ml      = ', &
   !     dhdt_lake_snow_ml(i,n)
   !IF (timestep_number == 350) STOP 'debug stop after timestep 350'

!END IF
!-----------------------------------------------------------------------------
! dhdt is an ice equivalent retreat rate eqn 14
! only allow Stefan melting if the lake supplies more heat to the
! interface than the firn conducts away below.
!-----------------------------------------------------------------------------
                       
      dhdt = 0.0
      IF (flux_lower > kdtdz_ml(i,n)) THEN
         dhdt = (flux_lower - kdtdz_ml(i,n)) / (rho_ice * lf)
      END IF

    
!-----------------------------------------------------------------------------
! dhdt (m of ice per sec) --> dh_water (m of water per timestep)
!-----------------------------------------------------------------------------
      
      dh_ice   = timestep * dhdt
      dh_water = dh_ice * rho_ice / rho_water
            
!-----------------------------------------------------------------------------
! output Stefan boundary as diagnostic
!-----------------------------------------------------------------------------
      dhdt_lake_snow_ml(i,n) = dh_ice
      
!-----------------------------------------------------------------------------
! Store the requested Stefan boundary retreat.
! The snow mass and geometry are adjusted later, after layersnow, so that
! sice, sliq and ds remain consistent with the JULES snow-layer geometry.
!-----------------------------------------------------------------------------
      IF (dh_ice > 0.0) THEN
         dhdt_lake_snow_ml(i,n) = dh_ice
      END IF ! Stefan dh_ice > 0

       
     ! IF (timestep_number < 340) THEN
     ! IF (timestep_number >= 333 .AND. timestep_number <= 335) THEN
         
      !      WRITE(*,*) '--- MELTLAKE_EVOLVE DEBUG ---'
      !      WRITE(*,'(A,I8)')    'timestep_number        = ', timestep_number
            !WRITE(*,'(A,I8)')    'i                      = ', i
       !     WRITE(*,'(A,L2)')    'has_lake               = ', has_lake(i,n)
       !     WRITE(*,'(A,L2)')    'exposed_water          = ', exposed_water(i,n)
       !     WRITE(*,'(A,L2)')    'has_vlid               = ', has_vlid(i,n)
       !     WRITE(*,'(A,L2)')    'has_lid                = ', has_lid(i,n)
            
       !     WRITE(*,'(A,F12.6)') 'lake_depth_ml          = ', lake_depth_ml(i,n)
       !     WRITE(*,'(A,F12.6)') 'dhdt_lake_snow_ml      = ', dhdt_lake_snow_ml(i,n)
       !     WRITE(*,'(A,F12.6)') 'lake_inflow (t-1)      = ', lake_inflow(i,n)
            
         !   WRITE(*,'(A,F12.6)') 'snow_surft             = ', snow_surft(i)
       !     WRITE(*,'(A,F12.6)') 'lake_temp_ml (C)       = ', lake_temp_ml(i) - 273.15

   
          !  WRITE(*,'(A,I8)')    'nsnow                 = ', nsnow(i)
          !  WRITE(*,'(A,F12.6)') 'ds_ml(1)              = ', ds_ml(i,1)
          !  WRITE(*,'(A,F12.6)') 'ds_ml(end)            = ', ds_ml(i,nsnow(i))
          !  WRITE(*,'(A,F12.6)') 'frac_melt             = ', frac_melt
          !  WRITE(*,'(A,F12.6)') 'snowdepth             = ', snowdepth(i)
          !  WRITE(*,'(A,F12.6)') 'sum_ds                = ', SUM(ds_ml(i,1:nsnow(i)))
          !  WRITE(*,'(A,F12.6)') 'difference            = ', snowdepth(i) - &
          !       SUM(ds_ml(i,1:nsnow(i)))

           ! WRITE(*,'(A,F12.6)') 'snow_surft            = ', snow_surft(i)
           ! WRITE(*,'(A,F12.6)') 'sum(sice+sliq)        = ', SUM(sice_ml(i,1:nsnow(i))) + SUM(sliq_ml(i,1:nsnow(i)))
           ! WRITE(*,'(A,F12.6)') 'sice_ml(1)            = ', sice_ml(i,1)
           ! WRITE(*,'(A,F12.6)') 'sliq_ml(1)            = ', sliq_ml(i,1)
        !    WRITE(*,'(A,F12.6)') 'tsnow_ml(1) (C)       = ', tsnow0_ml(i,n) - 273.15
        !    WRITE(*,'(A,F12.6)') 'ksnow0_ml             = ', ksnow0_ml(i)
        !    WRITE(*,'(A,F12.6)') 'kdtdz_ml              = ', kdtdz_ml(i)
   

         !   WRITE(*,'(A,F12.6)') 'flux_lower             = ', flux_lower
         !   WRITE(*,'(A,F12.6)') 'flux_upper             = ', flux_upper
         !   WRITE(*,'(A,F12.6)') 'sw_absorb              = ', sw_absorb
         !   WRITE(*,'(A,F12.6)') 'dTdt                   = ', dTdt
         !   WRITE(*,'(A,F12.6)') 'dhdt                   = ', dhdt
         !   WRITE(*,'(A,F12.6)') 'dh_ice                 = ', dh_ice
         !   WRITE(*,'(A,F12.6)') 'dh_water               = ', dh_water

            !IF (timestep_number == 335) STOP 'debug stop after timestep 340'

       !  END IF
     
           
   END IF !has_lake


   
!-----------------------------------------------------------------------------
! Re-calculate state after bottom Stefan melt
!-----------------------------------------------------------------------------
   CALL update_lake_states(i,n)

END DO ! land pts
END IF ! elev ice pts
END DO ! tile pts

!IF (timestep_number==2451) THEN
!   stop
!END IF

!$OMP END PARALLEL DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

SUBROUTINE update_lake_states(i,n)

   INTEGER, INTENT(IN) :: i, n

     
   has_lid(i,n)  = (lid_depth_ml(i,n)  > lid_min_depth)

   has_vlid(i,n) = (vlid_depth_ml(i,n) >= vlid_seed_depth) .AND.              &
              (vlid_depth_ml(i,n) <  lid_min_depth) .AND.                     &
              .NOT. has_lid(i,n)
   

   has_lake(i,n) = (lake_depth_ml(i,n) > lake_min_depth) .OR.                 &
                 ((lake_depth_ml(i,n) > 0.0) .AND.                            &
                  (has_lid(i,n) .OR. has_vlid(i,n)))

   exposed_water(i,n) = (lake_depth_ml(i,n) > lake_min_depth) .AND.           &
                      .NOT. has_lid(i,n) .AND. .NOT. has_vlid(i,n)

END SUBROUTINE update_lake_states


END SUBROUTINE calculate_lake_bottom_melt
END MODULE calculate_lake_bottom_melt_mod


         
