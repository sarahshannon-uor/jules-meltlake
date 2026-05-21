! ****************************rCOPYRIGHT*******************************

! (c) [University of Edinburgh]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE RELAYERSNOW-----------------------------------------------
! Description:
!     Calling routine for meltlake module
! Method:
!     Increment lake depth using excess water from snowpack
!     Get lake albedo for nc output
!     Get lake temp from tile surface temp using 4/3 law
!     Find boundary change between lake bottom and snowpack top using
!     Stefan condition.
!     Use Stefan condition to reduce the sice at the stop of the snowpack
!     Add the extra water from that sice to the lake depth.
!     When sice is reduced from the Stefan condition then also drain sliq
!     into the lake. Output adjusted sice, sliq and ds 
! Code Owner: s.r.shannon@reading.ac.uk
! Subroutine Interface:
MODULE meltlake_evolve_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MELTLAKE_EVOLVE_MOD'

CONTAINS

  SUBROUTINE meltlake_evolve(land_pts,          & !IN
                      timestep,                 & !IN
                      surft_pts,                & !IN
                      surft_index,              & !IN
                      nsnow,                    & !IN 
                      lake_inflow,              & !IN
                      ksnow0_ml,                & !IN
                      tsnow_ml,                 & !IN
                      lw_down_surft,            & !IN  
                      tstar_surft,              & !IN
                      sw_surft,                 & !IN  already uses lake albedo
                      lid_depth_ml,             & !IN
                      snow_surft,               & !IN/OUT
                      lake_depth_ml,            & !IN/OUT
                      lake_albedo_ml,           & !IN/OUT
                      lake_temp_ml,             & !IN/OUT
                      has_lake,                 & !IN/OUT
                      exposed_water,            & !IN/OUT
                      has_lid,                  & !IN/OUT
                      sice_ml,                  & !IN/OUT
                      sliq_ml,                  & !IN/OUT
                      ds_ml,                    & !IN/OUT 
                      kdtdz_ml)                   !OUT
                      
                     
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
   surft_pts
    ! Number of tile points.
  
    
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
   timestep
    ! Timestep length (s).


REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  sw_surft(land_pts),                                                          &                                                   
    ! Net shortwave radiation on tile (W/m2) using lake albedo
  lw_down_surft(land_pts),                                                     &
    ! Surface downward LW radiation on tiles (W/m2), jules_land_sf_implicit.F90
  tstar_surft(land_pts),                                                       &
    ! Tile surface temperature (K)
  lake_inflow(land_pts),                                                       &
       ! melt water from snowpack (kgm-2)
  tsnow_ml(land_pts),                                                          &
       ! snowpack top level temp (K)  
  ksnow0_ml(land_pts),                                                         &
       ! thermal conductivity of top snow layer (W m⁻¹ K⁻¹)
  lid_depth_ml(land_pts) 
       ! depth of frozen lid (m)

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_index(land_pts),                                                       &
    ! Index of tile points.
  nsnow(land_pts)                                                        
    ! Number of snow layers.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
   lake_depth_ml(land_pts),                                                    &
    ! Convective rainfall rate (kg/m2/s).
   lake_albedo_ml(land_pts),                                                   &
    ! Albedo of lake
   lake_temp_ml(land_pts),                                                     &
    ! Temperature of melt lake (K)
   snow_surft(land_pts),                                                       &
    ! Snow mass on tiles (kgm-2)
   sice_ml(land_pts,nsmax_ml),                                                 &
    ! Ice content of snow layers (kg/m2)
   sliq_ml(land_pts,nsmax_ml),                                                 &
    ! Liquid content of snow layers (kg/m2)
   ds_ml(land_pts, nsmax_ml)
      ! snowpack top level depth (m)
      
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
   kdtdz_ml(land_pts)
     ! Conductive heat flux from lake into snowpack. pass this to snow module
     ! so the heat flux will be adjusted if there is exposed water 

LOGICAL, INTENT(IN OUT) ::                                                     &
     has_lake(land_pts),                                                       &
     exposed_water(land_pts),                                                  &
     has_lid(land_pts)


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
  n                                                                           
    ! Tile loop counter.
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
  ds_old

! --- diagnostic to check how much lake depth is from inflow versus stefan
REAL, SAVE       :: cum_inflow_m = 0.0
REAL, SAVE       :: cum_stefan_m = 0.0

REAL, PARAMETER :: tau   = 1.0
                 ! shortwave attenuation coeff (m-1)
                 ! 0.5 light penatrates deeper,  2-3 shallow heating

REAL, PARAMETER :: emis_water  = 0.98

REAL, PARAMETER :: Jturb  = 1.907e-5
                 ! Turbulent heat flux factor (ms⁻¹ K⁻¹/3) Eqn 16 Buzzard 

REAL(KIND=real_jlslsm), PARAMETER :: lid_min_depth = 0.1
REAL(KIND=real_jlslsm), PARAMETER :: lake_min_depth = 0.1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE_EVOLVE'

INTEGER              :: errcode            ! Error code
CHARACTER(LEN=80)    :: ERRMSG             ! Error message

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! 
!-----------------------------------------------------------------------------
DO k = 1,surft_pts
   i = surft_index(k)

!-----------------------------------------------------------------------------
! Create an exposed water flag. Use this to modify snowpack upper boundary
!-----------------------------------------------------------------------------
   has_lake(i) = (lake_depth_ml(i) > lake_min_depth)
   has_lid(i)  = has_lake(i) .AND. (lid_depth_ml(i) > lid_min_depth)
   exposed_water(i) = has_lake(i) .AND. .NOT. has_lid(i)
   
!-----------------------------------------------------------------------------
! No lake yet grow one using meltwater 
!-----------------------------------------------------------------------------
   IF (.NOT. has_lake(i)) THEN
      IF (lake_inflow(i) > 0.0) THEN
         lake_depth_ml(i) = lake_depth_ml(i) +                          &
              lake_inflow(i) / rho_water
         print *, 'has_lake(i), has_lid(i), exposed_water(i)',has_lake(i), has_lid(i), exposed_water(i)
      END IF
   END IF

   cum_inflow_m = cum_inflow_m + lake_inflow(i) / rho_water

!-----------------------------------------------------------------------------
! Re-calculate flags after lake growth
!-----------------------------------------------------------------------------
   has_lake(i) = (lake_depth_ml(i) > lake_min_depth)
   has_lid(i)  = has_lake(i) .AND. (lid_depth_ml(i) > lid_min_depth)
   exposed_water(i) = has_lake(i) .AND. .NOT. has_lid(i)

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
     
   IF (exposed_water(i)) THEN
                     
      expon_term = 3.6 * lake_depth_ml(i)

      IF (expon_term < 50.0) THEN
         ex = EXP(-expon_term)
         lake_albedo_ml(i) = (9702.0 * ex + 1000.0)                 &
              / (-539.0 * ex + 20000.0)
      ELSE
         lake_albedo_ml(i) = 0.05
      END IF
               
   END IF

!-----------------------------------------------------------------------------
! Lower flux from lake --> snowpack eqn 16
!-----------------------------------------------------------------------------
   IF (has_lake(i)) THEN

      delta_t = lake_temp_ml(i)  - tm
                 
      flux_lower = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
           * ABS(delta_t)**(4.0/3.0)

   END IF

!-----------------------------------------------------------------------------
! Update lake temp. Assume the lake is turbulently mixed so use a bulk temp
!-----------------------------------------------------------------------------
            
!-----------------------------------------------------------------------------
! Upper flux, lake interior --> tile surface or lake interior to lid base eqn 16.
!-----------------------------------------------------------------------------
   IF (exposed_water(i) .OR. has_lid(i)) THEN

      IF (exposed_water(i)) THEN
         delta_t = lake_temp_ml(i) - tstar_surft(i)
      ELSE
         delta_t = lake_temp_ml(i) - tm
      END IF
                 
      flux_upper = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
           * ABS(delta_t)**(4.0/3.0)

!-----------------------------------------------------------------------------
! Shortwave absorbed into lake, use a simple Beer-law form.
! Not using lake layers like in monarchs, can be simple if we only need a
! bulk lake temperature.                  
!-----------------------------------------------------------------------------
      IF (exposed_water(i)) THEN
         sw_absorb = sw_surft(i) *                                        &
              (1.0 - EXP(-tau * lake_depth_ml(i)))
      ELSE
         sw_absorb = 0.0
      END IF

!-----------------------------------------------------------------------------
! Change in lake temperature. eqn 15
!-----------------------------------------------------------------------------
         
      dTdt = (- flux_upper - flux_lower + sw_absorb) /                 &
           (rho_water * hcapw * lake_depth_ml(i))
            
      lake_temp_ml(i) = lake_temp_ml(i) + timestep * dTdt

!-----------------------------------------------------------------------------
! Under a lid, keep bulk lake at melting point
!-----------------------------------------------------------------------------
      IF (has_lid(i)) THEN
         lake_temp_ml(i) = tm
      END IF
      
      IF (lake_temp_ml(i) > 273.15 + 27.0) THEN
         
         WRITE(*,*) '*** meltlake temp > 27C ***'
         WRITE(*,*) 'has_lid(i)                = ', has_lid(i)
         WRITE(*,*) 'has_lake(i)               = ', has_lake(i)
         WRITE(*,*) 'expsoed_water(i)          = ', exposed_water(i)
         WRITE(*,'(A,F16.8)') 'sw_absorb       = ', sw_absorb
         WRITE(*,'(A,F16.8)') 'tstar_surft     = ', tstar_surft(i) - 273.15
         WRITE(*,'(A,F16.8)') 'lake_temp_ml    = ', lake_temp_ml(i) - 273.15
         WRITE(*,'(A,F16.8)') 'lake_depth_ml   = ', lake_depth_ml(i)
         WRITE(*,'(A,F16.8)') 'delta_t         = ', delta_t
         WRITE(*,'(A,F16.8)') 'flux_upper      = ', flux_upper
         WRITE(*,'(A,F16.8)') 'flux_lower      = ', flux_lower
         WRITE(*,'(A,F16.8)') 'sw_absorb       = ', sw_absorb
         WRITE(*,'(A,F16.8)') 'dTdt            = ', dTdt
         STOP 'meltlake temp exceeded 27C'
      END IF
      
   END IF ! exposed_water or has_lid
    
!-----------------------------------------------------------------------------
! Refreezing will trigger lid formation, dont cap lake_temp but refer to
! lid_evolve subroutine
!-----------------------------------------------------------------------------
         IF (lake_temp_ml(i) < tm) lake_temp_ml(i) = tm 


!-----------------------------------------------------------------------------
! Bottom Stefan condition is active for lid and expsosed water
!-----------------------------------------------------------------------------         
   IF (has_lake(i)) THEN

!-----------------------------------------------------------------------------
! Re-calculate flux_upper before lake temp update for cross checking against 
! surf_ht_flux in sf_flux. flux_upper=surf_ht_flux in magnitude but have 
! opposite signs. They match all good :)
!-----------------------------------------------------------------------------
      
      IF (exposed_water(i)) THEN
         delta_t =  lake_temp_ml(i) - tstar_surft(i)

         flux_upper_diag = SIGN(1.0, delta_t) * rho_water * hcapw *  &
              Jturb * ABS(delta_t)**(4.0/3.0)
      ELSE
         flux_upper_diag = 0.0
      END IF

!-----------------------------------------------------------------------------
! heat conducted away from lake bottom into colder snowpack 
!  eqn 14. direction format:  from the T1 to the T2  
!-----------------------------------------------------------------------------                   
            
      kdtdz_ml(i) = ksnow0_ml(i) * (tm - tsnow_ml(i)) / ds_ml(i,1)

!-----------------------------------------------------------------------------
! dhdt is an ice equivalent retreat rate eqn 14
! only allow Stefan melting if the lake supplies more heat to the
! interface than the firn conducts away below.
!-----------------------------------------------------------------------------
                       
      dhdt = 0.0
      IF (flux_lower > kdtdz_ml(i)) THEN
         dhdt = (flux_lower - kdtdz_ml(i)) / (rho_ice * lf)
      END IF

!-----------------------------------------------------------------------------
! dhdt (m of ice per sec) --> dh_water (m of water per timestep)
!-----------------------------------------------------------------------------

      dh_ice   = timestep * dhdt
      dh_water = dh_ice * rho_ice / rho_water
            
!-----------------------------------------------------------------------------
! add water from boundary retreat to the lake.
! remove the equivalent sice from the snow mass. This is working if there is a lid 
!-----------------------------------------------------------------------------           
      !IF (dh_ice > 0.0) THEN
     
         IF (dh_ice > 0.0 .AND. lid_depth_ml(i) < 0.001) THEN ! to deactive if there is a lid

            dh_remain = dh_ice

!-----------------------------------------------------------------------------
! looping over snowpack levels and removing sice.  This might be overkill since 
! Stefan retreat will be small compared to ds(1) i.e. 0.1 m, but guarding against 
! a setup with smaller top layer depth . Removal is based on fraction of top 
! snow layer retreating i.e. stefan retreat/ds top snow layer
! remove ice mass fraction, and liq mass fraction
!-----------------------------------------------------------------------------                
         DO n = 1, nsnow(i)

            IF (dh_remain <= 0.0) EXIT
                
            ds_old = ds_ml(i,n)

            frac_melt = MIN(1.0, dh_remain / ds_old)
                      
            dsice = frac_melt * sice_ml(i,n)
            dsliq = frac_melt * sliq_ml(i,n)

!-----------------------------------------------------------------------------
! water from ice retreat + draining all sliq from layer
!-----------------------------------------------------------------------------
            lake_depth_ml(i) = lake_depth_ml(i) + (dsice + dsliq) / rho_water

!-----------------------------------------------------------------------------
! adjust sice, sliq and ds 
!-----------------------------------------------------------------------------                 
            sice_ml(i,n) = sice_ml(i,n) - dsice
            sliq_ml(i,n) = sliq_ml(i,n) - dsliq
            ds_ml(i,n)   = ds_old  * (1.0 - frac_melt)
                 
            dh_remain = dh_remain - frac_melt * ds_old

                                  
            ! WRITE(*,*) 'n                                  = ', n
            ! WRITE(*,'(A,F16.8)') 'frac_melt                = ', frac_melt
            ! WRITE(*,'(A,F16.8)') 'dsice                    = ', dsice
            ! WRITE(*,'(A,F16.8)') 'dsliq                    = ', dsliq
            ! WRITE(*,'(A,F16.8)') 'sice                     = ', sice_ml(i,n)
            ! WRITE(*,'(A,F16.8)') 'sliq                     = ', sliq_ml(i,n)
            ! WRITE(*,'(A,F16.8)') 'ds_old                   = ', ds_old
            ! WRITE(*,'(A,F16.8)') 'ds_ml(i,n)               = ', ds_ml(i,n) 
            ! WRITE(*,'(A,F16.8)') 'lake_add                 = ', (dsice + dsliq) / rho_water
                 
         END DO !nsnow

         
      END IF ! Stefan dh_ice > 0


      !WRITE(*,*) '--- LAKE ENERGY DIAGNOSTICS ---'
      !WRITE(*,*) '--------------------------------', timestep_number
      !WRITE(*,*) 'nsnow(i)                           = ', nsnow(i) 
      !WRITE(*,'(A,F16.8)') 'sw_absorb                = ', sw_absorb
       WRITE(*,'(A,F16.8)') 'flux_upper_diag          = ', flux_upper_diag
      !WRITE(*,'(A,F16.8)') 'flux_lower               = ', flux_lower
      !WRITE(*,'(A,F16.8)') 'lake_temp after          = ', lake_temp_ml(i) -273.15 
      !WRITE(*,'(A,F16.8)') 'tstar_surft              = ', tstar_surft(i) - 273.15
      !WRITE(*,'(A,F16.8)') 'lake_inflow              = ', lake_inflow(i)
      !WRITE(*,'(A,F16.8)') 'lake_temp - tstar_surft  = ', lake_temp_ml(i)-tstar_surft(i)
      !WRITE(*,'(A,F16.8)') 'dTdt                     = ', dTdt
         
                         
      !WRITE(*,'(A,F16.8)') 'ksnow0_ml(i)             =',  ksnow0_ml(i)
      !WRITE(*,'(A,F16.8)') 'tsnow(i,1)              = ', tsnow_ml(i)
      !WRITE(*,'(A,F16.8)') 'kdtdz(i)_ml              = ', kdtdz_ml(i)
      !WRITE(*,'(A,F16.8)') 'flux_lower               = ', flux_lower
      !WRITE(*,'(A,F16.8)') 'flux_lower - kdTdz       = ', flux_lower - kdtdz_ml(i)
      !WRITE(*,'(A,F16.8)')  'ds_ml(i),               = ', ds_ml(i,1)
      !WRITE(*,'(A,F16.8)') 'dh_water (m)             = ', dh_water
      !WRITE(*,'(A,F16.8)') 'snow_surft(i)            = ', snow_surft(i)
           
      !WRITE(*,*) '------------------------------------------------'
      
      
      cum_stefan_m = cum_stefan_m + dh_water
           
   END IF !has_lake

!-----------------------------------------------------------------------------
! Re-calculate flags after bottom Stefan melt
!-----------------------------------------------------------------------------
   has_lake(i) = (lake_depth_ml(i) > lake_min_depth)
   has_lid(i)  = has_lake(i) .AND. (lid_depth_ml(i) > lid_min_depth)
   exposed_water(i) = has_lake(i) .AND. .NOT. has_lid(i)

     
    WRITE(*,'(A,F16.8)') 'cum_inflow_m', cum_inflow_m
    WRITE(*,'(A,F16.8)') 'cum_stefan_m', cum_stefan_m

END DO ! land_pts

!IF (timestep_number==2451) THEN
!   stop
!END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE meltlake_evolve
END MODULE meltlake_evolve_mod


         
