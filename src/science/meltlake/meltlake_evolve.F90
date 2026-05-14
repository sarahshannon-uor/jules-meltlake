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
!     Stefan condition 
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
                      snow_surft,               & !IN/OUT
                      lake_depth_ml,            & !IN/OUT
                      lake_albedo_ml,           & !IN/OUT
                      lake_temp_ml,             & !IN/OUT
                      exposed_water,            & !IN/OUT
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


USE jules_snow_mod, ONLY:                                                      &
  snow_hcon
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
  timestep              ! Timestep length (s).


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
  ksnow0_ml(land_pts)

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
      
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                      &
     kdtdz_ml(land_pts)

LOGICAL, INTENT(IN OUT) :: exposed_water(land_pts)

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
  !kdTdz, &
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
REAL, SAVE :: cum_inflow_m = 0.0
REAL, SAVE :: cum_stefan_m = 0.0

REAL, PARAMETER :: tau   = 1.0
                 ! shortwave attenuation coeff (m-1)
                 ! 0.5 light penatrates deeper,  2-3 shallow heating

REAL, PARAMETER :: emis_water  = 0.98

REAL, PARAMETER :: Jturb  = 1.907e-5
                 ! Turbulent heat flux factor (ms⁻¹ K⁻¹/3) Eqn 16 Buzzard 


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE_EVOLVE'

INTEGER              :: errcode            ! Error code
CHARACTER(LEN=80) :: ERRMSG             ! Error message

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
   exposed_water(i) = (lake_depth_ml(i) > 0.1)
   
!-----------------------------------------------------------------------------
! Increment the lake depth
!-----------------------------------------------------------------------------
   IF (.NOT. exposed_water(i)) THEN
      IF (lake_inflow(i) > 0.0) THEN
         lake_depth_ml(i) = lake_depth_ml(i) +                          &
              lake_inflow(i) / rho_water
         
      END IF
   END IF

   cum_inflow_m = cum_inflow_m + lake_inflow(i) / rho_water

!-----------------------------------------------------------------------------
! Get albedo - doing this again to output albedo to nc. this is already in
! src/science/radiation/jules_land_albedo_jls_mod.F90
! assume albedo is constant 0.05 for deep lakes (> 14m approx)
! that depth is only likely for idealised tests eqn 13
!-----------------------------------------------------------------------------
      IF (exposed_water(i)) THEN
                     
         expon_term = 3.6 * lake_depth_ml(i)

         IF (expon_term < 50.0) THEN
            ex = EXP(-expon_term)
            lake_albedo_ml(i) = (9702.0 * ex + 1000.0)                 &
                 / (-539.0 * ex + 20000.0)
         ELSE
            lake_albedo_ml(i) = 0.05
         END IF
               
!-----------------------------------------------------------------------------
! Update lake temp. Assume the lake is turbulently mixed so use a bulk temp
!-----------------------------------------------------------------------------

         flux_upper = 0.0
         flux_lower = 0.0
            
!-----------------------------------------------------------------------------
! Upper flux, lake interior --> tile surface. eqn 16
!-----------------------------------------------------------------------------
                        
         delta_t = lake_temp_ml(i) - tstar_surft(i) 
                 
         flux_upper = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
              * ABS(delta_t)**(4.0/3.0)
            
!-----------------------------------------------------------------------------
! Lower flux, lake interior --> lens snowpack. eqn 16
!-----------------------------------------------------------------------------
                
         delta_t = lake_temp_ml(i)  - tm
                 
         flux_lower = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
              * ABS(delta_t)**(4.0/3.0)

!-----------------------------------------------------------------------------
! Shortwave absorbed into lake, use a simple Beer-law form.
! Not using lake layers like in monarchs, can be simple if we only need a
! bulk lake temperature.                  
!-----------------------------------------------------------------------------
                 
         sw_absorb = sw_surft(i) *                                        &
              (1.0 - EXP(-tau * lake_depth_ml(i)))

!-----------------------------------------------------------------------------
! Change in lake temperature. eqn 15
!-----------------------------------------------------------------------------
         
         dTdt = (- flux_upper - flux_lower + sw_absorb) /                 &
              (rho_water * hcapw * lake_depth_ml(i))
            
         lake_temp_ml(i) = lake_temp_ml(i) + timestep * dTdt
            
!-----------------------------------------------------------------------------
! Refreezing will trigger lid formation, return to this 
!-----------------------------------------------------------------------------
         IF (lake_temp_ml(i) < tm) lake_temp_ml(i) = tm 

!-----------------------------------------------------------------------------
! Re-calculate flux_upper after lake temp update for cross checking against 
! surf_ht_flux in sf_flux. the values should match but have opposite signs
! They match all good :)
!-----------------------------------------------------------------------------         
         
         delta_t =  lake_temp_ml(i) - tstar_surft(i)

         flux_upper_diag = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb *          &
              ABS(delta_t)**(4.0/3.0)

              
!-----------------------------------------------------------------------------
! heat conducted away from lake bottom into colder firn below
! use const value for thermal conduct for now (snow_hcon not ksnow0_ml(i))
! eqn 14 
!-----------------------------------------------------------------------------                   
            
            kdtdz_ml(i) = snow_hcon * (tm - tsnow_ml(i)) / ds_ml(i,1)

!-----------------------------------------------------------------------------
! dhdt is an ice equivalent retreat rate eqn 14
! only allow Stefan melting if the lake supplies more heat to the
! interface than the firn conducts away below.
!-----------------------------------------------------------------------------
                       
            dhdt = 0.0
            IF (flux_lower - kdtdz_ml(i) > 0.0) THEN
               dhdt = (flux_lower - kdtdz_ml(i)) / (rho_ice * lf)
            END IF

!-----------------------------------------------------------------------------
! dhdt (m of ice per sec) --> dh_water (m of water per timestep)
!-----------------------------------------------------------------------------

            dh_ice   = timestep * dhdt
            dh_water = dh_ice * rho_ice / rho_water
            
!-----------------------------------------------------------------------------
! add water from boundary retreat to the lake.
! remove the equivalent sice from the snow mass 
!-----------------------------------------------------------------------------           
           IF (dh_ice > 0.0) THEN
              
              dh_remain = dh_ice

!-----------------------------------------------------------------------------
! looping over snowpack levels and removing sice.  This might be overkill since 
! Stefan retreat will be small compared to ds(1) i.e. 0.1 m, but guarding against 
! a setup with smaller top layer depth . Removal is based on fraction
! stefan retreat/ds layer
! 
!-----------------------------------------------------------------------------                
              DO n = 1, nsnow(i)

                 IF (dh_remain <= 0.0) EXIT
                
                 ds_old = ds_ml(i,n)

                 frac_melt = MIN(1.0, dh_remain / ds_old)
                      
                 dsice = frac_melt * sice_ml(i,n)
                 dsliq = sliq_ml(i,n) ! this drains all sliq 

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

                                  
                 WRITE(*,*) 'n                                  = ', n
                 WRITE(*,'(A,F16.8)') 'frac_melt                = ', frac_melt
                 WRITE(*,'(A,F16.8)') 'dsice                    = ', dsice
                 WRITE(*,'(A,F16.8)') 'dsliq                    = ', dsliq
                 WRITE(*,'(A,F16.8)') 'sice                     = ', sice_ml(i,n)
                 WRITE(*,'(A,F16.8)') 'sliq                     = ', sliq_ml(i,n)
                 WRITE(*,'(A,F16.8)') 'ds_old                   = ', ds_old
                 WRITE(*,'(A,F16.8)') 'ds_ml(i,n)               = ', ds_ml(i,n) 
                 WRITE(*,'(A,F16.8)') 'lake_add                 = ', (dsice + dsliq) / rho_water
                 
              END DO !nsnow
              
             
              
           END IF ! Stefan dh_ice > 0


           WRITE(*,*) '--- LAKE ENERGY DIAGNOSTICS ---'
           WRITE(*,*) '--------------------------------', timestep_number
           WRITE(*,*) 'nsnow(i)', nsnow(i) 
           WRITE(*,'(A,F16.8)') 'sw_absorb                = ', sw_absorb
           WRITE(*,'(A,F16.8)') 'flux_upper_diag          = ', flux_upper_diag
           WRITE(*,'(A,F16.8)') 'flux_lower               = ', flux_lower
           WRITE(*,'(A,F16.8)') 'lake_temp after          = ', lake_temp_ml(i) -273.15 
           WRITE(*,'(A,F16.8)') 'tstar_surft              = ', tstar_surft(i) - 273.15
           WRITE(*,'(A,F16.8)') 'tstar_surft              = ', tstar_surft(i) + - 273.15
           WRITE(*,'(A,F16.8)') 'lake_inflow              = ', lake_inflow(i)
           WRITE(*,'(A,F16.8)') 'lake_temp - tstar_surft  = ', lake_temp_ml(i)-tstar_surft(i)
           WRITE(*,'(A,F16.8)') 'dTdt                     = ', dTdt
         
            
        
           WRITE(*,'(A,F16.8)') 'snow_hcon                = ', snow_hcon
           WRITE(*,'(A,F16.8)') 'ksnow0_ml(i)             =',  ksnow0_ml(i)
           WRITE(*,'(A,F16.8)') 'tsnow(i,1)              = ', tsnow_ml(i)
           ! WRITE(*,'(A,F16.8)') 'dm                      = ', dm
           WRITE(*,'(A,F16.8)') 'kdtdz(i)_ml              = ', kdtdz_ml(i)
           WRITE(*,'(A,F16.8)') 'flux_lower               = ', flux_lower
           WRITE(*,'(A,F16.8)') 'flux_lower - kdTdz       = ', flux_lower - kdtdz_ml(i)
           !WRITE(*,'(A,F16.8)')  'ds_ml(i),               = ', ds_ml(i,1)
           WRITE(*,'(A,F16.8)') 'dh_water (m)             = ', dh_water
           WRITE(*,'(A,F16.8)') 'snow_surft(i)            = ', snow_surft(i)
           
           WRITE(*,*) '------------------------------------------------'


          
           cum_stefan_m = cum_stefan_m + dh_water
           
           
           
        END IF ! exposed_water
        
        WRITE(*,'(A,F16.8)') 'cum_inflow_m', cum_inflow_m
        WRITE(*,'(A,F16.8)') 'cum_stefan_m', cum_stefan_m
   END DO ! land_pts

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE meltlake_evolve
END MODULE meltlake_evolve_mod

!IF (timestep_number==2451) THEN
            !   stop
         !END IF
         
