! *****************************COPYRIGHT*******************************

! (c) [University of Edinburgh]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE SNOW-----------------------------------------------------

! Description:
!     Calling routine for meltlake module
! Method:
!     Increment lake depth using excess water from snowpack
!     Get lake albedo for nc output
!     Get lake temp from tile surface temp using 4/3 law
!
! Code Owner: s.r.shannon@reading.ac.uk
!
! Subroutine Interface:
MODULE meltlake_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MELTLAKE_MOD'

CONTAINS

SUBROUTINE  meltlake (land_pts,                 & !IN
                      timestep,                 & !IN
                      nsurft,                   & !IN 
                      surft_pts,                & !IN
                      surft_index,              & !IN
                      lake_inflow,              & !IN
                      dt_elev_ml,               & !IN
                      lw_down_surft,            & !IN  
                      tstar_surft,              & !IN/OUT 
                      sw_surft,                 & !IN  already uses lake albedo
                      lake_depth_ml,            & !IN/OUT
                      lake_albedo_ml,           & !IN/OUT
                      lake_temp_ml,             & !IN/OUT
                      exposed_water,            & !IN/OUT
                      !Ancil info (IN)
                      l_lice_point,             & !IN (land_pts)
                      l_lice_surft)               !IN (ntype)  


!USE sf_diags_mod,            ONLY: strnewsfdiag
USE jules_surface_types_mod, ONLY: ntype
USE jules_surface_mod,       ONLY: l_elev_land_ice
USE theta_field_sizes,       ONLY: t_i_length, t_j_length

USE csigma,                  ONLY: sbcon

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

USE model_time_mod, ONLY: timestep_number

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
  timestep              ! Timestep length (s).


REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  sw_surft(land_pts,nsurft),                                                   &                                                   
    ! Net shortwave radiation on tile (W/m2) using lake albedo
  lw_down_surft(land_pts,nsurft),                                              &
    ! Surface downward LW radiation on tiles (W/m2), jules_land_sf_implicit.F90
  lake_inflow(land_pts,nsurft),                                                &
       ! melt water from snowpack (kgm-2)
  dt_elev_ml(land_pts,nsurft)
       ! orographic temp offset (elevated tile temp - gridbox mean temp (oK))  
     
!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_pts(nsurft),                                                           &
    ! Number of tile points.
  surft_index(land_pts,nsurft)
    ! Index of tile points.


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
   tstar_surft(land_pts,nsurft)
    ! Tile surface temperature (K)

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------

!ancil_info (IN)
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)
LOGICAL, INTENT(IN) :: l_lice_surft(ntype)

LOGICAL, INTENT(IN OUT) :: exposed_water(land_pts,nsurft)

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
  q_base,                                                                      &
    ! heat flux from lake into snowpack (W m⁻²)
  q_sw,                                                                        &
    ! Depth integrated shortwave absorbed in lake (W m-2) 
  q_lw,                                                                        &
    ! Net longwave flux at lake surface (W m-2)
  q_sens_lat,                                                                  &
    ! sensible + latent heat flux (W m-2)
  q_eff,                                                                       &
    ! Net energy flux into lake (W m-2)
  dh_dt,    &                   
  flux_lower, &
    ! flux lake and ice
  flux_upper, &
    ! flux between atmosphere and lake
  delta_t, &
  dt, &
  dTdt, &
  q_in, &
  sw_absorb, &
  lake_temp_eff_ml
    ! elevation adjust to lake temp to match sf_flux  

!REAL, DIMENSION(land_pts, nsurft), INTENT(IN) :: q_surface

REAL, PARAMETER :: ftl = 20.0   ! Sensible heat flux (W m-2), positive upward
REAL, PARAMETER :: le  = 50.0   ! Latent heat flux (W m-2), positive upward


REAL, PARAMETER :: tau   = 1.0
                 ! shortwave attenuation coeff (m-1)
                 ! 0.5 light penatrates deeper,  2-3 shallow heating

REAL, PARAMETER :: emis_water  = 0.98
REAL, PARAMETER :: hcon_snow  = 2.2 !thermal conductivity of ice (~2.2 W m⁻¹ K⁻¹)

REAL, PARAMETER :: dz_lake_base = 0.1 !assumed thickness (m) of the layer controlling heat exchange beneath the lake

REAL, PARAMETER :: tbase  = 273.15

REAL, PARAMETER :: Jturb  = 1.907e-5
                 ! Turbulent heat flux factor (ms⁻¹ K⁻¹/3) Eqn 16 Buzzard 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE'

! useful wildcard grep 
! grep -RInE 'dtstar_surft[[:space:]]*\(.*\)[[:space:]]*='

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO n = 1,nsurft
   DO k = 1,surft_pts(n)
      i = surft_index(k,n)
      
      IF (l_elev_land_ice .AND. l_lice_surft(n)) THEN
           
!-----------------------------------------------------------------------------
! Create an exposed water flag. Use this to modify snowpack upper boundary
!-----------------------------------------------------------------------------
               
         exposed_water(i,n) = (lake_depth_ml(i,n) > 0.1)
               
!-----------------------------------------------------------------------------
! Increment the lake depth
!-----------------------------------------------------------------------------
         IF (lake_inflow(i,n) > 0.0) THEN
            lake_depth_ml(i,n) = lake_depth_ml(i,n) +                          &
                 lake_inflow(i,n) / rho_water
         END IF
!-----------------------------------------------------------------------------
! Get albedo - doing this again to output albedo to nc. this is already in
! src/science/radiation/jules_land_albedo_jls_mod.F90
! assume albedo is constant 0.05 for deep lakes (> 14m approx)
! that depth is only likely for idealised tests eqn 13
!-----------------------------------------------------------------------------
         IF (exposed_water(i,n)) THEN
                     
            expon_term = 3.6 * lake_depth_ml(i,n)

            IF (expon_term < 50.0) THEN
               ex = EXP(-expon_term)
               lake_albedo_ml(i,n) = (9702.0 * ex + 1000.0)                 &
                    / (-539.0 * ex + 20000.0)
            ELSE
               lake_albedo_ml(i,n) = 0.05
            END IF
               
!-----------------------------------------------------------------------------
! Update lake temp. Assume the lake is turbulently mixed so use a bulk temp
!-----------------------------------------------------------------------------
            flux_upper = 0.0
            flux_lower = 0.0
            
!-----------------------------------------------------------------------------
! Upper flux, lake interior --> tile surface. eqn 16
! Use elevated tile effective lake temperature so this matches sf_flux
!-----------------------------------------------------------------------------
            lake_temp_eff_ml = lake_temp_ml(i,n) + dt_elev_ml(i,n)
            
            delta_t = lake_temp_eff_ml - tstar_surft(i,n) 
                 
            flux_upper = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
                 * ABS(delta_t)**(4.0/3.0)
            
!-----------------------------------------------------------------------------
! Lower flux, lake interior --> lens snowpack. eqn 16
!-----------------------------------------------------------------------------
                
            delta_t = lake_temp_ml(i,n)  - tm
                 
            flux_lower = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
                 * ABS(delta_t)**(4.0/3.0)

!-----------------------------------------------------------------------------
! Shortwave absorbed into lake, use a simple Beer-law form.
! Not using lake layers like in monarchs, can be simple if we only need a
! bulk lake temperature.                  
!-----------------------------------------------------------------------------
                 
            sw_absorb = sw_surft(i,n) *                                      &
                 (1.0 - EXP(-tau * lake_depth_ml(i,n)))

!-----------------------------------------------------------------------------
! Change in lake temperature. eqn 15
!-----------------------------------------------------------------------------
            dTdt = (- flux_upper - flux_lower + sw_absorb) /                &
                 (rho_water * hcapw * lake_depth_ml(i,n))

            lake_temp_ml(i,n) = lake_temp_ml(i,n) + timestep * dTdt
            
            IF (lake_temp_ml(i,n) < tm) lake_temp_ml(i,n) = tm ! refreezing 

            WRITE(*,*) '--- LAKE ENERGY DIAGNOSTICS ---'
            
            WRITE(*,'(A,F16.8)') 'sw_absorb                = ', sw_absorb
            WRITE(*,'(A,F16.8)') 'flux_upper               = ', flux_upper
            WRITE(*,'(A,F16.8)') 'flux_lower               = ', flux_lower
            WRITE(*,'(A,F16.8)') 'lake_temp                = ', lake_temp_ml(i,n) + dt_elev_ml(i,n) 
            WRITE(*,'(A,F16.8)') 'tstar_surft              = ', tstar_surft(i,n)
            WRITE(*,'(A,F16.8)') 'lake_inflow              = ', lake_inflow(i,n)
            WRITE(*,'(A,F16.8)') 'lake_temp - tstar_surft  = ', lake_temp_ml(i,n)-tstar_surft(i,n)
            WRITE(*,'(A,F16.8)') 'dTdt                     = ', dTdt
            
            WRITE(*,*) '--------------------------------', timestep_number
            !IF (timestep_number==2451) THEN
            !   stop
            !END IF
!-----------------------------------------------------------------------------
! Lake-snow boundary change from Stefan condition
!-----------------------------------------------------------------------------  
               
               !--- net SW absorbed on tile (sw_surft from calc_downward_rad_mod.F90)
               !--- uses albedo (over wrote snow with lake albedo)
              ! q_sw = sw_surft(i,n) *  &
              !(1.0 - EXP(-tau * lake_depth_ml(i,n))) / tau
                            

                 
              !---Stefan condition
              !q_base = hcon_snow * (lake_temp_ml(i,n) - t_base) /  dz_lake_base

              !delta_t = lake_temp_ml(i,n) - 273.15

              !flux_lower = SIGN(1.0, delta_t) * rho_water * hcapw * 1.907e-5 *  &
              !     ABS(delta_t)**(4.0/3.0)

              !flux_upper =  -(q_lw + q_sens_lat)

              ! --- net energy into the lake
              !dt = (-flux_lower - flux_upper + q_sw) /  &
              !     (rho_water * hcapw * lake_depth_ml(i,n))

              !lake_temp_ml(i,n) = lake_temp_ml(i,n) + timestep * dt
              
             ! lake_depth_ml(i,n) = lake_depth_ml(i,n) - dh_dt * timestep
              
             END IF ! exposed_water
          END IF ! elev land ice tile
       END DO
    END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE meltlake
END MODULE meltlake_mod
