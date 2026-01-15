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

! Subroutine Interface:
MODULE meltlake_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MELTLAKE_MOD'

CONTAINS

SUBROUTINE  meltlake (land_pts,                 & !IN
                      timestep,                 & !IN
                      nsurft,                   & !IN 
                      surft_pts,                & !IN
                      surft_index,              & !IN
                      snice_runoff_surft,       & !IN/OUT 
                      !lw_down_surft,            & !IN  can't highjack these without using switches (set in sf_diag.F90) 
                      !lw_up_surft,              & !IN     
                      !tstar_surft,              & !IN/OUT
                      !sw_surft,                 & !IN  
                      lake_depth_ml,            & !IN/OUT
                      lake_albedo_ml,           & !IN/OUT
                      !Ancil info (IN)
                      l_lice_point,             & !IN (land_pts)
                      l_lice_surft)               !IN (ntype)  


!USE sf_diags_mod,            ONLY: strnewsfdiag
USE jules_surface_types_mod, ONLY: ntype
USE jules_surface_mod,       ONLY: l_elev_land_ice
USE water_constants_mod,     ONLY: rho_water 
USE csigma,                  ONLY: sbcon
USE theta_field_sizes,       ONLY: t_i_length, t_j_length
  
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts              ! Total number of land points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep              ! Timestep length (s).

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  nsurft,                                                                      &
    ! Number of land tiles.
  surft_pts(nsurft),                                                           &
    ! Number of tile points.
  surft_index(land_pts,nsurft)!,                                                &
    ! Index of tile points.
  !sw_surft(land_pts,nsurft),                                                   &
    ! Net shortwave radiation on tile (W/m2)
    ! lake albedo was applied in calc_downward_radiation_mod.F90
  !lw_down_surft(land_pts,nsurft),                                              &                                                
    ! Surface downward LW radiation on tiles (W/m2), jules_land_sf_implicit.F90
  !lw_up_surft(land_pts,nsurft)                                              
    ! Surface upward LW radiation on tiles (W/m2), jules_land_sf_implicit.F90 
!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
!TYPE (strnewsfdiag), INTENT(IN OUT) :: sf_diag

!INTEGER, INTENT(IN OUT) ::                                                     &
!  nsnow(land_pts,nsurft)   ! Number of snow layers.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
   snice_runoff_surft(land_pts,nsurft),                                        &
    ! Excess water out of snowpack kg/m/s
   lake_depth_ml(land_pts,nsurft),                                             & 
    ! Convective rainfall rate (kg/m2/s).
   lake_albedo_ml(land_pts,nsurft)!,                                            &
    ! Albedo of lake
   !tstar_surft(land_pts,nsurft)

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
!REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
!  ds_ml(land_pts,nsurft,nsmax_ml),                                             &
    ! Snow layer thicknesses (m).
  
!-----------------------------------------------------------------------------
! New arguments to replace USE statements
!-----------------------------------------------------------------------------

!ancil_info (IN)
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)
LOGICAL, INTENT(IN) :: l_lice_surft(ntype)

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
!REAL(KIND=real_jlslsm) ::                                                      &
!  csnow(land_pts,nsmax_ml),                                                       &
    ! Areal heat capacity of layers (J/K/m2).
  
! Snow quantities on a single surft, i.e. over snow layers (sl)
!REAL(KIND=real_jlslsm) ::                                                      &
!  ds_sl(land_pts,nsmax_ml),                                                       &
    ! Snow layer thicknesses (m).
  
REAL(KIND=real_jlslsm) ::                                                      &
  sw_net_lake,                                                                 &
    ! net shortwave radiation on meltlake (W/m2)
  lw_net_lake
    ! net longwave radiation on meltlake (W/m2)

REAL, PARAMETER :: rho_w = 1000.0
REAL, PARAMETER :: cp_w  = 4186.0
REAL, PARAMETER :: cp_a  = 1004.0
REAL, PARAMETER :: Rd    = 287.05
!REAL, PARAMETER :: sbcon = 5.670374419e-8
REAL, PARAMETER :: Lv    = 2.5e6
REAL, PARAMETER :: Lf    = 3.34e5
REAL, PARAMETER :: alpha_w = 0.07
REAL, PARAMETER :: emis_w  = 0.97
REAL, PARAMETER :: CH = 1.5e-3
REAL, PARAMETER :: CE = 1.5e-3
REAL, PARAMETER :: Tm = 273.15


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE'

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO n = 1,nsurft
    DO k = 1,surft_pts(n)
        i = surft_index(k,n)
          
            IF (l_elev_land_ice .AND. l_lice_surft(n)) THEN

              IF (snice_runoff_surft(i,n) > 0.0) THEN

               lake_depth_ml(i,n) = lake_depth_ml(i,n) + &
              (snice_runoff_surft(i,n)  * timestep) / rho_water

              lake_albedo_ml(i,n) = (9702.0 + 1000.0 * EXP(3.6 * lake_depth_ml(i,n))) &
             / (-539.0 + 20000.0 * EXP(3.6 * lake_depth_ml(i,n)))

            
            !print *, 'sw_surft(i,n)',sw_surft(i,n)   
			!emis_surft(i,n) from physiol_jls.F90 from nveg_params
            !lake_h = lake_depth_ml(i,n)
            !lake_temp = lake_temp_ml(i,n) 
            !fluxes%sw_surft
              
            !lw_net_lake = emis_w * (lw_down - sbcon * tstar_surft(i,n)**4)
            

            !H   = rho_a * cp_a * CH * U * (Tw - Ta)
            !LE  = rho_a * Lv   * CE * U * (qsat - qa)

            !lake_temp_ml(i,n) = lake_temp + timestep * (sw_net_lake + lw_net_lake - H - LE) / (rho_w * cp_w * lake_h)


            !if (i.eq.1) then
            !print *, 'runoff,depth,albedo ' , sf_diag%snice_runoff_surft(i,n), lake_depth_ml(i,n), lake_albedo_ml(i,n)  
            !end if 
            END IF 
        END IF 
    END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE meltlake
END MODULE meltlake_mod
