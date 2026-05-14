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
                      nsnow,                    & !IN
                      lake_inflow,              & !IN
                      kdtdz_ml,                 & !OUT
                      lw_down_surft,            & !IN  
                      tstar_surft,              & !IN/OUT 
                      sw_surft,                 & !IN  already adjusted for exposed water albedo
                      tsnow_ml,                 & !IN/OUT
                      ds_ml,                    & !IN/OUT
                      sice_ml,                  & !IN/OUT
                      sliq_ml,                  & !IN/OUT
                      snow_surft,               & !IN/OUT
                      lake_depth_ml,            & !IN/OUT
                      lake_albedo_ml,           & !IN/OUT
                      lake_temp_ml,             & !IN/OUT
                      exposed_water,            & !IN/OUT
                      ksnow0_ml,                &! IN
                      !Ancil info (IN)
                      l_lice_point,             & !IN (land_pts)
                      l_lice_surft)               !IN (ntype)  

USE meltlake_evolve_mod, ONLY: meltlake_evolve
USE relayersnow_mod, ONLY: relayersnow
  
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

USE jules_meltlake_mod, ONLY: nsmax_ml

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
  !dt_elev_ml(land_pts,nsurft),                                                 &
       ! orographic temp offset (elevated tile temp - gridbox mean temp (oK))  
  ksnow0_ml(land_pts,nsurft)
       ! Thermal conductivity of top snowpack layer 
!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_pts(nsurft),                                                           &
    ! Number of tile points.
  surft_index(land_pts,nsurft),                                                &
    ! Index of tile points.
  nsnow(land_pts,nsurft)
    ! number of snow layers

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
   tstar_surft(land_pts,nsurft),                                               &
    ! Tile surface temperature (K)
   snow_surft(land_pts,nsurft),                                               &
    ! Snow mass on tiles (kg m-2)
   tsnow_ml(land_pts,nsurft,nsmax_ml),                                         &
    ! Snowpack layer temperatures (K).
   ds_ml(land_pts,nsurft,nsmax_ml),                                            &
    ! Depth of snowpack layers (m)
   sice_ml(land_pts,nsurft,nsmax_ml),                                          &
    ! Snowpack layer ice mass (kg m-2).
   sliq_ml(land_pts,nsurft,nsmax_ml)
    ! Snowpack layer liquid mass (kg m-2)

!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                      &
     kdtdz_ml(land_pts,nsurft)

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
   dh(land_pts,nsurft)
  
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MELTLAKE'

! useful wildcard grep 
! grep -RInE 'dtstar_surft[[:space:]]*\(.*\)[[:space:]]*='

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


DO n = 1,nsurft

   IF (l_elev_land_ice .AND. l_lice_surft(n)) THEN
      
     
      CALL meltlake_evolve(land_pts,      & !IN
           timestep,                      & !IN
           surft_pts(n),                  & !IN
           surft_index(:,n),              & !IN
           nsnow(:,n),                    & !IN
           lake_inflow(:,n),              & !IN
           ksnow0_ml(:,n),                & !IN
           tsnow_ml(:,n,1),               & !IN
           lw_down_surft(:,n),            & !IN  
           tstar_surft(:,n),              & !IN 
           sw_surft(:,n),                 & !IN already uses lake albedo
           snow_surft(:,n),               & !IN/OUT
           lake_depth_ml(:,n),            & !IN/OUT
           lake_albedo_ml(:,n),           & !IN/OUT
           lake_temp_ml(:,n),             & !IN/OUT
           exposed_water(:,n),            & !IN/OUT
           sice_ml(:,n,:),                & !IN/OUT
           sliq_ml(:,n,:),                & !IN/OUT
           ds_ml(:,n,:),                  & !IN/OUT
           kdtdz_ml(:,n))                  !OUT
           
          
                   
      !CALL relayersnow (                                                     
      !     land_pts,                     & !in                                                     
      !     surft_pts(n),                 & !in                                                     
      !     n_wtrac_jls,                  & !in                                                     
      !     surft_index(:,n),             & !in                                                 
      !     nsmax_ml,                     & !in                                                 
      !     dzsnow_ml,                    & !in                                                 
      !     rgrain0,                      & !in                                                 
      !     rho0,                         & !in                                                 
      !     sice0,                        & !in                                                 
      !     snowfall,                     & !in
      !     snowmass,                     & !in                                                   
      !     tsnow0,                       & !in                                                   
      !     wtrac_sn%sice0,               & !in                                                   
      !     nsnow(:,n),                   & !in/out                                               
      !     ds_sl_ml,                     & !in/out                                               
      !     rgrain(:,n),                  & !in/out                                               
      !     rgrainl_sl_ml,                & !in/out                                               
      !     sice_sl_ml,                   & !in/out                                               
      !     rho_snow_grnd(:,n),           & !in/out                                              
      !     sliq_sl_ml,                   & !in/out                                               
      !     tsnow_sl_ml,                  & !in/out                                               
      !     wtrac_sn%sice_sl,             & !in/out                                               
      !     wtrac_sn%sliq_sl,             & !in/out                                               
      !     rho_snow_sl_ml,               & !out                                                  
      !     snowdepth(:,n) )              & !out

                          
         
          END IF ! elev land ice tile
      
    END DO !nsurft




    
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE meltlake
END MODULE meltlake_mod
