! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
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
                      n_wtrac_jls,              & !IN                           
                      surft_pts,                & !IN
                      surft_index,              & !IN
                      nsnow,                    & !IN
                      lake_inflow,              & !IN
                      kdtdz_ml,                 & !OUT
                      lw_down_surft,            & !IN  
                      tstar_surft,              & !IN/OUT 
                      sw_surft,                 & !IN  
                      tsnow_ml,                 & !IN/OUT
                      ds_ml,                    & !IN/OUT
                      sice_ml,                  & !IN/OUT
                      sliq_ml,                  & !IN/OUT
                      ls_snow,                  & !IN/OUT
                      con_snow,                 & !IN/OUT
                      ls_rain,                  & !IN/OUT
                      con_rain,                 & !IN/OUT
                      lake_depth_ml,            & !IN/OUT
                      lake_albedo_ml,           & !IN/OUT
                      lake_temp_ml,             & !IN/OUT
                      lid_temp_ml,              & !IN/OUT
                      lid_depth_ml,             & !IN/OUT
                      vlid_depth_ml,            & !IN/OUT
                      has_lake,                 & !IN/OUT
                      exposed_water,            & !IN/OUT
                      has_lid,                  & !IN/OUT
                      has_vlid,                 & !IN/OUT
                      did_insert_lid,           & !IN/OUT
                      ksnow0_ml,                & !IN
                      lake_state_ml,            & !OUT
                      snow_surft,               & !IN/OUT
                      snowdepth,                & !IN/OUT
                      rho_snow_grnd,            & !IN/OUT
                      rho_snow_ml,              & !OUT
                      dhdt_lake_snow_ml,        & !OUT
                      dhdt_lid_lake_ml,         & !OUT
                      !Ancil info (IN)
                      l_lice_point,             & !IN (land_pts)
                      l_lice_surft)               !IN (ntype)  

USE meltlake_evolve_mod,     ONLY: meltlake_evolve
USE lid_evolve_mod,          ONLY: lid_evolve
USE relayersnow_mod,         ONLY: relayersnow

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

USE jules_meltlake_mod, ONLY: nsmax_ml, dzsnow_ml

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
   nsurft,                                                                     &
    ! Number of land tiles.
   n_wtrac_jls
    ! Number of water tracers

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep              ! Timestep length (s).



REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  sw_surft(land_pts,nsurft),                                                   &                            ! Net shortwave radiation on tile (W/m2) using lake albedo
  lw_down_surft(land_pts,nsurft),                                              &
    ! Surface downward LW radiation on tiles (W/m2), jules_land_sf_implicit.F90
  lake_inflow(land_pts,nsurft),                                                &
    ! melt water from snowpack (kgm-2)
  ksnow0_ml(land_pts,nsurft)
    ! Thermal conductivity of top snowpack layer
  
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
INTEGER, INTENT(IN OUT) ::                                                     &
   nsnow(land_pts,nsurft)
    ! number of snow layers

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
   con_rain(land_pts),                                                         &
    ! Convective rainfall rate (kg/m2/s).
   ls_rain(land_pts),                                                          &
    ! Large-scale rainfall fall rate (kg/m2/s).
   con_snow(land_pts),                                                         &
    ! Convective frozen rainfall rate (kg/m2/s).
   ls_snow(land_pts),                                                          &
    ! Large-scale frozen precip fall rate (kg/m2/s).
   lake_depth_ml(land_pts,nsurft),                                             &
    ! Convective rainfall rate (kg/m2/s).
   lake_albedo_ml(land_pts,nsurft),                                            &
    ! Albedo of lake
   lake_temp_ml(land_pts,nsurft),                                              &
    ! Temperature of melt lake (K)
   tstar_surft(land_pts,nsurft),                                               &
    ! Tile surface temperature (K)
   snow_surft(land_pts,nsurft),                                                &
    ! Snow mass on tiles (kg m-2)
   snowdepth(land_pts,nsurft),                                                 &
     ! Snow depth (m).
   tsnow_ml(land_pts,nsurft,nsmax_ml),                                         &
    ! Snowpack layer temperatures (K).
   ds_ml(land_pts,nsurft,nsmax_ml),                                            &
    ! Depth of snowpack layers (m)
   sice_ml(land_pts,nsurft,nsmax_ml),                                          &
    ! Snowpack layer ice mass (kg m-2).
   sliq_ml(land_pts,nsurft,nsmax_ml),                                          &
    ! Snowpack layer liquid mass (kg m-2)
   lid_depth_ml(land_pts,nsurft),                                              &
    ! Lid depth (m)
   vlid_depth_ml(land_pts,nsurft),                                             &
    ! Virtual lid depth (m)
   lid_temp_ml(land_pts,nsurft),                                               &
    ! Lid temp (K)
   rho_snow_grnd(land_pts,nsurft)
    ! Snowpack bulk density (kg/m3).
!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
   kdtdz_ml(land_pts,nsurft),                                                 &
     ! conductive heat flux from lake into snowpack 
   lake_state_ml(land_pts,nsurft),                                            & 
     ! lake states
   rho_snow_ml(land_pts,nsurft,nsmax_ml),                                     &
    ! Snow layer densities for meltlake(kg/m3).
   dhdt_lake_snow_ml(land_pts,nsurft),                                        & 
    ! Stefan boundary movement lake bottom and snowpack top (ms-1 ice equiv) 
   dhdt_lid_lake_ml(land_pts,nsurft)
    ! Stefan boundary movement lid bottom and lake top (ms-1 ice equiv)
!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                     &
  snowfall(land_pts),                                                         &
   ! Total frozen precip reaching the ground in timestep
   ! (kg/m2) - includes any canopy unloading.
  sice0(land_pts),                                                            &
   ! Ice content of fresh snow (kg/m2).
   ! Where nsnow=0, sice0 is the mass of the snowpack.
  tsnow0(land_pts),                                                           & 
    ! Temperature of fresh snow (K).
  rgrain(land_pts,nsurft),                                                    &
    ! Snow surface grain size (microns).
  rgrainl_ml(land_pts,nsurft,nsmax_ml),                                       & 
    ! Snow layer grain size for meltlake (microns).
  rgrain0(land_pts),                                                          &
    ! Fresh snow grain size (microns).
  rho0(land_pts),                                                             & 
    ! Density of fresh snow (kg/m3).
    ! Where nsnow=0, rho0 is the density of the snowpack.
  sice0_wtrac(land_pts,n_wtrac_jls),                                          &     ! Water tracer ice content of fresh snow (kg/m2).
    ! Where nsnow=0, sice0 is the mass of the snowpack.
  sice_wtrac(land_pts,nsmax_ml,n_wtrac_jls),                                     &
    ! Water tracer ice content of snow layers (kg/m2).
  sliq_wtrac(land_pts,nsmax_ml,n_wtrac_jls)
    ! Water tracer liquid content of snow layers (kg/m2).
       
!ancil_info (IN)
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)
LOGICAL, INTENT(IN) :: l_lice_surft(ntype)

LOGICAL, INTENT(IN OUT) ::                                                    &
 exposed_water(land_pts,nsurft),                                              &
 has_lid(land_pts,nsurft), &
 has_vlid(land_pts,nsurft), &
 has_lake(land_pts,nsurft), &
 did_insert_lid(land_pts, nsurft)

!TYPE(wtrac_sn_type) :: wtrac_sn         ! Water tracer working arrays

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
!   dh(land_pts,nsurft)



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
           sw_surft(:,n),                 & !IN 
           lid_depth_ml(:,n),             & !IN
           vlid_depth_ml(:,n),            & !IN
           snow_surft(:,n),               & !IN/OUT
           ls_snow,                       & !IN/OUT
           con_snow,                      & !IN/OUT
           ls_rain,                       & !IN/OUT
           con_rain,                      & !IN/OUT
           lake_depth_ml(:,n),            & !IN/OUT
           lake_albedo_ml(:,n),           & !IN/OUT
           lake_temp_ml(:,n),             & !IN/OUT
           has_lake(:,n),                 & !IN/OUT
           exposed_water(:,n),            & !IN/OUT
           has_lid(:,n),                  & !IN
           has_vlid(:,n),                 & !IN
           sice_ml(:,n,:),                & !IN/OUT
           sliq_ml(:,n,:),                & !IN/OUT
           ds_ml(:,n,:),                  & !IN/OUT
           kdtdz_ml(:,n),                 & !OUT
           dhdt_lake_snow_ml(:,n))          !OUT
           
          
      CALL lid_evolve(land_pts,           & !IN
           timestep,                      & !IN
           surft_pts(n),                  & !IN
           surft_index(:,n),              & !IN
           tstar_surft(:,n),              & !IN 
           ls_snow,                       & !IN/OUT
           con_snow,                      & !IN/OUT
           ls_rain,                       & !IN/OUT
           con_rain,                      & !IN/OUT
           lake_temp_ml(:,n),             & !IN/OUT
           lake_depth_ml(:,n),            & !IN/OUT
           lid_temp_ml(:,n),              & !IN/OUT
           lid_depth_ml(:,n),             & !IN/OUT
           vlid_depth_ml(:,n),            & !IN/OUT
           has_lake(:,n),                 & !IN/OUT
           exposed_water(:,n),            & !IN/OUT
           has_lid(:,n),                  & !IN/OUT
           has_vlid(:,n),                 & !IN/OUT
           did_insert_lid(:,n),           & !IN/OUT
           nsnow(:,n),                    & !IN/OUT
           ds_ml(:,n,:),                  & !IN/OUT
           sice_ml(:,n,:),                & !IN/OUT
           sliq_ml(:,n,:),                & !IN/OUT
           tsnow_ml(:,n,:),               & !IN/OUT
           snow_surft(:,n),               & !OUT
           lake_state_ml(:,n),            & !OUT
           dhdt_lid_lake_ml(:,n))           !OUT


      
      ! This is runtime wasteful - it runs though all land_pts ever though I know
      ! which ones have did_insert_lid=true. Come back to this.  
      IF (ANY(did_insert_lid(surft_index(1:surft_pts(n),n),n))) THEN


         WRITE(*,*) '--- BEFORE RELAYERSNOW INSERT LID ---'
         WRITE(*,*) 'timestep_number   = ', timestep_number
         WRITE(*,*) 'has_lake               = ', has_lake(1,9)
         WRITE(*,*) 'has_lid                = ', has_lid(1,9)
         WRITE(*,*) 'has_vlid               = ', has_vlid(1,9)
         WRITE(*,*) 'exposed_water          = ', exposed_water(1,9)
         WRITE(*,*) 'did_insert_lid         = ', did_insert_lid(1,9)
         
         
         WRITE(*,*) 'nsnow(1,9)        = ', nsnow(1,9)

         
         
         WRITE(*,'(A,F16.8)') 'snowdepth(1,9) = ', snowdepth(1,9)
         WRITE(*,'(A,F16.8)') 'snow_surft(1,9)= ', snow_surft(1,9)
         WRITE(*,'(A,F16.8)') 'ds top         = ', ds_ml(1,9,1)
         WRITE(*,'(A,F16.8)') 'sice top       = ', sice_ml(1,9,1)
         WRITE(*,'(A,F16.8)') 'sliq top       = ', sliq_ml(1,9,1)
         WRITE(*,'(A,F16.8)') 'tsnow top oC   = ', tsnow_ml(1,9,1) - 273.15
         WRITE(*,'(A,F16.8)') 'ds(1,9,2)        = ', ds_ml(1,9,2)
         WRITE(*,'(A,F16.8)') 'sice(1,9,2)      = ', sice_ml(1,9,2)
         WRITE(*,'(A,F16.8)') 'sliq(1,9,2)      = ', sliq_ml(1,9,2)
         WRITE(*,'(A,F16.8)') 'tsnow(1,9,2) oC  = ', tsnow_ml(1,9,2) - 273.15
         WRITE(*,'(A,F16.8)') 'ds(1,9,3)        = ', ds_ml(1,9,3)
         WRITE(*,'(A,F16.8)') 'sice(1,9,3)      = ', sice_ml(1,9,3)
         WRITE(*,'(A,F16.8)') 'sliq(1,9,3)      = ', sliq_ml(1,9,3)
         WRITE(*,'(A,F16.8)') 'tsnow(1,9,3) oC  = ', tsnow_ml(1,9,3) - 273.15


         
         snowfall(:)         = 0.0
         sice0(:)            = 0.0
         tsnow0(:)           = tm
         !snowmass(:)         = snow_surft(:,n)
         sice0_wtrac(:,:)    = 0.0
         sice_wtrac(:,n,:)   = 0.0
         sliq_wtrac(:,n,:)   = 0.0

      CALL relayersnow ( land_pts,       &
           surft_pts(n),                 &
           n_wtrac_jls,                  &
           surft_index(:,n),             &
           nsmax_ml,                     &
           dzsnow_ml,                    &
           rgrain0,                      &
           rho0,                         &
           sice0,                        &
           snowfall,                     &
           snow_surft(:,n),              &
           tsnow0,                       &
           sice0_wtrac,                  &
           nsnow(:,n),                   &
           ds_ml(:,n,:),                 &
           rgrain(:,n),                  &
           rgrainl_ml(:,n,:),            &
           sice_ml(:,n,:),               &
           rho_snow_grnd(:,n),           &
           sliq_ml(:,n,:),               &
           tsnow_ml(:,n,:),              &
           sice_wtrac(:,n,:),            &
           sliq_wtrac(:,n,:),            &
           rho_snow_ml(:,n,:),           &
           snowdepth(:,n) )

      
      WRITE(*,*) '--- AFTER RELAYERSNOW INSERT LID---'
      
      WRITE(*,*) 'timestep_number            = ', timestep_number
      WRITE(*,*) 'nsnow(1,9)                 = ', nsnow(1,9)
      WRITE(*,'(A,F16.8)') 'snowdepth(1,9) = ', snowdepth(1,9)
      WRITE(*,'(A,F16.8)') 'snow_surft(1,9)= ', snow_surft(1,9)
     
      WRITE(*,'(A,F16.8)') 'ds(1,9,2)        = ', ds_ml(1,9,2)
      WRITE(*,'(A,F16.8)') 'sice(1,9,2)      = ', sice_ml(1,9,2)
      WRITE(*,'(A,F16.8)') 'sliq(1,9,2)      = ', sliq_ml(1,9,2)
      WRITE(*,'(A,F16.8)') 'tsnow(1,9,2) oC  = ', tsnow_ml(1,9,2) - 273.15
      WRITE(*,'(A,F16.8)') 'ds(1,9,3)        = ', ds_ml(1,9,3)
      WRITE(*,'(A,F16.8)') 'sice(1,9,3)      = ', sice_ml(1,9,3)
      WRITE(*,'(A,F16.8)') 'sliq(1,9,3)      = ', sliq_ml(1,9,3)
      WRITE(*,'(A,F16.8)') 'tsnow(1,9,3) oC  = ', tsnow_ml(1,9,3) - 273.15
      WRITE(*,*) 'did_insert_lid             = ', did_insert_lid(1,9)
      
      
   END IF ! did_insert_lid
         
   END IF ! elev land ice tile
      
END DO !nsurft

!WRITE(*,'(A,F18.8)') 'ds_ml        = ', ds_ml(1,9,1)
      !WRITE(*,'(A,F18.8)') 'tsnow_ml     = ', tsnow_ml(1,9,1)-273.15
      !WRITE(*,'(A,F18.8)') 'kdtdz_ml     = ', kdtdz_ml(1,9)
      !WRITE(*,'(A,F18.8)') 'tstar        = ', tstar_surft(1,9)-273.15
      !WRITE(*,'(A,F18.8)') 'ksnow0_ml    = ', ksnow0_ml(1,9)
      !WRITE(*,'(A,F18.8)') 'lake_temp_ml = ', lake_temp_ml(1,9)-273.15
      !WRITE(*,'(A,F18.8)') 'lid_temp_ml  = ', lid_temp_ml(1,9)-273.15
      !WRITE(*,'(A,F18.8)') 'lake_depth_ml= ', lake_depth_ml(1,9)
      !WRITE(*,'(A,F18.8)') 'lid_depth_ml = ', lid_depth_ml(1,9)


    
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE meltlake
END MODULE meltlake_mod
