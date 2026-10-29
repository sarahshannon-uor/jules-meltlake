! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE LID_UPDATE----------------------------------------------

! Description:
!     Calling routine for lid_evolve and relayersnow after insertion
!     of a permanent lake lid into the snowpack.
! Method:
!     Evolve virtual and permanent lake lids.
!     If a permanent lid is inserted into the snowpack, call relayersnow
!     to produce a consistent snow layer structure.
!
! Code Owner: s.r.shannon@reading.ac.uk
!
! Subroutine Interface:
MODULE lid_update_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LID_UPDATE_MOD'

CONTAINS

SUBROUTINE lid_update(land_pts,               & !IN
                      timestep,               & !IN
                      nsurft,                 & !IN
                      n_wtrac_jls,            & !IN
                      surft_pts,              & !IN
                      surft_index,            & !IN
                      frac,                   & !IN/OUT
                      nsnow,                  & !IN/OUT
                      tstar_surft,            & !IN
                      tsnow_ml,               & !IN/OUT
                      ds_ml,                  & !IN/OUT
                      sice_ml,                & !IN/OUT
                      sliq_ml,                & !IN/OUT
                      ls_snow,                & !IN/OUT
                      con_snow,               & !IN/OUT
                      ls_rain,                & !IN/OUT
                      con_rain,               & !IN/OUT
                      lake_depth_ml,          & !IN/OUT
                      lake_temp_ml,           & !IN/OUT
                      lid_temp_ml,            & !IN/OUT
                      lid_depth_ml,           & !IN/OUT
                      vlid_depth_ml,          & !IN/OUT
                      has_lake,               & !IN/OUT
                      exposed_water,          & !IN/OUT
                      has_lid,                & !IN/OUT
                      has_vlid,               & !IN/OUT
                      did_insert_lid,         & !IN/OUT
                      snow_on_lid,            & !IN/OUT
                      lid_snow_melt_flux_ml,  & !IN
                      lake_state_ml,          & !IN/OUT
                      snow_surft,             & !IN/OUT
                      melt_surft,             & !IN/OUT
                      snowdepth,              & !IN/OUT
                      rho_snow_grnd,          & !IN/OUT
                      rho_snow_ml,            & !IN/OUT
                      rgrain,                 & !IN/OUT
                      rgrainl_ml,             & !IN/OUT
                      sice_wtrac,             & !IN/OUT
                      sliq_wtrac,             & !IN/OUT
                      dhdt_lid_lake_ml,       & !IN/OUT
                      lid_snow_depth_ml,      & !IN/OUT
                      lid_snow_temp_ml,       & !IN/OUT
                      lid_snowmelt_water_ml,  & !IN/OUT
                      lid_ice_melt_flux_ml,   & !IN
                      lid_ice_meltwater_ml,   & !IN/OUT
                      lid_rain_water_ml,      & !IN/OUT
                      dhdt_lake_snow_ml,      & !IN
                      lake_inflow,            & !IN
                      ice_lens_index,         & !IN/OUT
                      ice_lens_depth,         & !IN/OUT
                      ei_surft_ml,            & !IN/OUT
                      lake_area_ml,           & !OUT
                      sice_imbal_ml,          & !IN/OUT
                      sliq_imbal_ml,          & !IN/OUT
                      en_imbal_ml,            & !IN/OUT
                      l_lice_point,           & !IN (land_pts)
                      l_lice_surft)             !IN (ntype)


USE lid_evolve_mod,          ONLY: lid_evolve
USE expand_area_mod,         ONLY: expand_area
USE merge_snowpacks_mod,     ONLY: merge_snowpacks
USE relayersnow_mod,         ONLY: relayersnow

USE jules_surface_types_mod, ONLY: ntype
USE jules_surface_mod,       ONLY: l_elev_land_ice

USE water_constants_mod,     ONLY: rho_ice, tm, rho_water

USE jules_snow_mod, ONLY:                                                    &
 rho_snow_const,                                                             &
  ! Constant density of lying snow (kg m-3) = 350
 r0
   ! Grain size for fresh snow (microns).

USE jules_meltlake_mod, ONLY: nsmax_ml, dzsnow_ml

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                       &
   land_pts,                                                                 &
     ! Total number of land points.
   nsurft,                                                                   &
    ! Number of land tiles.
   n_wtrac_jls
    ! Number of water tracers

REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
   timestep
    ! Timestep length (s).


REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
  lid_snow_melt_flux_ml(land_pts,nsurft),                                    &
    ! Snowmelt flux from snow on the lake lid (kg/m/s)
  lid_ice_melt_flux_ml(land_pts,nsurft),                                     &
    ! Melt flux from the upper surface of bare lake ice (kg/m/s)
  dhdt_lake_snow_ml(land_pts,nsurft),                                        &
    ! Lake-snow Stefan boundary movement (m ice equivalent per timestep).
  lake_inflow(land_pts,nsurft)
    ! Water supplied from snow/firn hydrology (kg m-2).
  
  
!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                       &
  surft_pts(nsurft),                                                         &
    ! Number of tile points.
  surft_index(land_pts,nsurft)
    ! Index of tile points.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN OUT) ::                                                   &
   nsnow(land_pts,nsurft)
    ! number of snow layers

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
   con_rain(land_pts),                                                       &
    ! Convective rainfall rate (kg/m2/s).
   ls_rain(land_pts),                                                        &
    ! Large-scale rainfall fall rate (kg/m2/s).
   con_snow(land_pts),                                                       &
    ! Convective frozen rainfall rate (kg/m2/s).
   ls_snow(land_pts),                                                        &
    ! Large-scale frozen precip fall rate (kg/m2/s).
   lake_depth_ml(land_pts,nsurft),                                           &
    ! Convective rainfall rate (kg/m2/s).
   lake_temp_ml(land_pts,nsurft),                                            &
    ! Temperature of melt lake (K)
   tstar_surft(land_pts,nsurft),                                             &
    ! Tile surface temperature (K)
   snow_surft(land_pts,nsurft),                                              &
    ! Snow mass on tiles (kg m-2)
   melt_surft(land_pts,nsurft),                                              &
    ! Surface snowmelt on tiles (kg/m2/s).
   snowdepth(land_pts,nsurft),                                               &
    ! Snow depth (m).
   tsnow_ml(land_pts,nsurft,nsmax_ml),                                       &
    ! Snowpack layer temperatures (K).
   ds_ml(land_pts,nsurft,nsmax_ml),                                          &
    ! Depth of snowpack layers (m)
   sice_ml(land_pts,nsurft,nsmax_ml),                                        &
    ! Snowpack layer ice mass (kg m-2).
   sliq_ml(land_pts,nsurft,nsmax_ml),                                        &
    ! Snowpack layer liquid mass (kg m-2)
   lid_depth_ml(land_pts,nsurft),                                            &
    ! Lid depth (m)
   vlid_depth_ml(land_pts,nsurft),                                           &
    ! Virtual lid depth (m)
   lid_temp_ml(land_pts,nsurft),                                             &
    ! Lid temp (K)
   rho_snow_grnd(land_pts,nsurft),                                           &
    ! Snowpack bulk density (kg/m3).
   lid_snow_depth_ml(land_pts,nsurft),                                       &
    ! Snow depth on virtual or permanent lid (m)
   lid_snow_temp_ml(land_pts,nsurft),                                        &
    ! Temp of snow on virtual or permanent lid (K)
   lid_snowmelt_water_ml(land_pts,nsurft),                                   &
    ! snowmelt from zero-layer snow scheme on top of the lid
   lid_ice_meltwater_ml(land_pts,nsurft),                                    &
    ! Liquid water produced by melting the upper surface of lake ice (m).
   lid_rain_water_ml(land_pts,nsurft),                                       &
    ! Rain water on lid (m)
   ei_surft_ml(land_pts,nsurft),                                             &
    ! Sublimation of snow (kg/m2/s).
   lake_state_ml(land_pts,nsurft),                                           &
     ! lake states (maybe should make this integer)
   frac(land_pts,nsurft),                                                    &
     ! tile frac
   lake_area_ml(land_pts,nsurft),                                            &
    ! lake area (m2) from volume area scaling
   rho_snow_ml(land_pts,nsurft,nsmax_ml),                                    &
    ! Snow layer densities for meltlake(kg/m3).
   sice_imbal_ml(land_pts),                                                  &
    ! Cumulative compensating snow ice mass adjustment (kg m-2 gridbox).
   sliq_imbal_ml(land_pts),                                                  &
    ! Cumulative compensating snow liquid water mass adjustment (kg m-2 gridbox).
   en_imbal_ml(land_pts),                                                    &
    ! Cumulative compensating snowpack enthalpy adjustment (J m-2 gridbox).
   ice_lens_index(land_pts,nsurft),                                          &
    ! Snow layer that ice lens is located in
  ice_lens_depth(land_pts,nsurft)
    ! Lens depth from surface (m)
   
   !-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                       &
   !rho_snow_ml(land_pts,nsurft,nsmax_ml),                                    &
    ! Snow layer densities for meltlake(kg/m3).
   dhdt_lid_lake_ml(land_pts,nsurft)
    ! Stefan boundary movement lid bottom and lake top (ms-1 ice equiv)
   
INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
   ! Melt-lake tile index.
  nonlake_ml = 10
   ! Non-lake elevated-ice tile index.

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
  sice0_wtrac(land_pts,n_wtrac_jls),                                          &
    ! Where nsnow=0, sice0 is the mass of the snowpack.
  sice_wtrac(land_pts,nsmax_ml,n_wtrac_jls),                                  &
    ! Water tracer ice content of snow layers (kg/m2).
  sliq_wtrac(land_pts,nsmax_ml,n_wtrac_jls)
    ! Water tracer liquid content of snow layers (kg/m2).

!ancil_info (IN)
LOGICAL, INTENT(IN) :: l_lice_point(land_pts)
LOGICAL, INTENT(IN) :: l_lice_surft(ntype)

LOGICAL, INTENT(IN OUT) ::                                                    &
 exposed_water(land_pts,nsurft),                                              &
 has_lid(land_pts,nsurft),                                                    &
 has_vlid(land_pts,nsurft),                                                   &
 has_lake(land_pts,nsurft),                                                   &
 did_insert_lid(land_pts,nsurft),                                             &
 snow_on_lid(land_pts,nsurft)

!TYPE(wtrac_sn_type) :: wtrac_sn         ! Water tracer working arrays
REAL(KIND=real_jlslsm) ::                                                    &
  frac_lake_old(land_pts)
! Lake tile fraction before area adjustment, used to conserve
   ! snow mass and energy when frac changes.

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  i,                                                                          &
    ! Land point index and loop counter.
  j,                                                                          &
    ! Tile pts loop counter.
  k,                                                                          &
    ! Tile number.
  ns

!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
!REAL(KIND=real_jlslsm) ::                                                     &
!   dh(land_pts,nsurft)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LID_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_elev_land_ice .AND. l_lice_surft(lake_ml) .AND. &
    surft_pts(lake_ml) > 0) THEN
        
!-----------------------------------------------------------------------------
! Add lateral runoff from nonlake tile to lake before lid evolution. 
! Laternal water transport is only allowed if there is exposed water
!-----------------------------------------------------------------------------
   DO k = 1, surft_pts(lake_ml)

      i = surft_index(k,lake_ml)

      IF (exposed_water(i,lake_ml) .AND.                     &
           frac(i,lake_ml) > 0.0) THEN

         !lake_depth_ml(i,lake_ml) = lake_depth_ml(i,lake_ml) &
         !     + frac(i,nonlake_ml) / frac(i,lake_ml)         &
         !     * lake_inflow(i,nonlake_ml) / rho_water

      END IF

   END DO

!-----------------------------------------------------------------------------
! grow or melt the lid
!-----------------------------------------------------------------------------
   CALL lid_evolve(land_pts,                        &
        timestep,                                   &
        surft_pts(lake_ml),                         &
        surft_index(:,lake_ml),                     &
        tstar_surft(:,lake_ml),                     &
        ls_snow,                                    &
        con_snow,                                   &
        ls_rain,                                    &
        con_rain,                                   &
        lake_temp_ml(:,lake_ml),                    &
        lake_depth_ml(:,lake_ml),                   &
        lid_temp_ml(:,lake_ml),                     &
        lid_depth_ml(:,lake_ml),                    &
        vlid_depth_ml(:,lake_ml),                   &
        has_lake(:,lake_ml),                        &
        exposed_water(:,lake_ml),                   &
        has_lid(:,lake_ml),                         &
        has_vlid(:,lake_ml),                        &
        did_insert_lid(:,lake_ml),                  &
        snow_on_lid(:,lake_ml),                     &
        lid_snow_melt_flux_ml(:,lake_ml),           &
        lake_state_ml(:,lake_ml),                   &
        dhdt_lid_lake_ml(:,lake_ml),                &
        lid_snow_depth_ml(:,lake_ml),               &
        lid_snow_temp_ml(:,lake_ml),                &
        lid_snowmelt_water_ml(:,lake_ml),           &
        lid_ice_melt_flux_ml(:,lake_ml),            &
        lid_ice_meltwater_ml(:,lake_ml),            &
        ei_surft_ml(:,lake_ml))
   
      
!-----------------------------------------------------------------------------
! Use volume area scaling to expand lake area and modify lake frac 
! and nonlake frac
!-----------------------------------------------------------------------------

   !CALL expand_area(                       &
   !     land_pts,                          &
   !     nsurft,                            &
   !     surft_pts(lake_ml),                &
   !     surft_index(:,lake_ml),            &
   !     nsnow,                             &
   !     sice_ml,                           &
   !     sliq_ml,                           &
   !     tsnow_ml,                          &
   !     exposed_water(:,lake_ml),         &
   !     frac,                              &
   !     lake_depth_ml(:,lake_ml),          &
   !     lake_area_ml(:,lake_ml),           &
   !     sice_imbal_ml,                     &
   !     sliq_imbal_ml,                     &
   !     en_imbal_ml)
   



IF (ANY(did_insert_lid(surft_index(1:surft_pts(lake_ml),lake_ml),lake_ml))) THEN

   print *, '------------------------------------------------------------'
   print *, ':-)  Inserting lid ice into snowpack '
   print *, '------------------------------------------------------------'

   snowfall(:)       = 0.0
   sice0(:)          = 0.0
   rho0(:)           = rho_ice
   tsnow0(:)         = tm
   rgrain0(:)        = 2000.0
   sice0_wtrac(:,:)  = 0.0
   

   DO j = 1, surft_pts(lake_ml)
      i = surft_index(j,lake_ml)

      IF (did_insert_lid(i,lake_ml)) THEN

         sice0(i)   = rho_ice * lid_depth_ml(i,lake_ml)
         rho0(i)    = rho_ice
         tsnow0(i)  = lid_temp_ml(i,lake_ml)
         rgrain0(i) = 2000.0

         snow_surft(i,lake_ml) = snow_surft(i,lake_ml) + sice0(i)

      END IF

   END DO

   CALL relayersnow(                                      &
        land_pts,                                         &
        surft_pts(lake_ml),                               &
        n_wtrac_jls,                                      &
        surft_index(:,lake_ml),                           &
        nsmax_ml,                                         &
        dzsnow_ml,                                        &
        rgrain0,                                          &
        rho0,                                             &
        sice0,                                            &
        snowfall,                                         &
        snow_surft(:,lake_ml),                            &
        tsnow0,                                           &
        sice0_wtrac,                                      &
        nsnow(:,lake_ml),                                 &
        ds_ml(:,lake_ml,:),                               &
        rgrain(:,lake_ml),                                &
        rgrainl_ml(:,lake_ml,:),                          &
        sice_ml(:,lake_ml,:),                             &
        rho_snow_grnd(:,lake_ml),                         &
        sliq_ml(:,lake_ml,:),                             &
        tsnow_ml(:,lake_ml,:),                            &
        sice_wtrac,                          &
        sliq_wtrac,                          &
        rho_snow_ml(:,lake_ml,:),                         &
        snowdepth(:,lake_ml) )

   
   print *, '------------------------------------------------------------'
   print *, ':-)  Inserting snow on lid into snowpack '
   print *, '------------------------------------------------------------'


   snowfall(:)       = 0.0
   sice0(:)          = 0.0
   rho0(:)           = rho_snow_const
   tsnow0(:)         = tm
   rgrain0(:)        = r0
   sice0_wtrac(:,:)  = 0.0

   DO j = 1, surft_pts(lake_ml)
      i = surft_index(j,lake_ml)

      IF (did_insert_lid(i,lake_ml) .AND.                 &
           lid_snow_depth_ml(i,lake_ml) > 0.0) THEN

         sice0(i) = rho_snow_const *                      &
              lid_snow_depth_ml(i,lake_ml)

         rho0(i)    = rho_snow_const
         tsnow0(i)  = lid_snow_temp_ml(i,lake_ml)
         rgrain0(i) = r0

         snow_surft(i,lake_ml) = snow_surft(i,lake_ml) + sice0(i)
         
      END IF

   END DO

   IF (ANY(sice0(surft_index(1:surft_pts(lake_ml),lake_ml)) > 0.0)) THEN
      
      CALL relayersnow(                                    &
           land_pts,                                       &
           surft_pts(lake_ml),                             &
           n_wtrac_jls,                                    &
           surft_index(:,lake_ml),                         &
           nsmax_ml,                                       &
           dzsnow_ml,                                      &
           rgrain0,                                        &
           rho0,                                           &
           sice0,                                          &
           snowfall,                                       &
           snow_surft(:,lake_ml),                          &
           tsnow0,                                         &
           sice0_wtrac,                                    &
           nsnow(:,lake_ml),                               &
           ds_ml(:,lake_ml,:),                             &
           rgrain(:,lake_ml),                              &
           rgrainl_ml(:,lake_ml,:),                        &
           sice_ml(:,lake_ml,:),                           &
           rho_snow_grnd(:,lake_ml),                       &
           sliq_ml(:,lake_ml,:),                           &
           tsnow_ml(:,lake_ml,:),                          &
           sice_wtrac,                        &
           sliq_wtrac,                        &
           rho_snow_ml(:,lake_ml,:),                       &
           snowdepth(:,lake_ml) )
      
   END IF

!----------------------------------------------------------------
! Lake is gone. This is the time to adjust the snowpack variables
! to conserve mass and enthalpy caused by changing frac. 
! Up until now these gridbox imbalances have been accumulinating
! Merge snowpack for lake+nonlake tile, apply imbalance
! correction. Set merged tile to nonlake tile
!----------------------------------------------------------------
   !CALL merge_snowpacks(                 &
   !     land_pts,                        &
   !     nsurft,                          &
   !     surft_pts(lake_ml),              &
   !     surft_index(:,lake_ml),          &
   !     did_insert_lid(:,lake_ml),       &
   !     frac,                            &
   !     nsnow,                           &
   !     sice_ml,                         &
   !     sliq_ml,                         &
   !     tsnow_ml,                        &
   !     ds_ml,                           &
   !     snowdepth,                       &
   !     snow_surft,                      &
   !     sice_imbal_ml,                   &
   !     sliq_imbal_ml,                   &
   !     en_imbal_ml)



!-----------------------------------------------------------------------
! Relayer merged snowpack on the nonlake tile
! No additional snow or ice is inserted.
!-----------------------------------------------------------------------

   !snowfall(:)       = 0.0
   !sice0(:)          = 0.0
   !rho0(:)           = rho_snow_const
   !tsnow0(:)         = tm
   !rgrain0(:)        = r0
   !sice0_wtrac(:,:)  = 0.0
   
   !CALL relayersnow(                                      &
   !     land_pts,                                         &
   !     surft_pts(nonlake_ml),                            &
   !     n_wtrac_jls,                                      &
   !     surft_index(:,nonlake_ml),                        &
   !     nsmax_ml,                                         &
   !     dzsnow_ml,                                        &
   !     rgrain0,                                          &
    !    rho0,                                             &
    !    sice0,                                            &
    !    snowfall,                                         &
    !    snow_surft(:,nonlake_ml),                         &
    !    tsnow0,                                           &
    !    sice0_wtrac,                                      &
     !   nsnow(:,nonlake_ml),                              &
     !   ds_ml(:,nonlake_ml,:),                            &
     !   rgrain(:,nonlake_ml),                             &
     !   rgrainl_ml(:,nonlake_ml,:),                       &
      !  sice_ml(:,nonlake_ml,:),                          &
      !  rho_snow_grnd(:,nonlake_ml),                      &
      !  sliq_ml(:,nonlake_ml,:),                          &
      !  tsnow_ml(:,nonlake_ml,:),                         &
       ! sice_wtrac,                                       &
       ! sliq_wtrac,                                       &
       ! rho_snow_ml(:,nonlake_ml,:),                      &
       ! snowdepth(:,nonlake_ml) )
         

!-----------------------------------------------------------------------
! Copy the merged and relayered nonlake snowpack to the lake tile.
! Both tiles have the same snowpack after complete freeze over.
!-----------------------------------------------------------------------

   !DO j = 1, surft_pts(lake_ml)
      
    !  i = surft_index(j,lake_ml)

     ! IF (did_insert_lid(i,lake_ml)) THEN

      !   nsnow(i,lake_ml)          = nsnow(i,nonlake_ml)

       !  sice_ml(i,lake_ml,:)      = sice_ml(i,nonlake_ml,:)
       !  sliq_ml(i,lake_ml,:)      = sliq_ml(i,nonlake_ml,:)
       !  tsnow_ml(i,lake_ml,:)     = tsnow_ml(i,nonlake_ml,:)
       !  ds_ml(i,lake_ml,:)        = ds_ml(i,nonlake_ml,:)

        ! rho_snow_ml(i,lake_ml,:)  = rho_snow_ml(i,nonlake_ml,:)
        ! rho_snow_grnd(i,lake_ml)  = rho_snow_grnd(i,nonlake_ml)

         !rgrainl_ml(i,lake_ml,:)   = rgrainl_ml(i,nonlake_ml,:)
         !rgrain(i,lake_ml)         = rgrain(i,nonlake_ml)

         !snowdepth(i,lake_ml)      = snowdepth(i,nonlake_ml)
         !snow_surft(i,lake_ml)     = snow_surft(i,nonlake_ml)

      !END IF

   !END DO
   
 !----------------------------------------------------------------
 ! Clear lake and lid states after both insertions
 !----------------------------------------------------------------

   DO j = 1, surft_pts(lake_ml)
      i = surft_index(j,lake_ml)

      IF (did_insert_lid(i,lake_ml)) THEN

         did_insert_lid(i,lake_ml)    = .FALSE.

         lake_depth_ml(i,lake_ml)     = 0.0
         lid_depth_ml(i,lake_ml)      = 0.0
         vlid_depth_ml(i,lake_ml)     = 0.0
         lid_temp_ml(i,lake_ml)       = tm

         lid_snow_depth_ml(i,lake_ml) = 0.0
         lid_snow_temp_ml(i,lake_ml)  = tm

         lake_state_ml(i,lake_ml)     = 0
         lake_area_ml(i,lake_ml)      = 0.0
         has_lake(i,lake_ml)          = .FALSE.
         exposed_water(i,lake_ml)     = .FALSE.
         has_lid(i,lake_ml)           = .FALSE.
         has_vlid(i,lake_ml)          = .FALSE.

         snow_on_lid(i,lake_ml)       = .FALSE.

    !     frac(i,lake_ml)              = 0.0
    !     frac(i,nonlake_ml)           = 1.0

         sice_imbal_ml(i)             = 0.0
         sliq_imbal_ml(i)             = 0.0
         en_imbal_ml(i)               = 0.0

         ice_lens_depth(i,lake_ml)    = -1.0
         ice_lens_index(i,lake_ml)    = -1.0

         ice_lens_depth(i,nonlake_ml) = -1.0
         ice_lens_index(i,nonlake_ml) = -1.0

      END IF

   END DO

END IF ! ANY(did_insert_lid)

END IF ! elevated lake ice tile

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE lid_update

END MODULE lid_update_mod
