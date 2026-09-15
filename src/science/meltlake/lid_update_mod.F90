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
                      lid_ice_melt_flux_ml,   & ! IN
                      lid_ice_meltwater_ml,   & ! IN/OUT
                      ei_surft_ml,            & !IN/OUT
                      l_lice_point,           & !IN (land_pts)
                      l_lice_surft)             !IN (ntype)  
                     

USE lid_evolve_mod,          ONLY: lid_evolve
USE relayersnow_mod,         ONLY: relayersnow

USE jules_surface_types_mod, ONLY: ntype
USE jules_surface_mod,       ONLY: l_elev_land_ice

USE water_constants_mod,     ONLY: rho_ice, tm

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
INTEGER, INTENT(IN) ::                                                         &
   land_pts,                                                                   &
     ! Total number of land points.
   nsurft,                                                                     &
    ! Number of land tiles.
   n_wtrac_jls
    ! Number of water tracers

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
   timestep
    ! Timestep length (s).


REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  lid_snow_melt_flux_ml(land_pts, nsurft),                                    &
    ! Snowmelt flux from snow on the lake lid (kg/m/s)
  lid_ice_melt_flux_ml(land_pts,nsurft)
    ! Melt flux from the upper surface of bare lake ice (kg/m/s)
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
   lake_temp_ml(land_pts,nsurft),                                              &
    ! Temperature of melt lake (K)
   tstar_surft(land_pts,nsurft),                                               &
    ! Tile surface temperature (K)
   snow_surft(land_pts,nsurft),                                                &
    ! Snow mass on tiles (kg m-2)
   melt_surft(land_pts,nsurft),                                               & 
    ! Surface snowmelt on tiles (kg/m2/s).
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
   rho_snow_grnd(land_pts,nsurft),                                             &
    ! Snowpack bulk density (kg/m3).
   lid_snow_depth_ml(land_pts, nsurft),                                        &
    ! Snow depth on virtual or permanent lid (m) 
   lid_snow_temp_ml(land_pts, nsurft),                                         &
    ! Temp of snow on virtual or permanent lid (K) 
   lid_snowmelt_water_ml(land_pts, nsurft),                                  &
    ! snowmelt from zero-layer snow scheme on top of the lid 
   lid_ice_meltwater_ml(land_pts,nsurft),                                      &
    ! Liquid water produced by melting the upper surface of lake ice (m).
   ei_surft_ml(land_pts, nsurft),                                              &
    ! Sublimation of snow (kg/m2/s).
   lake_state_ml(land_pts,nsurft)
     ! lake states (maybe should make this integer)
!-----------------------------------------------------------------------------
! Array arguments with intent(out)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
   rho_snow_ml(land_pts,nsurft,nsmax_ml),                                     &
    ! Snow layer densities for meltlake(kg/m3).
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
 has_lid(land_pts,nsurft), &
 has_vlid(land_pts,nsurft), &
 has_lake(land_pts,nsurft), &
 did_insert_lid(land_pts,nsurft), &
 snow_on_lid(land_pts,nsurft)

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
  n, &                                                                           
    ! Tile loop counter.
  ns
!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
!REAL(KIND=real_jlslsm) ::                                                      &
!   dh(land_pts,nsurft)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LID_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO n = 1,nsurft

   IF (l_elev_land_ice .AND. l_lice_surft(n) .AND. surft_pts(n) > 0) THEN

      sice0_wtrac(:,:)  = 0.0
      sice_wtrac(:,n,:) = 0.0
      sliq_wtrac(:,n,:) = 0.0
         
      CALL lid_evolve(land_pts,           &
           timestep,                      &
           surft_pts(n),                  &
           surft_index(:,n),              &
           tstar_surft(:,n),              &
           ls_snow,                       &
           con_snow,                      &
           ls_rain,                       &
           con_rain,                      &
           lake_temp_ml(:,n),             &
           lake_depth_ml(:,n),            &
           lid_temp_ml(:,n),              &
           lid_depth_ml(:,n),             &
           vlid_depth_ml(:,n),            &
           has_lake(:,n),                 &
           exposed_water(:,n),            &
           has_lid(:,n),                  &
           has_vlid(:,n),                 &
           did_insert_lid(:,n),           &
           snow_on_lid(:,n),              &
           lid_snow_melt_flux_ml(:,n),    &
           lake_state_ml(:,n),            &
           dhdt_lid_lake_ml(:,n),         &
           lid_snow_depth_ml(:,n),        &
           lid_snow_temp_ml(:,n),         &
           lid_snowmelt_water_ml(:,n),    &
           lid_ice_melt_flux_ml(:,n),     & 
           lid_ice_meltwater_ml(:,n),     & 
           ei_surft_ml(:,n))
          
      
      IF (ANY(did_insert_lid(surft_index(1:surft_pts(n),n),n))) THEN

         print *, '------------------------------------------------------------'
         print *, ':-)  Inserting lid ice into snowpack '
         print *, '------------------------------------------------------------'
         
         snowfall(:)       = 0.0
         sice0(:)          = 0.0
         rho0(:)           = rho_ice
         tsnow0(:)         = tm
         rgrain0(:)        = 2000.0
         sice0_wtrac(:,:)  = 0.0

         DO j = 1, surft_pts(n)
            i = surft_index(j,n)

            IF (did_insert_lid(i,n)) THEN

               sice0(i)   = rho_ice * lid_depth_ml(i,n)
               rho0(i)    = rho_ice
               tsnow0(i)  = lid_temp_ml(i,n)
               rgrain0(i) = 2000.0

               snow_surft(i,n) = snow_surft(i,n) + sice0(i)

            END IF

         END DO
         
         CALL relayersnow(                                      &
              land_pts,                                         &
              surft_pts(n),                                     &
              n_wtrac_jls,                                      &
              surft_index(:,n),                                 &
              nsmax_ml,                                         &
              dzsnow_ml,                                        &
              rgrain0,                                          &
              rho0,                                             &
              sice0,                                            &
              snowfall,                                         &
              snow_surft(:,n),                                  &
              tsnow0,                                           &
              sice0_wtrac,                                      &
              nsnow(:,n),                                       &
              ds_ml(:,n,:),                                     &
              rgrain(:,n),                                      &
              rgrainl_ml(:,n,:),                                &
              sice_ml(:,n,:),                                   &
              rho_snow_grnd(:,n),                               &
              sliq_ml(:,n,:),                                   &
              tsnow_ml(:,n,:),                                  &
              sice_wtrac(:,n,:),                                &
              sliq_wtrac(:,n,:),                                &
              rho_snow_ml(:,n,:),                               &
              snowdepth(:,n) )

         print *, '------------------------------------------------------------'
         print *, ':-)  Inserting snow on lid into snowpack '
         print *, '------------------------------------------------------------'

         
         snowfall(:)       = 0.0
         sice0(:)          = 0.0
         rho0(:)           = rho_snow_const
         tsnow0(:)         = tm
         rgrain0(:)        = r0
         sice0_wtrac(:,:)  = 0.0

         DO j = 1, surft_pts(n)
            i = surft_index(j,n)

            IF (did_insert_lid(i,n) .AND.                         &
                 lid_snow_depth_ml(i,n) > 0.0) THEN

               sice0(i) = rho_snow_const *                        &
                    lid_snow_depth_ml(i,n)

               rho0(i)    = rho_snow_const
               tsnow0(i)  = lid_snow_temp_ml(i,n)
               rgrain0(i) = r0

               snow_surft(i,n) = snow_surft(i,n) + sice0(i)

            END IF

         END DO

         IF (ANY(sice0(surft_index(1:surft_pts(n),n)) > 0.0)) THEN
            
            CALL relayersnow(                                    &
                 land_pts,                                       &
                 surft_pts(n),                                   &
                 n_wtrac_jls,                                    &
                 surft_index(:,n),                               &
                 nsmax_ml,                                       &
                 dzsnow_ml,                                      &
                 rgrain0,                                        &
                 rho0,                                           &
                 sice0,                                          &
                 snowfall,                                       &
                 snow_surft(:,n),                                &
                 tsnow0,                                         &
                 sice0_wtrac,                                    &
                 nsnow(:,n),                                     &
                 ds_ml(:,n,:),                                   &
                 rgrain(:,n),                                    &
                 rgrainl_ml(:,n,:),                              &
                 sice_ml(:,n,:),                                 &
                 rho_snow_grnd(:,n),                             &
                 sliq_ml(:,n,:),                                 &
                 tsnow_ml(:,n,:),                                &
                 sice_wtrac(:,n,:),                              &
                 sliq_wtrac(:,n,:),                              &
                 rho_snow_ml(:,n,:),                             &
                 snowdepth(:,n) )
            
         END IF
         !----------------------------------------------------------------
         ! Clear lake and lid states after both insertions 
         !----------------------------------------------------------------
   
         DO j = 1, surft_pts(n)
            i = surft_index(j,n)

            IF (did_insert_lid(i,n)) THEN

               did_insert_lid(i,n)    = .false.
               
               lake_depth_ml(i,n)     = 0.0
               lid_depth_ml(i,n)      = 0.0
               vlid_depth_ml(i,n)     = 0.0
               lid_temp_ml(i,n)       = tm

               lid_snow_depth_ml(i,n) = 0.0
               lid_snow_temp_ml(i,n)  = tm
               
               lake_state_ml(i,n)     = 0
               has_lake(i,n)          = .false.
               exposed_water(i,n)     = .false.
               has_lid(i,n)           = .false.
               has_vlid(i,n)          = .false.

               snow_on_lid(i,n)       = .false.
               

            END IF

         END DO

      END IF
      
   END IF
   
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE lid_update

END MODULE lid_update_mod
