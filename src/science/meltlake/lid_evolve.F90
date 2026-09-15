! ****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
! SUBROUTINE LID_EVOLVE-----------------------------------------------
! Description:
! Calculates the formation and evolution of virtual and permanent
! ice lids on melt lakes, including snow accumulation and melt on
! the snow on the lid.
!
! This is part of the JULES meltlake scheme based on Buzzard
! et al. (2018) and the MONARCHS model.  
! Buzzard et al. (2018):
!    https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017MS001155
!
! MONARCHS model:
!    https://github.com/monarchs-ice/monarchs
!
! Method:
! A thin virtual lid is seeded when exposed lake water experiences
! freezing surface conditions. The virtual lid subsequently grows or
! melts at its lower boundary according to a Stefan condition using
! conductive heat flux through the lid and turbulent heat exchange
! with the underlying lake. When the virtual lid reaches
! lid_min_depth it is converted to a permanent lid.
!
! Conductive heat transfer through the lid uses a zero-layer
! treatment. If snow is present on the lid, its thermal resistance is
! included when calculating the snow-lid interface temperature and
! conductive heat flux.
!
! Surface melt of bare virtual or permanent lake ice is calculated 
! in sf_melt and removes mass from the upper
! surface of the lid. This process is separate from melting at the
! lower lake-lid boundary.
!
! Snowfall may accumulate on virtual or permanent lids. Snow on the
! lid is represented with constant density. Surface melt removes snow
! mass and produces liquid water on top of the lid.
!
! Meltwater produced at the lid surface is kept separate from the
! underlying lake water because the ice lid is treated as impermeable.
!
! This routine maintains two separate surface meltwater stores:
!
! lid_snowmelt_water_ml
! liquid water produced by melting snow lying on the lid
!
! lid_ice_meltwater_ml
! liquid water produced by melting the upper surface of the ice lid
!
! Neither water store is added directly to lake_depth_ml. The water is
! retained on top of the lid for subsequent lateral routing.
!
! Rainwater falling on the lid is handled separately outside this
! routine and stored in lid_rain_water_ml. The rainwater and the two
! meltwater stores can later be combined in water_on_lid_depth_ml for
! lateral routing and melt lake area evolution.

! 
! Code Owner: s.r.shannon@reading.ac.uk
! Subroutine Interface:
MODULE lid_evolve_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LID_EVOLVE_MOD'

CONTAINS

  SUBROUTINE lid_evolve(land_pts,               & !IN
                      timestep,                 & !IN
                      surft_pts,                & !IN
                      surft_index,              & !IN
                      tstar_surft,              & !IN
                      ls_snow,                  & !IN/OUT
                      con_snow,                 & !IN/OUT
                      ls_rain,                  & !IN/OUT
                      con_rain,                 & !IN/OUT
                      lake_temp_ml,             & !IN/OUT
                      lake_depth_ml,            & !IN/OUT
                      lid_temp_ml,              & !IN/OUT
                      lid_depth_ml,             & !IN/OUT
                      vlid_depth_ml,            & !IN/OUT
                      has_lake,                 & !IN/OUT
                      exposed_water,            & !IN/OUT
                      has_lid,                  & !IN/OUT
                      has_vlid,                 & !IN/OUT
                      did_insert_lid,           & !IN/OUT
                      snow_on_lid,              & !IN/OUT
                      lid_snow_melt_flux_ml,  & !IN
                      lake_state_ml,            & !OUT
                      dhdt_lid_lake_ml,         & !OUT
                      lid_snow_depth_ml,        & !IN/OUT
                      lid_snow_temp_ml,         & !IN/OUT
                      lid_snowmelt_water_ml,    & !IN/OUT
                      lid_ice_melt_flux_ml,     & !IN
                      lid_ice_meltwater_ml,     & !IN/OUT
                      ei_surft_ml)
                      
                     
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


USE jules_snow_mod, ONLY:                                                    &
 snow_hcon,                                                                  & 
  ! Thermal conductivity of lying snow (Watts per m per K) = 0.265
 snow_hcap, &
  ! Thermal capacity of lying snow (J/K/m3) = 0.63e6
 rho_snow_const
  ! Constant density of lying snow (kg m-3) = 350

USE jules_meltlake_mod, ONLY: l_meltlake, nsmax_ml

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
     ! Total number of land points.
  surft_pts
    ! Number of tile points.
  
    
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep
    ! Timestep length (s).

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_index(land_pts)
    ! Index of tile points.
  
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  tstar_surft(land_pts),                                                       &
    ! Tile surface temperature (K)
  lid_snow_melt_flux_ml(land_pts),                                             &
    ! Snowmelt flux from snow on the lake lid (kg m-2 s-1).
  lid_ice_melt_flux_ml(land_pts)
    ! Melt flux from the upper surface of bare lake ice (kg m-2 s-1).

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  con_rain(land_pts),                                                          &
    ! Convective rainfall rate (kg/m2/s).
  ls_rain(land_pts),                                                           &
    ! Large-scale rainfall fall rate (kg/m2/s).
  con_snow(land_pts),                                                          &
    ! Convective frozen rainfall rate (kg/m2/s).
  ls_snow(land_pts),                                                           &
    ! Large-scale frozen precip fall rate (kg/m2/s).
  lake_temp_ml(land_pts),                                                      &
   ! Temperature of melt lake (K)
  lake_depth_ml(land_pts),                                                     &
    ! Melt lake depth (m).
  lid_temp_ml(land_pts),                                                       &
   ! Temperature of lid (K)
  lid_depth_ml(land_pts),                                                      &
    ! Depth of lid (m)
  vlid_depth_ml(land_pts),                                                     &
   ! Depth of virtual lid (m)
  lid_snow_depth_ml(land_pts),                                                 & 
    ! Depth of snow on virtual or permanent lid (m)
  lid_snow_temp_ml(land_pts),                                                  &
    ! Temp of zero layer snow on lid or vlid (K)
  lid_snowmelt_water_ml(land_pts),                                             &
    !  Liquid water from melting snow on the lake lid (m).
  lid_ice_meltwater_ml(land_pts),                                              &
    ! Liquid water from melting the upper surface of the lake ice lid (m).
  ei_surft_ml(land_pts)
    ! Sublimation of snow (kg/m2/s).

LOGICAL, INTENT(IN OUT) ::                                                     &
  has_lid(land_pts),                                                           &
  has_vlid(land_pts),                                                          &
  exposed_water(land_pts),                                                     &
  has_lake(land_pts), &
  did_insert_lid(land_pts), &
  snow_on_lid(land_pts) 

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  lake_state_ml(land_pts),                                                     &
    ! lake state flags 
   dhdt_lid_lake_ml(land_pts)
    ! Stefan boundary movement lid bottom and lake top (m T-1 ice equiv)
     
 
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  i,                                                                           &
    ! Land point index and loop counter.
  k,                                                                           &
    ! Tile number.
  n                                                                           
    ! Tile loop counter.

!-----------------------------------------------------------------------------
! Local arrays
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  kdtdz_ml, &
    ! Conductive heat flux from lid into lake
  flux_lower, &
    ! flux lake and ice
  delta_t, &
  t_lid_top, &
  dhdt, &
   ! Lake bottom-snowpack surface boundary change (m per sec ice equivalent)
  dh_ice, &
   ! Boundary change from Stefan condition (m of ice per timestep)
  dh_water, &
   ! Water equiv of boundary change from Stefan condition (m of water per timestep)
  ! mass to be removed from snowpack from Stefan boundary retreat (kgm-2)
  flux_upper_diag, &
  lid_depth_eff, &
  rain_add, &
  snow_add,&
  r_snow, &
  r_ice,&
  snowmass_lid, &
   ! Mass of snow on lid (kg/m2)
  snow_remove, &
    ! Mass of snow on lid to be removed by melting (kg/m2)
  lid_ice_melt_mass
    ! Mass removed from the upper surface of the ice lid (kg m-2).

LOGICAL ::                                                                    &
  seeded_vlid
    ! true if a virtual lid has formed
  

! Move some of these constants to jules_meltlake.nml
REAL(KIND=real_jlslsm), PARAMETER ::                                          &
  Jturb  = 1.907e-5,                                                          &
    ! Turbulent heat-transfer coefficient (m s-1 K-1/3), Eq. 16 Buzzard 
  ice_hcon = 2.2,                                                             &
    ! Thermal conductivity of ice (W m-1 K-1).
  vlid_seed_depth = 0.001,                                                    &
    ! Minimum initial thickness assigned to a virtual lid (m)
    ! Note: setting a smaller value will cause a high value outlier 
    ! for kdtdh_lake_lid_ml because the conduction path will be very small 
  lid_min_depth = 0.1,                                                        &
    !  Virtual lid thickness threshold for formation of a permanent lid (m)
    !  value from Buzzard 
  lake_min_depth = 0.1
    ! Minimum liquid lake depth used to define an exposed melt lake (m)
    ! value from Buzzard 


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LID_EVOLVE'

INTEGER              :: errcode            ! Error code
CHARACTER(LEN=80)    :: ERRMSG                ! Error message

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Lid growth by Stefan condition at the lake top and lid bottom interface
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED)                                            &
!$OMP PRIVATE(k,i,n,kdtdz_ml,flux_lower,delta_t,t_lid_top,dhdt,dh_ice,       &
!$OMP         dh_water,lid_depth_eff,rain_add,snow_add,seeded_vlid,          &
!$OMP         snowmass_lid,snow_remove,lid_ice_melt_mass)
DO k = 1,surft_pts
   i = surft_index(k)

   dhdt_lid_lake_ml(i) = 0.0
   rain_add           = 0.0
   snow_add           = 0.0
   flux_lower         = 0.0
   kdtdz_ml           = 0.0
   dhdt               = 0.0
   dh_ice             = 0.0
   dh_water           = 0.0
   did_insert_lid(i)  = .false.
   seeded_vlid        = .false.
   snow_on_lid(i)     = lid_snow_depth_ml(i) > 0.0

  
!-----------------------------------------------------------------------------
! Seed a virtual lid if conditions are freezing. Reduce the lake depth
! by the equivalent seeded depth
!-----------------------------------------------------------------------------
   IF (has_lake(i) .AND. exposed_water(i) .AND. tstar_surft(i) < tm) THEN
         
      IF (vlid_depth_ml(i) < vlid_seed_depth .AND. &
           .NOT. has_lid(i) .AND. .NOT. has_vlid(i)) THEN
            
         vlid_depth_ml(i) = vlid_seed_depth
         seeded_vlid = .true.
         
         lake_depth_ml(i) = lake_depth_ml(i) - vlid_seed_depth * rho_ice / rho_water
      END IF
      
   END IF

!-----------------------------------------------------------------------------
! update states after possible lid seeding
!-----------------------------------------------------------------------------
   CALL update_lake_states(i)
   
   
!-----------------------------------------------------------------------------
! Refreeze or melt bottom of lid using Stefan condition 
!-----------------------------------------------------------------------------
   IF (has_vlid(i) .OR. has_lid(i)) THEN

!-----------------------------------------------------------------------
! Skip Stefan growth on the timestep when a new vlid is first seeded.
! The seed is only used to initiate the lid state; Stefan growth starts
! from the next timestep. this is avoid a large jump in kdtdz_ml due to
! a very small lid thickness. This is not really necessary. 
!-----------------------------------------------------------------------      
      IF (seeded_vlid) THEN

         dhdt_lid_lake_ml(i) = 0.0
         lid_temp_ml(i) = tm
         
      ELSE

      IF (has_lid(i)) THEN

!-----------------------------------------------------------------------------
! Calculate lid conductive heat flux and temperature using the zero layer
! lid treatment, including the insulating effect of snow on the lid.
!-----------------------------------------------------------------------------
         CALL get_lid_thermo_zero_layer(                                        &
              tstar_surft(i),                                                   &
              lid_depth_ml(i),                                                  &
              has_lid(i),                                                       &
              has_vlid(i),                                                      &
              lid_snow_depth_ml(i),                                             &
              kdtdz_ml,                                                         &
              lid_temp_ml(i),                                                   &
              lid_snow_temp_ml(i))

      ELSE IF (has_vlid(i)) THEN

         CALL get_lid_thermo_zero_layer(                                        &
              tstar_surft(i),                                                   &
              vlid_depth_ml(i),                                                 &
              has_lid(i),                                                       &
              has_vlid(i),                                                      &
              lid_snow_depth_ml(i),                                             &
              kdtdz_ml,                                                         &
              lid_temp_ml(i),                                                   &
              lid_snow_temp_ml(i))
      END IF
   
      
!-----------------------------------------------------------------------------
! Lower flux, lake top --> lid bottom eqn 16
! lid bottom (tm)
!
!   ↑
!   │ flux_lower
!   │
!
! lake water (lake_temp_ml) 
!
!-----------------------------------------------------------------------------
      delta_t = lake_temp_ml(i) - tm

      flux_lower = SIGN(1.0, delta_t) * rho_water * hcapw * Jturb       &
           * ABS(delta_t)**(4.0/3.0)

      
!-----------------------------------------------------------------------------
!  Stefan boundary movement eqn 14
!  dhdt can be +ve (refreezing) or -ve (melting) 
!-----------------------------------------------------------------------------
      dhdt = (kdtdz_ml - flux_lower) / (rho_ice * lf)

      dh_ice   = timestep * dhdt
      dh_water = dh_ice * rho_ice / rho_water

     
!-----------------------------------------------------------------------------
! output Stefan boundary as diagnostic
!-----------------------------------------------------------------------------
      dhdt_lid_lake_ml(i) = dh_ice
      
!-----------------------------------------------------------------------------
!  +ve  = freezing, permanent or vlid grows 
! reduce lake depth if lid grows
! Do not refreeze more than lake depth available
!-----------------------------------------------------------------------------
      IF (dh_ice > 0.0) THEN

         dh_water = MIN(dh_water, lake_depth_ml(i))
         dh_ice   = dh_water * rho_water / rho_ice

         IF (has_lid(i)) THEN
            lid_depth_ml(i) = lid_depth_ml(i) + dh_ice
         ELSE
            vlid_depth_ml(i) = vlid_depth_ml(i) + dh_ice
         END IF

         lake_depth_ml(i) = lake_depth_ml(i) - dh_water 

             
!-----------------------------------------------------------------------------
! -ve  = melting, lid shrinks. 
! increase lake depth when virtual lid shrinks
! Note: virtual lid melts from the bottom not the top by ablation. Only 
! virtual lid can melt not the permanent lid. This is what happens in Buzzard 
!-----------------------------------------------------------------------------
      ELSE IF (dh_ice < 0.0 .AND. has_vlid(i)) THEN

         dh_ice   = ABS(dh_ice)
         dh_water = dh_ice * rho_ice / rho_water

         dh_ice   = MIN(dh_ice, vlid_depth_ml(i))
         dh_water = dh_ice * rho_ice / rho_water

         vlid_depth_ml(i) = vlid_depth_ml(i) - dh_ice
         lake_depth_ml(i) = lake_depth_ml(i) + dh_water

 !-----------------------------------------------------------------------------
 ! If the remaining virtual lid is thinner than the seed depth
 ! melt the remainder and add this tiny mass to the lake depth to conserve mass
 !-----------------------------------------------------------------------------
         IF (vlid_depth_ml(i) > 0.0 .AND. &
              vlid_depth_ml(i) < vlid_seed_depth) THEN

            lake_depth_ml(i) = lake_depth_ml(i) + &
                 vlid_depth_ml(i) * rho_ice / rho_water

            vlid_depth_ml(i) = 0.0

         END IF
         
      END IF

   END IF ! seeded_lid
   
END IF ! has_lid or has_vlid Stefan condition


!-----------------------------------------------------------------------------
! The virtual lid has grown enough in depth so convert to a permanent lid
!-----------------------------------------------------------------------------
   IF (vlid_depth_ml(i) >= lid_min_depth) THEN
      lid_depth_ml(i)  = vlid_depth_ml(i)
      vlid_depth_ml(i) = 0.0
   END IF

!-----------------------------------------------------------------------------
! There might be a permanent lid now, so update the states
!-----------------------------------------------------------------------------
   CALL update_lake_states(i)

!-----------------------------------------------------------------------------
! Melt the upper surface of a bare ice lid. This was not in Buzzard model
! but needed so tstar_surt is set to tm for a bare lid
!-----------------------------------------------------------------------------
   IF ((has_lid(i) .OR. has_vlid(i)) .AND. .NOT. snow_on_lid(i)) THEN

      lid_ice_melt_mass = MAX(lid_ice_melt_flux_ml(i), 0.0) * timestep

     
      IF (has_lid(i)) THEN

         lid_ice_melt_mass = MIN(lid_ice_melt_mass,                           &
           rho_ice * lid_depth_ml(i))

         lid_depth_ml(i) = lid_depth_ml(i) -                                  &
           lid_ice_melt_mass / rho_ice

      ELSE IF (has_vlid(i)) THEN

         lid_ice_melt_mass = MIN(lid_ice_melt_mass,                            &
              rho_ice * vlid_depth_ml(i))

         vlid_depth_ml(i) = vlid_depth_ml(i) -                                 &
              lid_ice_melt_mass / rho_ice

      END IF

      lid_ice_meltwater_ml(i) = lid_ice_meltwater_ml(i) +                      &
           lid_ice_melt_mass / rho_water

   END IF


   
   IF (has_lid(i) .OR. has_vlid(i)) THEN

!-----------------------------------------------------------------------------
! Apply surface mass loss by sublimation.
!-----------------------------------------------------------------------------
      CALL apply_lid_sublimation(i)
      
!-----------------------------------------------------------------------------
! Remove snow on lid mass by melting.  
! Use melt_surft that comes from energy balance (sf_melt) not coldsnow 
! which is for the multi layer snow scheme. 
!-----------------------------------------------------------------------------
      IF (snow_on_lid(i)) THEN
         
         ! Snow mass on the lid (kg/m2)
         snowmass_lid =  lid_snow_depth_ml(i) * rho_snow_const 
         
         ! Snow mass melted during this timestep 
         snow_remove = MAX(lid_snow_melt_flux_ml(i), 0.0) * timestep
         
         ! Only remove what exists
         snow_remove = MIN(snow_remove, snowmass_lid)

         ! Remove melted snow mass
         snowmass_lid = snowmass_lid - snow_remove

         ! Update snow on lid depth
         lid_snow_depth_ml(i) = snowmass_lid / rho_snow_const

         ! What to do with the melted water on the lid ???
         ! Could this be a lake, on top of a lid, on top of a lake
         ! Don't add the water to the lake depth because the lid is impermeable. 
         ! Keep the water depth on lid in a seperate variable. 
        
         IF (snow_remove > 0.0) THEN
            lid_snowmelt_water_ml(i) = lid_snowmelt_water_ml(i)           &
                 + snow_remove / rho_water
         END IF
            
   END IF
      
END IF ! has_lid or has_vlid 
      
!-----------------------------------------------------------------------------
! Set a flag if lake has fully frozen over. Flag is used later to insert the 
! lid ice and snow on lid, into the snowpack. 
!-----------------------------------------------------------------------------
   IF (lake_depth_ml(i) <= 1.0e-5 .AND. lid_depth_ml(i) > 0.0) THEN
      lake_depth_ml(i)  = 0.0
      did_insert_lid(i) = .true.
   END IF


!-----------------------------------------------------------------------------
! update states now the lid and lake depths have changed
!-----------------------------------------------------------------------------   
   CALL update_lake_states(i)

  
!-----------------------------------------------------------------------------
! Reset inactive temperatures and depths 
!-----------------------------------------------------------------------------

IF (.NOT. has_lake(i)) THEN
   lake_temp_ml(i) = tm
END IF

IF (.NOT. has_lid(i)) THEN
   lid_temp_ml(i)  = tm
   lid_depth_ml(i) = 0.0
END IF

IF (.NOT. has_vlid(i)) THEN
   vlid_depth_ml(i) = 0.0
END IF

IF ((.NOT. has_lid(i)) .AND. (.NOT. has_vlid(i))) THEN
   lid_snow_depth_ml(i) = 0.0
   snow_on_lid(i)       = .false.
END IF

!-----------------------------------------------------------------------------
! Lake state flag, just for plotting output. To do: move this to an seperate 
! subroutine
!0 = no lake
!1 = exposed water
!2 = bare virtual lid
!3 = bare permanent lid
!4 = snow on virtual lid
!5 = snow on permanent lid
!-----------------------------------------------------------------------------
lake_state_ml(i) = 0

IF (exposed_water(i)) THEN
   lake_state_ml(i) = 1
ELSE IF (has_vlid(i) .AND. .NOT. snow_on_lid(i)) THEN
   lake_state_ml(i) = 2
ELSE IF (has_lid(i) .AND. .NOT. snow_on_lid(i)) THEN
   lake_state_ml(i) = 3
ELSE IF (has_vlid(i) .AND. snow_on_lid(i)) THEN
   lake_state_ml(i) = 4
ELSE IF (has_lid(i) .AND. snow_on_lid(i)) THEN
   lake_state_ml(i) = 5
END IF

   
END DO ! land_pts
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


CONTAINS

SUBROUTINE update_lake_states(i)

   INTEGER, INTENT(IN) :: i

   has_lid(i)  = (lid_depth_ml(i)  >= lid_min_depth)

   has_vlid(i) = (vlid_depth_ml(i) >= vlid_seed_depth) .AND.                 &
              (vlid_depth_ml(i) <  lid_min_depth) .AND.                      &
              .NOT. has_lid(i)

   
   has_lake(i) = (lake_depth_ml(i) > lake_min_depth) .OR.                    &
                 ((lake_depth_ml(i) > 0.0) .AND.                             &
                  (has_lid(i) .OR. has_vlid(i)))

   exposed_water(i) = (lake_depth_ml(i) > lake_min_depth) .AND.              &
                      .NOT. has_lid(i) .AND. .NOT. has_vlid(i)


   
END SUBROUTINE update_lake_states

SUBROUTINE get_lid_thermo_zero_layer(                                        &
    tstar_surft,                                                             &
    lid_thickness,                                                           &
    has_lid,                                                                 &
    has_vlid,                                                                &
    lid_snow_depth_ml,                                                       &
    kdtdz_ml,                                                                &
    lid_temp_ml,                                                             &
    lid_snow_temp_ml)
!--------------------------------------------------------------------
! Get lid conductive flux (kdtdz_ml) and lid temperature with/without snow on top 
! Bare lid: lid-top temperature = tile surface temperature, capped at tm.
! Snow-covered lid: lid-top temperature = snow-lid interface temperature,
! found by splitting the temperature drop between snow surface and lid
! bottom (tm) according to snow and lid resistances in series.
!
! Without snow on lid                ! With snow on lid
!                                    !
!                                    ! snow surface (tstar_surft)
!                                    !
!                                    !   ↑
!                                    !   │ snow and lid resistances
!                                    !   │ used to calculate t_lid_top
!                                    !
! lid top (tstar_surft)              ! lid top / snow-lid interface (t_lid_top)
!                                    !
!   ↑                                !   ↑
!   │ kdtdz_ml through lid           !   │ kdtdz_ml through lid
!   │                                !   │
!                                    !
! lid bottom (tm)                    ! lid bottom (tm)
!--------------------------------------------------------------------
  
  REAL(KIND=real_jlslsm), INTENT(IN) ::                                      &
    tstar_surft,                                                             &
      ! Tile surface temperature (K).
    lid_thickness,                                                           &
      ! Thickness of permanent lid or virtual lid (m).
    lid_snow_depth_ml
      ! Snow depth on top of lid or virtual lid (m).

  LOGICAL, INTENT(IN) ::                                                     &
    has_lid,                                                                 &
      ! True if permanent lid exists.
    has_vlid
      ! True if virtual lid exists.

  REAL(KIND=real_jlslsm), INTENT(OUT) ::                                     &
    kdtdz_ml,                                                                &
      ! Conductive heat flux through lid or virtual lid (W m-2).
    lid_temp_ml,                                                             &
      ! Mean lid temperature (K). For virtual lid this is returned as tm.
    lid_snow_temp_ml
      ! Mean temperature of snow on lid (K)
  
  REAL(KIND=real_jlslsm) ::                                                  &
    t_lid_top,                                                               &
      ! Temperature at top of lid, or snow-lid interface if snow is present.
    r_snow,                                                                  &
      ! Thermal resistance of snow layer (m2 K W-1).
    r_lid
      ! Thermal resistance of lid or virtual lid (m2 K W-1).

  !--------------------------------------------------------------------
  ! Snow-covered lid or virtual lid, zero-layer treatment
  ! Snow provides insulation only, no prognostic snow temperature.
  !--------------------------------------------------------------------
  IF (lid_snow_depth_ml > 0.0) THEN

     r_snow = lid_snow_depth_ml / snow_hcon
     r_lid  = lid_thickness     / ice_hcon

     ! Snow-lid interface temperature from resistances in series
     t_lid_top = tm - (tm - MIN(tstar_surft, tm)) * r_lid / (r_snow + r_lid)

     lid_snow_temp_ml = 0.5 * (MIN(tstar_surft, tm) +                       &
                               t_lid_top)

  !--------------------------------------------------------------------
  ! Bare lid or bare virtual lid
  !--------------------------------------------------------------------
  ELSE

     t_lid_top = MIN(tstar_surft, tm)
     lid_snow_temp_ml = tm

  END IF

  kdtdz_ml = ice_hcon * (tm - t_lid_top) / lid_thickness

  IF (has_lid) THEN
     lid_temp_ml = 0.5 * (t_lid_top + tm)
  ELSE IF (has_vlid) THEN
     lid_temp_ml = tm
  END IF

  
END SUBROUTINE get_lid_thermo_zero_layer

SUBROUTINE apply_lid_sublimation(i)

!-----------------------------------------------------------------------------
! Remove surface mass associated with ei_surft_ml.
! Mass is removed from the material exposed to the atmosphere in the
! following order: snow on lid, permanent lid, virtual lid.
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  i
    ! Land point index.

REAL(KIND=real_jlslsm) ::                                                     &
  sublim_remove,                                                              &
    ! Mass to be removed by sublimation during timestep (kg m-2).
  snowmass_lid
    ! Mass of snow on lid (kg m-2).

sublim_remove = MAX(ei_surft_ml(i), 0.0) * timestep

IF (sublim_remove > 0.0) THEN

  IF (snow_on_lid(i)) THEN

    snowmass_lid = lid_snow_depth_ml(i) * rho_snow_const

    sublim_remove = MIN(sublim_remove, snowmass_lid)

    snowmass_lid = snowmass_lid - sublim_remove

    lid_snow_depth_ml(i) = snowmass_lid / rho_snow_const

    snow_on_lid(i) = (lid_snow_depth_ml(i) > 0.0)

  ELSE IF (has_lid(i)) THEN

    sublim_remove = MIN(sublim_remove,                                        &
         rho_ice * lid_depth_ml(i))

    lid_depth_ml(i) = lid_depth_ml(i) -                                       &
         sublim_remove / rho_ice

  ELSE IF (has_vlid(i)) THEN

    sublim_remove = MIN(sublim_remove,                                        &
         rho_ice * vlid_depth_ml(i))

    vlid_depth_ml(i) = vlid_depth_ml(i) -                                     &
         sublim_remove / rho_ice

  END IF

END IF

END SUBROUTINE apply_lid_sublimation

END SUBROUTINE lid_evolve
END MODULE lid_evolve_mod
