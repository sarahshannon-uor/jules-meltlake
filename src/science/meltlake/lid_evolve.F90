! ****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE LID_EVOLVE-----------------------------------------------
! Description:
!    virtual lid: can grow and melt
!    permanent lid: can grow, but not shrink by this Stefan melt term
! Method:

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
                      nsnow,                    & !IN/OUT
                      ds_ml,                    & !IN/OUT
                      sice_ml,                  & !IN/OUT
                      sliq_ml,                  & !IN/OUT
                      tsnow_ml,                 & !IN/OUT
                      snow_surft,               & !OUT
                      lake_state_ml,            & !OUT
                      dhdt_lid_lake_ml,         & !OUT
                      lid_snow_depth_ml,        & !IN/OUT
                      lid_snow_temp_ml,         & !IN/OUT
                      snowfall,                 & !OUT
                      tsnow0,                   & !OUT
                      rho0,                     & !OUT
                      rgrain0,                  & !OUT
                      sice0)                      !OUT
                     
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
 rho_snow_fresh,                                                             &
    ! Density of fresh snow (kg per m**3) = 100
 snow_hcon,                                                                  & 
  ! Thermal conductivity of lying snow (Watts per m per K) = 0.265
 snow_hcap
 ! Thermal capacity of lying snow (J/K/m3) = 0.63e6
 
USE jules_meltlake_mod, ONLY: l_meltlake, nsmax_ml

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

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_index(land_pts)
    ! Index of tile points.
  
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  tstar_surft(land_pts)
    ! Tile surface temperature (K)

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN OUT) ::                                                     &
   nsnow(land_pts)                                                        
    ! Number of snow layers.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  con_rain(land_pts),                                                         &
    ! Convective rainfall rate (kg/m2/s).
  ls_rain(land_pts),                                                          &
    ! Large-scale rainfall fall rate (kg/m2/s).
  con_snow(land_pts),                                                         &
    ! Convective frozen rainfall rate (kg/m2/s).
  ls_snow(land_pts),                                                          &
    ! Large-scale frozen precip fall rate (kg/m2/s).
  lake_temp_ml(land_pts),                                                      &
   ! Temperature of melt lake (K)
  lake_depth_ml(land_pts),                                                     &
    ! Convective rainfall rate (kg/m2/s).
  lid_temp_ml(land_pts),                                                       &
   ! Temperature of lid (K)
  lid_depth_ml(land_pts),                                                      &
    ! Depth of lid (m)
  vlid_depth_ml(land_pts),                                                     &
   ! Depth of virtual lid (m)
  sice_ml(land_pts,nsmax_ml),                                                  &
    ! Ice content of snow layers (kg/m2)
  sliq_ml(land_pts,nsmax_ml),                                                  &
    ! Liquid content of snow layers (kg/m2)
  ds_ml(land_pts, nsmax_ml),                                                   &
    ! snowpack top level depth (m)
  tsnow_ml(land_pts, nsmax_ml),                                                & 
    ! snowpack level temperatures (K)
  lid_snow_depth_ml(land_pts),                                                 & 
    ! Depth of snow on virtual or permanent lid (m)
  lid_snow_temp_ml(land_pts),                                                  &
    ! Temp of zero layer snow on lid or vlid (K)
  snow_surft(land_pts)
! snow mass can change because lid is injected into snowpack

LOGICAL, INTENT(IN OUT) ::                                                    &
  has_lid(land_pts),                                                          &
  has_vlid(land_pts),                                                         &
  exposed_water(land_pts),                                                    &
  has_lake(land_pts), &
  did_insert_lid(land_pts), &
  snow_on_lid(land_pts) 

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  lake_state_ml(land_pts),                                                    &
    ! logicals for lake state
  !snow_surft(land_pts),                                                       & 
    ! snow mass can change because lid is injected into snowpack 
  dhdt_lid_lake_ml(land_pts),                                                 & 
    ! Stefan boundary movement lid bottom and lake top (ms-1 ice equiv)
  snowfall(land_pts), &
  tsnow0(land_pts), &
  rho0(land_pts),&
  rgrain0(land_pts),& 
  sice0(land_pts)
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  i,                                                                           &
    ! Land point index and loop counter.
!  j,                                                                           &
    ! Tile pts loop counter.
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
  !dm, &
   ! mass to be removed from snowpack from Stefan boundary retreat (kgm-2)
  !flux_upper_diag, &
  lid_depth_eff, &
  rain_add, &
  snow_add
  

REAL, PARAMETER :: Jturb  = 1.907e-5
                 ! Turbulent heat flux factor (ms⁻¹ K⁻¹/3) Eqn 16 Buzzard 

LOGICAL :: seeded_vlid

REAL(KIND=real_jlslsm), PARAMETER :: ice_hcon = 2.2

REAL(KIND=real_jlslsm), PARAMETER :: vlid_seed_depth = 0.001!1.0e-3 ! use a seed for now 

REAL(KIND=real_jlslsm), PARAMETER :: lid_min_depth = 0.1
REAL(KIND=real_jlslsm), PARAMETER :: lake_min_depth = 0.1

REAL(KIND=real_jlslsm) :: r_snow, r_lid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LID_EVOLVE'

INTEGER              :: errcode            ! Error code
CHARACTER(LEN=80) :: ERRMSG                ! Error message

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Lid growth by Stefan condition at the lake top / lid bottom interface
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED)                                            &
!$OMP PRIVATE(k,i,n,kdtdz_ml,flux_lower,delta_t,t_lid_top,dhdt,dh_ice,       &
!$OMP         dh_water,lid_depth_eff,rain_add,snow_add,seeded_vlid)
DO k = 1,surft_pts
   i = surft_index(k)

   dhdt_lid_lake_ml(i) = 0.0
   rain_add   = 0.0
   snow_add   = 0.0
   flux_lower = 0.0
   kdtdz_ml   = 0.0
   dhdt       = 0.0
   dh_ice     = 0.0
   dh_water   = 0.0
   did_insert_lid(i) = .false.
   snow_on_lid(i)    = .false.
   seeded_vlid       = .false.

!-----------------------------------------------------------------------------
! If lake exists and is freezing, seed a thin lid. Taking the seed depth 
! from no where and not adjusting the lake depth. It's small though
!-----------------------------------------------------------------------------
   IF (has_lake(i).AND.tstar_surft(i) < tm) THEN

      IF (vlid_depth_ml(i) < vlid_seed_depth .AND. .NOT. has_lid(i)) THEN
         vlid_depth_ml(i) = vlid_seed_depth
         seeded_vlid = .true.
      END IF

   END IF

!-----------------------------------------------------------------------------
! update states again after possible lid seeding
!-----------------------------------------------------------------------------
   CALL update_lake_states(i)

!-----------------------------------------------------------------------------
! Do Stefan condition for virtual and permanent lid 
!-----------------------------------------------------------------------------
   IF (has_vlid(i) .OR. has_lid(i)) THEN


!-----------------------------------------------------------------------
! Skip Stefan growth on the timestep when a new vlid is first seeded.
! The seed is only used to initiate the lid state; Stefan growth starts
! from the next timestep. this is avoid a large jump kdtdz_ml due to
! a very small lid thickness 
!-----------------------------------------------------------------------      
      IF (seeded_vlid) THEN

         dhdt_lid_lake_ml(i) = 0.0
         lid_temp_ml(i) = tm

      ELSE

    
!-----------------------------------------------------------------------------
! Get temp at top of lid. If no snow then top of lid is tstar_surft
!-----------------------------------------------------------------------------
      t_lid_top = MIN(tstar_surft(i), tm)

      IF (has_lid(i)) THEN  
         lid_temp_ml(i) = 0.5 * (t_lid_top + tm)
      ELSE IF (has_vlid(i)) THEN
         lid_temp_ml(i) = tm
      END IF

     ! IF (has_lid(i)) THEN

      !   CALL get_lid_thermo_zero_layer(                                        &
       !       tstar_surft(i),                                                   &
       !       lid_depth_ml(i),                                                  &
       !       has_lid(i),                                                       &
       !       has_vlid(i),                                                      &
       !       lid_snow_depth_ml(i),                                             &
       !       kdtdz_ml,                                                         &
       !       lid_temp_ml(i))

      !ELSE IF (has_vlid(i)) THEN

      !   CALL get_lid_thermo_zero_layer(                                        &
      !        tstar_surft(i),                                                   &
      !        vlid_depth_ml(i),                                                 &
      !        has_lid(i),                                                       &
      !        has_vlid(i),                                                      &
      !        lid_snow_depth_ml(i),                                             &
      !        kdtdz_ml,                                                         &
      !        lid_temp_ml(i))
         
      !END IF
   
      
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
! Heat conducted between the lid bottom and lid top (atmosphere)
! direction format: from the T1 to the T2
! setting up to work for virtual and permantent lid. Eqn 14.
!      
! lid top (t_lid_top i.e. tstar_surft)
!
!   ↑
!   │ kdtdz
!   │
!
! lid bottom (tm)
!-----------------------------------------------------------------------------
      IF (has_lid(i)) THEN
         lid_depth_eff = lid_depth_ml(i)
      ELSE IF (has_vlid(i)) THEN
         lid_depth_eff = MAX(vlid_depth_ml(i), 0.001)
      END IF

      kdtdz_ml = ice_hcon * (tm - t_lid_top) / lid_depth_eff

!-----------------------------------------------------------------------------
!  Stefan boundary movement Eqn 14
!  dhdt can be +ve (refreezing) or -ve (melting) 
! 
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
! -ve  = melting, lid shrinks - !!! This only applies to vlid
! increase lake depth is lid shrinks
! Do not melt more vlid than exists 
!-----------------------------------------------------------------------------
      ELSE IF (dh_ice < 0.0 .AND. has_vlid(i)) THEN

         dh_ice   = ABS(dh_ice)
         dh_water = dh_ice * rho_ice / rho_water

         dh_ice   = MIN(dh_ice, vlid_depth_ml(i))
         dh_water = dh_ice * rho_ice / rho_water

         vlid_depth_ml(i) = vlid_depth_ml(i) - dh_ice
         lake_depth_ml(i) = lake_depth_ml(i) + dh_water

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
! Snow and rain on lid
!-----------------------------------------------------------------------------
   IF (has_lid(i) .OR. has_vlid(i)) THEN
      IF (ls_rain(i) > 0.0 .OR. con_rain(i) > 0.0 .OR.                      &
           ls_snow(i) > 0.0 .OR. con_snow(i) > 0.0) THEN

         rain_add = ls_rain(i) + con_rain(i)
         snow_add = ls_snow(i) + con_snow(i)


         WRITE(*,*) '--- SNOW ON LID UPDATE ---'
         WRITE(*,'(A,I8)')    'timestep                = ', timestep_number
         WRITE(*,'(A,I8)')    'i                       = ', i
         WRITE(*,'(A,L2)')    'has_lid                 = ', has_lid(i)
         WRITE(*,'(A,L2)')    'has_vlid                = ', has_vlid(i)
         WRITE(*,'(A,L2)')    'snow_on_lid             = ', snow_on_lid(i)
         WRITE(*,'(A,F12.6)') 'ls_snow                 = ', ls_snow(i)
         WRITE(*,'(A,F12.6)') 'con_snow                = ', con_snow(i)
         WRITE(*,'(A,F12.6)') 'snow_add                = ', snow_add
         WRITE(*,'(A,F12.6)') 'lid_snow_depth before   = ', lid_snow_depth_ml(i)

         IF (has_lid(i)) THEN
          
            lid_snow_depth_ml(i) = lid_snow_depth_ml(i) +                      &
                 snow_add * timestep /rho_snow_fresh 
            
         ELSE IF (has_vlid(i)) THEN
          
            lid_snow_depth_ml(i) = lid_snow_depth_ml(i) +                      &
                 snow_add * timestep /rho_snow_fresh 
          
         END IF

         !snow_on_lid(i) = (lid_snow_depth_ml(i) > 0.0)
         
         WRITE(*,'(A,F12.6)') 'lid_snow_depth after    = ', lid_snow_depth_ml(i)
         WRITE(*,'(A,L2)')    'snow_on_lid after      = ', snow_on_lid(i)

         
         
         lake_depth_ml(i) = lake_depth_ml(i)                                &
              + rain_add * timestep / rho_water

         ls_rain(i)  = 0.0
         con_rain(i) = 0.0
         ls_snow(i)  = 0.0
         con_snow(i) = 0.0

         

      END IF


!-----------------------------------------------------------------------------
! Remove sublimation and melt from snow on lid
!-----------------------------------------------------------------------------
      !IF (snow_on_lid(i)) THEN

       !  snowmass_lid = rho_snow_fresh * lid_snow_depth_ml(i)

        ! snow_remove = (ei_surft(i) + melt_surft(i)) * timestep
        ! snow_remove = MIN(snow_remove, snowmass_lid)

        ! snow_melt = MIN(melt_surft(i) * timestep, snowmass_lid)

        ! snowmass_lid = snowmass_lid - snow_remove
        ! lid_snow_depth_ml(i) = snowmass_lid / rho_snow_fresh

   ! melted snow becomes lake water
         !lake_depth_ml(i) = lake_depth_ml(i) + snow_melt / rho_water

         !snow_on_lid(i) = (lid_snow_depth_ml(i) > 0.0)

      !END IF


      
   END IF ! has_lid or has_vlid snow and rain on lid 
      
!-----------------------------------------------------------------------------
! Get inputs to pass to relayersnow to insert fully refrozen permanent lid into
! the top of the snowpack
!-----------------------------------------------------------------------------
   CALL prepare_lid_insertion_for_relayer(i)
      
!-----------------------------------------------------------------------------
! update states now the lid and lake depths have changed
!-----------------------------------------------------------------------------   
   CALL update_lake_states(i)

   IF (did_insert_lid(i)) THEN
      WRITE(*,'(A,I8)')    'timestep = ', timestep_number
      WRITE(*,'(A,L2)')    'did_insert_lid = ',did_insert_lid(i)
      WRITE(*,'(A,F12.6)') 'lake_depth_ml = ', lake_depth_ml(i)
      WRITE(*,'(A,F12.6)') 'lid_depth_ml  = ', lid_depth_ml(i)
      WRITE(*,'(A,F12.6)') 'vlid_depth_ml = ', vlid_depth_ml(i)
      WRITE(*,'(A,L2)')    'has_lake      = ', has_lake(i)
      WRITE(*,'(A,L2)')    'has_lid       = ', has_lid(i)
      WRITE(*,'(A,L2)')    'has_vlid      = ', has_vlid(i)
      WRITE(*,'(A,L2)')    'exposed_water = ', exposed_water(i)
      
   END IF
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
END IF

!-----------------------------------------------------------------------------
! Lake state flag, just for plotting output, move this to extra subroutine
! after CALL meltlake_evolve, CALL lid_evolve 
! 0 = no lake
! 1 = exposed water
! 2 = virtual lid
! 3 = permanent lid
!-----------------------------------------------------------------------------
   lake_state_ml(i) = 0

   IF (exposed_water(i)) THEN
      lake_state_ml(i) = 1
   END IF

   IF (has_vlid(i)) THEN
      lake_state_ml(i) = 2
   END IF

   IF (has_lid(i)) THEN
      lake_state_ml(i) = 3
   END IF
   
   
END DO ! land_pts
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


CONTAINS

SUBROUTINE update_lake_states(i)

   INTEGER, INTENT(IN) :: i

   has_lid(i)  = (lid_depth_ml(i)  > lid_min_depth)

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
    lid_temp_ml)

!--------------------------------------------------------------------
! Get lid conductive flux (kdtdz_ml) and lid temperature with/without snow on top 
! Bare lid: lid-top temperature = tile surface temperature, capped at tm.
! Snow-covered lid: lid-top temperature = snow-lid interface temperature,
! found by splitting the temperature drop between snow surface and lid
! bottom (tm) according to snow and lid resistances in series.
! 
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
    lid_temp_ml
      ! Mean lid temperature (K). For virtual lid this is returned as tm.

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
     r_lid  = lid_thickness      / ice_hcon

     ! Snow-lid interface temperature from resistances in series
     t_lid_top = tm - (tm - MIN(tstar_surft, tm)) * r_lid / (r_snow + r_lid)

  !--------------------------------------------------------------------
  ! Bare lid or bare virtual lid
  !--------------------------------------------------------------------
  ELSE

     t_lid_top = MIN(tstar_surft, tm)

  END IF

  kdtdz_ml = ice_hcon * (tm - t_lid_top) / lid_thickness

  IF (has_lid) THEN
     lid_temp_ml = 0.5 * (t_lid_top + tm)
  ELSE IF (has_vlid) THEN
     lid_temp_ml = tm
  END IF

END SUBROUTINE get_lid_thermo_zero_layer





!SUBROUTINE get_lid_top_temp(i, lid_thickness, t_lid_top)
  !--------------------------------------------------------------------
  ! Get lid top temperature with or without snow on top of lid
  ! If no snow on lid then lid top temp is tile surface temp
  ! If snow on lid then get lid top temp by partitionaing the total temp drop
  ! between tile surface temp and lid bottom (tm) using the
  ! resistances in series.  
  !--------------------------------------------------------------------
  
 !  INTEGER, INTENT(IN) :: i
 !  REAL(KIND=real_jlslsm), INTENT(IN)  :: lid_thickness
 !  REAL(KIND=real_jlslsm), INTENT(OUT) :: t_lid_top

   ! thermal resistance of lid and snow on top of lid
 !  REAL(KIND=real_jlslsm) :: r_snow, r_lid 

 !  WRITE(*,*) '--- GET_LID_TOP_TEMP ---'
 !  WRITE(*,'(A,I8)')    'timestep                = ', timestep_number
 !  WRITE(*,'(A,I8)')    'i                       = ', i
 !  WRITE(*,'(A,F12.6)') 'lid_thickness           = ', lid_thickness
 !  WRITE(*,'(A,F12.6)') 'lid_snow_depth_ml       = ', lid_snow_depth_ml(i)
 !  WRITE(*,'(A,F12.6)') 'tstar_surft (C)         = ', tstar_surft(i) - 273.15
 !  WRITE(*,'(A,F12.6)') 'tm (C)                  = ', tm - 273.15

   ! if snow on lid then lid temp then 
 !  IF (lid_snow_depth_ml(i) > 0.0) THEN

      ! thicker snow larger resistance 
  !    r_snow = lid_snow_depth_ml(i) / snow_hcon
      
   !   r_lid  = lid_thickness / ice_hcon

   !   t_lid_top = tm - (tm - MIN(tstar_surft(i), tm)) *                      &
   !        r_lid / (r_snow + r_lid)
      
   !   WRITE(*,'(A)')       'branch                  = snow on lid'
   !   WRITE(*,'(A,F12.6)') 'snow_hcon               = ', snow_hcon
   !   WRITE(*,'(A,F12.6)') 'ice_hcon                = ', ice_hcon
   !   WRITE(*,'(A,F12.6)') 'r_snow                  = ', r_snow
   !   WRITE(*,'(A,F12.6)') 'r_lid                   = ', r_lid
   !   WRITE(*,'(A,F12.6)') 't_lid_top (C)           = ', t_lid_top - 273.15
  
   !ELSE

   !   t_lid_top = MIN(tstar_surft(i), tm)
   !   WRITE(*,'(A)')       'branch                  = bare lid'
   !   WRITE(*,'(A,F12.6)') 't_lid_top (C)           = ', t_lid_top - 273.15

   !END IF

!END SUBROUTINE get_lid_top_temp

SUBROUTINE prepare_lid_insertion_for_relayer(i)
  !--------------------------------------------------------------------
  ! Get inputs to pass to relayersnow to insert the permanent lid 
  ! into the snowpack. Relayersnow is called in meltlake_mod
  !--------------------------------------------------------------------

  INTEGER, INTENT(IN) :: i

  
  !did_insert_lid(i) = .false.
  
  IF (lake_depth_ml(i) <= 1.0e-5 .AND. lid_depth_ml(i) > 0.0) THEN

     
     !--------------------------------------------------------------------
     ! Modified arguments to pass to relayersnow 
     !--------------------------------------------------------------------
     snowfall(i)   = 0.0
     snow_surft(i) = snow_surft(i) + rho_ice * lid_depth_ml(i)
     sice0(i)      = lid_depth_ml(i) * rho_ice
     tsnow0(i)     = lid_temp_ml(i)
     rho0(i)       = rho_ice
     rgrain0(i)    = 2000.0                                               
     
         
     !--------------------------------------------------------------------
     ! Lid has now been transferred into the snowpack
     !--------------------------------------------------------------------
     lid_depth_ml(i)   = 0.0
     lid_temp_ml(i)    = tm
     lake_depth_ml(i)  = 0.0
     did_insert_lid(i) = .true.
     snow_on_lid(i)    = .false.
     
  END IF

END SUBROUTINE prepare_lid_insertion_for_relayer

  

END SUBROUTINE lid_evolve
END MODULE lid_evolve_mod
