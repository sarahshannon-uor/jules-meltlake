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
                      nsnow,                    & !IN/OUT
                      ds_ml,                    & !IN/OUT
                      sice_ml,                  & !IN/OUT
                      sliq_ml,                  & !IN/OUT
                      tsnow_ml,                 & !IN/OUT
                      snow_surft,               & !OUT
                      lake_state_ml,            & !OUT
                      dhdt_lid_lake_ml)           !OUT
                     
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
  tsnow_ml(land_pts, nsmax_ml)
    ! snowpack level temperatures (K)


LOGICAL, INTENT(IN OUT) ::                                                    &
  has_lid(land_pts),                                                          &
  has_vlid(land_pts),                                                         &
  exposed_water(land_pts),                                                    &
  has_lake(land_pts), &
  did_insert_lid(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  lake_state_ml(land_pts),                                                    &
    ! logicals for lake state
  snow_surft(land_pts),                                                       & 
    ! snow mass can change because lid is injected into snowpack 
  dhdt_lid_lake_ml(land_pts)
    ! Stefan boundary movement lid bottom and lake top (ms-1 ice equiv)

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
  !expon_term,                                                                  &
    ! exponential term in albedo eqn (m)
  !ex,                                                                          &
    ! exponential term in albedo eqn (m)
  kdtdz_ml, &
    ! Conductive heat flux from lid into lake
  flux_lower, &
    ! flux lake and ice
  !flux_upper, &
    ! flux between atmosphere and lake
  delta_t, &
  !dt, &
  !dTdt, &
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


REAL(KIND=real_jlslsm), PARAMETER :: hcon_lid = 2.2

REAL(KIND=real_jlslsm), PARAMETER :: lid_seed_depth = 0.01!1.0e-3 ! use a seed for now (should add vlid)

REAL(KIND=real_jlslsm), PARAMETER :: lid_min_depth = 0.1
REAL(KIND=real_jlslsm), PARAMETER :: lake_min_depth = 0.1


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LID_EVOLVE'

INTEGER              :: errcode            ! Error code
CHARACTER(LEN=80) :: ERRMSG             ! Error message

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Lid growth by Stefan condition at the lake top / lid bottom interface
!-----------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------
! If lid or virtual lid is present, add precipitation to the lid/lake system
! snow -> lid or virtual lid
! rain -> lake
! Refactor: Needs to happen before CALL snow instead of setting win=0, snowfall=0
! in snow 
!-----------------------------------------------------------------------------
   IF (has_lid(i) .OR. has_vlid(i)) THEN

      IF (ls_rain(i) > 0.0 .OR. con_rain(i) > 0.0 .OR.                      &
           ls_snow(i) > 0.0 .OR. con_snow(i) > 0.0) THEN

         rain_add = ls_rain(i) + con_rain(i)
         snow_add = ls_snow(i) + con_snow(i)

         IF (has_lid(i)) THEN
            lid_depth_ml(i) = lid_depth_ml(i)                               &
                 + snow_add * timestep / rho_ice
            
            
         ELSE IF (has_vlid(i)) THEN
            vlid_depth_ml(i) = vlid_depth_ml(i)                             &
                 + snow_add * timestep / rho_ice
           
         END IF

         
         lake_depth_ml(i) = lake_depth_ml(i)                                &
              + rain_add * timestep / rho_water

         ls_rain(i)  = 0.0
         con_rain(i) = 0.0
         ls_snow(i)  = 0.0
         con_snow(i) = 0.0

      END IF
   END IF
!-----------------------------------------------------------------------------
! If lake exists and is freezing, seed a thin lid. Taking the seed depth 
! from no where and not adjusting the lake depth. It's small though
!-----------------------------------------------------------------------------
   IF (has_lake(i).AND.tstar_surft(i) < tm) THEN

      IF (vlid_depth_ml(i) < lid_seed_depth .AND. .NOT. has_lid(i)) THEN
         vlid_depth_ml(i) = lid_seed_depth
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
      
!-----------------------------------------------------------------------------
! Lid bulk temp is mean of tile surface temp and tm
!-----------------------------------------------------------------------------
      t_lid_top = MIN(tstar_surft(i), tm)
      lid_temp_ml(i) = 0.5 * (t_lid_top + tm)

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
! Heat conducted between the lid bottom and lid top
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
         lid_depth_eff = MAX(vlid_depth_ml(i), 0.01)!vlid_depth_ml(i) prevent seed shock
      END IF

      kdtdz_ml = hcon_lid * (tm - t_lid_top) / lid_depth_eff

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

 END IF ! has_lid or has_vlid Stefan condition
     
!-----------------------------------------------------------------------------
! The virtual lid has grown enough in depth so convert to a permanent lid
!-----------------------------------------------------------------------------
   IF (vlid_depth_ml(i) >= lid_min_depth) THEN
      lid_depth_ml(i)  = vlid_depth_ml(i)
      vlid_depth_ml(i) = 0.0
   END IF

!-----------------------------------------------------------------------------
! Insert fully refrozen permanent lid into the top of the snowpack
!-----------------------------------------------------------------------------
   CALL insert_permanent_lid_into_top_snow(i)
      
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


CONTAINS

SUBROUTINE update_lake_states(i)

   INTEGER, INTENT(IN) :: i

   has_lid(i)  = (lid_depth_ml(i)  > lid_min_depth)

   has_vlid(i) = (vlid_depth_ml(i) >= lid_seed_depth) .AND.                    &
              (vlid_depth_ml(i) <  lid_min_depth) .AND.                        &
              .NOT. has_lid(i)
   
   !has_vlid(i) = (vlid_depth_ml(i) >= lid_seed_depth) .AND. .NOT. has_lid(i)

   has_lake(i) = (lake_depth_ml(i) > lake_min_depth) .OR.                    &
                 ((lake_depth_ml(i) > 0.0) .AND.                             &
                  (has_lid(i) .OR. has_vlid(i)))

   exposed_water(i) = (lake_depth_ml(i) > lake_min_depth) .AND.              &
                      .NOT. has_lid(i) .AND. .NOT. has_vlid(i)

END SUBROUTINE update_lake_states
  
SUBROUTINE insert_permanent_lid_into_top_snow(i)
  ! Add the fully refrozen permanent lid into the top of the snowpack.
  ! This is called before relayersnow (again) so it can reorganise the
  ! column 

  INTEGER, INTENT(IN) :: i

  REAL(KIND=real_jlslsm) :: sice_add
  INTEGER :: n

  did_insert_lid(i) = .false.
  
  IF (lake_depth_ml(i) <= 1.0e-5 .AND. lid_depth_ml(i) > 0.0) THEN

     WRITE(*,*) '--- INSERT PERMANENT LID INTO TOP SNOW ---'
     WRITE(*,*) 'i                              = ', i
     WRITE(*,'(A,F16.8)') 'lake_depth_ml        = ', lake_depth_ml(i)
     WRITE(*,'(A,F16.8)') 'lid_depth_ml         = ', lid_depth_ml(i)
     WRITE(*,'(A,F16.8)') 'lid_temp_ml          = ', lid_temp_ml(i) - 273.15
     WRITE(*,*) 'nsnow before                   = ', nsnow(i)

     !--------------------------------------------------------------------
     ! Add lid depth and mass into the existing top snow layer
     !--------------------------------------------------------------------
     sice_add = rho_ice * lid_depth_ml(i)

     ds_ml(i,1)    = ds_ml(i,1)    + lid_depth_ml(i)
     sice_ml(i,1)  = sice_ml(i,1)  + sice_add
     tsnow_ml(i,1) = lid_temp_ml(i)

     WRITE(*,'(A,F16.8)') 'sice_add            = ', sice_add
     WRITE(*,'(A,F16.8)') 'ds_ml(i,1)          = ', ds_ml(i,1)
     WRITE(*,'(A,F16.8)') 'sice_ml(i,1)        = ', sice_ml(i,1)
     WRITE(*,'(A,F16.8)') 'sliq_ml(i,1)        = ', sliq_ml(i,1)
     WRITE(*,'(A,F16.8)') 'tsnow_ml(i,1)       = ', tsnow_ml(i,1) - 273.15

     !--------------------------------------------------------------------
     ! Update bulk snow diagnostics
     !--------------------------------------------------------------------
     
     snow_surft(i) = 0.0

     DO n = 1, nsnow(i)
        snow_surft(i) = snow_surft(i) + sice_ml(i,n) + sliq_ml(i,n)
     END DO

     !--------------------------------------------------------------------
     ! Lid has now been transferred into the snowpack
     !--------------------------------------------------------------------
     lid_depth_ml(i) = 0.0
     lid_temp_ml(i)  = tm

     did_insert_lid(i) = .true.
     
     
     WRITE(*,'(A,F16.8)') 'snow_surft after    = ', snow_surft(i)
     WRITE(*,'(A,F16.8)') 'lid_depth_ml after  = ', lid_depth_ml(i)
     WRITE(*,'(A,F16.8)') 'lid_temp_ml after   = ', lid_temp_ml(i) - 273.15
     WRITE(*,*) '------------------------------------------'

  END IF

END SUBROUTINE insert_permanent_lid_into_top_snow  

  

END SUBROUTINE lid_evolve
END MODULE lid_evolve_mod
