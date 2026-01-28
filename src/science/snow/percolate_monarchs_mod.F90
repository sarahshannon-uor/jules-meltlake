! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!-----------------------------------------------------------------------------
! Description:
!   MONARCHS-style percolation and perching for a single snow column.
!
!   This routine is called from SNOWPACK after before percolation & refreezing 
!   
!
!   Refactor
!
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

MODULE percolate_monarchs_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PERCOLATE_MONARCHS_MOD'

CONTAINS

  SUBROUTINE percolate_monarchs( i, nsnow, nsmax, n_wtrac_jls, timestep,     &
                               csnow, ds, tsnow,                             &
                               sice, sliq,                                   &
                               sice_wtrac, sliq_wtrac,                       &
                               win, win_wtrac,                               &
                               sf_diag, surft_n, ice_lens_depth )


!-----------------------------------------------------------------------------
! Inputs and supporting modules
!-----------------------------------------------------------------------------

USE jules_water_tracers_mod, ONLY: l_wtrac_jls, wtrac_calc_ratio_fn_jules
  ! l_wtrac_jls enables tracer bookkeeping.
  ! wtrac_calc_ratio_fn_jules provides tracer-to-water ratio in a store.

USE jules_snow_mod,          ONLY: rho_firn_pore_closure
  ! Density threshold at which pores are assumed closed (kg m-3).
  ! Used here as an "ice lens / impermeable layer" criterion.

USE water_constants_mod, ONLY:                                                 &
  lf,                                                                          &
   ! latent heat of fusion (J kg-1)
  rho_ice,                                                                     &
   ! density of ice (kg m-3)
  rho_water,                                                                   &
   ! density of water (kg m-3)
  tm
   ! melting temperature (K)

!USE meltlake_vars_mod, ONLY: meltlake_vars_type  

USE sf_diags_mod, ONLY: strnewsfdiag
  ! Surface diagnostics, updated for freezing if enabled.

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  i,                                                                           &  
    ! Land point index for diagnostics arrays inside sf_diag.
  nsnow,                                                                       &
    ! Number of active snow layers in this column (<= nsmax).
  nsmax,                                                                       &
    ! Maximum possible snow layers, sets array bounds.
  surft_n,                                                                     &
    ! Surface tile index for diagnostics.
  n_wtrac_jls                                                                 
    ! Number of water tracers.
  
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
    ! Model timestep (s). Used to convert masses to flux diagnostics.
  csnow(:,:)
    ! Areal heat capacity of each snow layer (J K-1 m-2).

!-----------------------------------------------------------------------------
! Scalar arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  win
    ! Liquid water mass in transit downward (kg m-2).
    ! On entry: water available to enter layer 1 this timestep (rain+melt etc).
    ! On exit: remaining outflow from the column (kg m-2), diagnosed by caller.

!TYPE(meltlake_vars_type), INTENT(IN OUT) :: meltlake_vars

TYPE (strnewsfdiag), INTENT(IN OUT) :: sf_diag
  ! Diagnostics structure, updated here for refreezing flux.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  ds(:,:),                                                                   &
    ! Layer thicknesses (m).
  tsnow(:,:),                                                                &
    ! Layer temperatures (K).
  sice(:,:),                                                                 &
    ! Ice mass per unit area in each layer (kg m-2).
  sliq(:,:),                                                                 &
    ! Liquid mass per unit area in each layer (kg m-2).
  sice_wtrac(:,:,:),                                                         &
    ! Tracer mass associated with ice store (kg m-2 tracer units).
  sliq_wtrac(:,:,:),                                                         &
    ! Tracer mass associated with liquid store.
  win_wtrac(:),                                                              &
    ! Tracer mass in the in-transit water (same units as tracer store).
  ice_lens_depth(:)
    ! Depth of uppermost lens (m)

!-----------------------------------------------------------------------------
! Local indices
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  n,                                                                           &
    ! Layer index, top-down.
  m,                                                                           &
    ! Reverse index used for upward redistribution during perching.
  i_wt,                                                                        &
    ! Water tracer index.
  n_lens
    ! index of first lens (0 = none yet)

!-----------------------------------------------------------------------------
! Local logicals
!-----------------------------------------------------------------------------
LOGICAL :: ice_lens
  ! True if the *current layer* satisfies the instantaneous pore-closure
  ! criterion (sice >= rho_firn_pore_closure * ds).
  ! This indicates that a lens is *forming here now*.
  !
  ! NOTE: This is a local, diagnostic flag only.
  ! It does NOT represent a persistent lens state.

LOGICAL :: blocked
  ! True once an effective impermeable barrier has been encountered above
  ! (lens depth or saturation barrier).
  !
  ! When .TRUE., all downward hydraulic processes are disabled for deeper
  ! layers in this timestep (no percolation, no capillary drainage).
  ! Refreezing may still occur below.

LOGICAL :: saturated_here
  ! True if the liquid water content in the current layer exceeds its
  ! full pore-space capacity:
  !   sliq > cap_full = porosity * rho_water * ds
  !
  ! When .TRUE., the layer is hydraulically saturated and any excess
  ! liquid must be redistributed upward (MONARCHS-style saturation handling).

LOGICAL :: at_lens_depth
  ! True if the *persistent uppermost ice lens depth* (ice_lens_depth)
  ! lies within the vertical extent of the current layer:
  !   z_top <= ice_lens_depth_m < z_bot
  !
  ! This flag enforces the impermeable barrier even if the instantaneous
  ! pore-closure criterion is not met in this layer due to regridding. 

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  coldsnow,                                                                    &
    ! Cold content relative to melting point (J m-2).
    ! Positive means layer is below tm and can refreeze liquid.
  dsice,                                                                       &
    ! Mass of liquid refrozen into ice this timestep in a layer (kg m-2).
  pfrac,                                                                       &
    ! Pore space fraction (dimensionless). Approximated here as 1 - ice volume fraction.
  cap_full,                                                                    &
    ! Pore space available for refreezing/percolation (kg m-2) 
  cap_mass,                                                                    &
    ! Capillary retained liquid mass in layer (kg m-2) 5% cap_full
  excess,                                                                      &
    ! Liquid mass (kg m-2) in excess of capillary retention, available to percolate.
  w_up,                                                                        &
    ! Perched mass (kg m-2) that must be stored above a lens.
  pfrac_m,                                                                     &
    ! Pore fraction in layer m when redistributing upward.
  space,                                                                       &
    ! Additional capacity (kg m-2) available in layer m pore space.
  addm,                                                                        &
    ! Amount (kg m-2) added to layer m during upward redistribution.
  ratio_sliq,                                                                  &
    ! Tracer-to-liquid ratio used for refreezing bookkeeping.
  wout,                                                                        &
    ! liquid water out after percolation and refreezing (kgm-2)
  z_top,                                                                       &
    ! depth from surface to top of current layer (m)
  z_bot,                                                                       &
   ! depth from surface to bottom of current layer (m)
  z_mid
   ! depth from surface to bottom of current layer (m)
 ! z_lens
    ! depth of upper most lens (m)



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PERCOLATE_MONARCHS'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

blocked = .FALSE.
wout    = 0.0

! running depth to top of current layer (m)
z_top   = 0.0

DO n = 1, nsnow

  z_bot = z_top + ds(i,n)
  z_mid = z_top + 0.5 * ds(i,n)

  !---------------------------------------------------------------------------
  ! 1) Apply incoming in-transit liquid to this layer ONLY if not blocked
  !---------------------------------------------------------------------------
  IF (.NOT. blocked) THEN
    sliq(i,n) = sliq(i,n) + win
    win = 0.0

    IF (l_wtrac_jls) THEN
      DO i_wt = 1, n_wtrac_jls
        sliq_wtrac(i,n,i_wt) = sliq_wtrac(i,n,i_wt) + win_wtrac(i_wt)
        win_wtrac(i_wt) = 0.0
      END DO
    END IF
  END IF

  !---------------------------------------------------------------------------
  ! 2) Refreezing based on cold content (ALWAYS, even below a lens)
  !---------------------------------------------------------------------------
  coldsnow = csnow(i,n) * (tm - tsnow(i,n))

  IF (coldsnow > 0.0 .AND. sliq(i,n) > 0.0) THEN

    dsice = MIN(sliq(i,n), coldsnow / lf)

    IF (l_wtrac_jls) THEN
      DO i_wt = 1, n_wtrac_jls
        IF (dsice == sliq(i,n)) THEN
          sice_wtrac(i,n,i_wt) = sice_wtrac(i,n,i_wt) + sliq_wtrac(i,n,i_wt)
          sliq_wtrac(i,n,i_wt) = 0.0
        ELSE
          ratio_sliq = wtrac_calc_ratio_fn_jules(i_wt, sliq_wtrac(i,n,i_wt), sliq(i,n))
          sice_wtrac(i,n,i_wt) = sice_wtrac(i,n,i_wt) + ratio_sliq * dsice
          sliq_wtrac(i,n,i_wt) = sliq_wtrac(i,n,i_wt) - ratio_sliq * dsice
        END IF
      END DO
    END IF

    sliq(i,n)  = sliq(i,n)  - dsice
    sice(i,n)  = sice(i,n)  + dsice
    tsnow(i,n) = tsnow(i,n) + lf * dsice / csnow(i,n)

    IF (sf_diag%l_snice) THEN
      sf_diag%snice_freez_surft(i,surft_n) = sf_diag%snice_freez_surft(i,surft_n) &
           + dsice / timestep
    END IF

  END IF

  !---------------------------------------------------------------------------
  ! If already blocked by an upstream lens, skip hydraulics but keep refreezing
  ! percolation has already been handled for this timestep above the lens,
  ! so there is no need to do it again, continue with refreezing only
  !---------------------------------------------------------------------------
  IF (blocked) THEN
     z_top = z_bot
    CYCLE
  END IF

  !-----------------------------------------------------------------------
  ! pore fraction and full pore capacity (MONARCHS saturation uses full cap)
  !-----------------------------------------------------------------------
  pfrac = 1.0 - sice(i,n) / (rho_ice * ds(i,n))
  pfrac = MAX(0.0, MIN(1.0, pfrac))

  cap_full = pfrac * rho_water * ds(i,n) ! FULL pore capacity (kg m-2)

  !---------------------------------------------------------------------------
  ! Instantaneous lens formation test (pore-closure density)
  !---------------------------------------------------------------------------
  ice_lens = .FALSE.
  IF ( sice(i,n) >= rho_firn_pore_closure * ds(i,n) ) ice_lens = .TRUE.

  ! Update persistent uppermost lens depth (MONARCHS: min old/new)
  IF (ice_lens) THEN
     IF (ice_lens_depth(i) < 0.0) THEN
        ice_lens_depth(i) = z_mid
     ELSE
        ice_lens_depth(i) = MIN(ice_lens_depth(i), z_mid)
     END IF
     !print *, 'ice lens depth', ice_lens_depth(i)
  END IF
 
  !---------------------------------------------------------------------------
  ! Are we at the persistent lens depth?
  !---------------------------------------------------------------------------
  at_lens_depth = .FALSE.
  IF (ice_lens_depth(i) >= 0.0) THEN
    IF (ice_lens_depth(i) >= z_top .AND. ice_lens_depth(i) < z_bot) at_lens_depth = .TRUE.
  END IF

  !---------------------------------------------------------------------------
  ! Saturation trigger (MONARCHS: if L exceeds pore space, push excess upward)
  !---------------------------------------------------------------------------
  saturated_here = (sliq(i,n) > cap_full)

  !---------------------------------------------------------------------------
  ! If lens exists here OR we are at stored lens depth OR saturated, do upward
  ! redistribution and then block deeper hydraulics.
  !---------------------------------------------------------------------------
  IF (ice_lens .OR. at_lens_depth .OR. saturated_here) THEN

     ! Use full pore capacity at this layer (MONARCHS saturation cap)
     cap_mass = cap_full

     w_up      = MAX(0.0, sliq(i,n) - cap_mass)
     sliq(i,n) = MIN(sliq(i,n), cap_mass)

     DO m = n-1, 1, -1
        pfrac_m = 1.0 - sice(i,m) / (rho_ice * ds(i,m))
        pfrac_m = MAX(0.0, MIN(1.0, pfrac_m))

        space = pfrac_m * rho_water * ds(i,m) - sliq(i,m)
        space = MAX(0.0, space)

        addm = MIN(space, w_up)
        sliq(i,m) = sliq(i,m) + addm
        w_up = w_up - addm

        IF (w_up <= 0.0) EXIT
     END DO

     ! leftover after filling upward
     wout = wout + w_up

     !PRINT *, 'BLOCK TRIGGER at n=', n, ' z_top=', z_top, ' z_mid=', z_mid, ' z_bot=', z_bot, &
     !    ' ice_lens=', ice_lens, ' at_lens_depth=', at_lens_depth, ' saturated=', saturated_here, &
     !    ' stored_ice_lens_depth=', ice_lens_depth(i)

     blocked = .TRUE.

     win = 0.0
     IF (l_wtrac_jls) THEN
       DO i_wt = 1, n_wtrac_jls
         win_wtrac(i_wt) = 0.0
       END DO
     END IF

  ELSE
     !-----------------------------------------------------------------------
     ! No lens and not saturated: downward percolation with 5% capillary retain
     !-----------------------------------------------------------------------
     cap_mass = 0.05 * cap_full

     excess    = MAX(0.0, sliq(i,n) - cap_mass)
     sliq(i,n) = sliq(i,n) - excess

     win = win + excess

     ! NOTE: tracer transport with excess is still missing (as in your original)
  END IF

  z_top = z_bot

END DO

! As in your original, add leftover after upward fill back into win
win = win + wout


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE percolate_monarchs
END MODULE percolate_monarchs_mod
