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
!   Conceptual picture:
!     - win is a scalar "in-transit" liquid mass (kg m-2) moving downward.
!     - Each layer can refreeze some liquid depending on cold content.
!     - Each layer retains a small amount of liquid by capillary forces
!       (fraction of pore space, MONARCHS-like).
!     - Any remaining excess attempts to percolate downward as win.
!     - If an ice lens is encountered (pore closure), excess is perched and
!       redistributed upward into available pore space above (saturation-upward).
!     - Any liquid that cannot be accommodated above the lens exits as runoff
!       (returned in win to the caller.
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
                               sf_diag, surft_n )


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
  win_wtrac(:)
    ! Tracer mass in the in-transit water (same units as tracer store).

!-----------------------------------------------------------------------------
! Local indices
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  n,                                                                           &
    ! Layer index, top-down.
  m,                                                                           &
    ! Reverse index used for upward redistribution during perching.
  i_wt
    ! Water tracer index.

!-----------------------------------------------------------------------------
! Local logicals
!-----------------------------------------------------------------------------
LOGICAL :: ice_lens
  ! True when the current layer is treated as impermeable (pore-closed).

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
  cap_mass,                                                                    &
    ! Capillary retained liquid mass in layer (kg m-2), proportional to pore space.
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
  ratio_sliq
    ! Tracer-to-liquid ratio used for refreezing bookkeeping.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PERCOLATE_MONARCHS'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Main loop over snow layers (top to bottom).
! This is a single-pass "streaming" algorithm:
!   win enters the layer, is absorbed/refrozen/retained, and any remainder is
!   passed downward as a new win.
!-----------------------------------------------------------------------------

DO n = 1, nsnow

    !-----------------------------------------------------------------------
    ! 1) Apply incoming in-transit liquid to this layer.
    !-----------------------------------------------------------------------
    sliq(i,n) = sliq(i,n) + win
    win = 0.0

    ! Tracer bookkeeping: incoming tracer mass is added to the layer liquid,
    ! and the in-transit tracer store is reset to zero for this layer.
    IF (l_wtrac_jls) THEN
      DO i_wt = 1, n_wtrac_jls
        sliq_wtrac(i,n,i_wt) = sliq_wtrac(i,n,i_wt) + win_wtrac(i_wt)
        win_wtrac(i_wt) = 0.0
      END DO
    END IF

    !-----------------------------------------------------------------------
    ! 2) Refreezing based on cold content.
    !
    ! coldsnow = csnow*(tm - T). If positive, layer is below melting point.
    ! Maximum refreezing is limited by:
    !   - available liquid in the layer
    !   - available energy to warm the layer to tm (cold content / lf)
    !-----------------------------------------------------------------------
    coldsnow = csnow(i,n) * (tm - tsnow(i,n))
    IF (coldsnow > 0.0 .AND. sliq(i,n) > 0.0) THEN
       dsice = MIN(sliq(i,n), coldsnow / lf)
       
       IF (sliq(i,n) < coldsnow / lf) THEN
       
          print *, n, 'n'
          print *, 'coldsnow', coldsnow
          print *, 'dsice', dsice
          print *, 'sliq', sliq(i,n)
          print *, 'tsnow_old', tsnow(i,n) - 273.15
          print *, 'deltaT', lf * dsice / csnow(i,n)
          print *, 'tsnow_new', tsnow(i,n) + lf * dsice / csnow(i,n) - 273.15
       END IF
       
      ! Tracer bookkeeping for freezing:
      !   - If all liquid freezes, move all tracer from liquid store to ice store.
      !   - If partial freezing, move tracer in proportion to the liquid that froze.
      IF (l_wtrac_jls) THEN
        DO i_wt = 1, n_wtrac_jls
          IF (dsice == sliq(i,n)) THEN
            ! All liquid freezes, so transfer all associated tracer.
            sice_wtrac(i,n,i_wt) = sice_wtrac(i,n,i_wt) + sliq_wtrac(i,n,i_wt)
            sliq_wtrac(i,n,i_wt) = 0.0
          ELSE
            ! Partial freezing: move tracer according to tracer/liquid ratio.
            ratio_sliq = wtrac_calc_ratio_fn_jules(i_wt, sliq_wtrac(i,n,i_wt), sliq(i,n))
            sice_wtrac(i,n,i_wt) = sice_wtrac(i,n,i_wt) + ratio_sliq * dsice
            sliq_wtrac(i,n,i_wt) = sliq_wtrac(i,n,i_wt) - ratio_sliq * dsice
          END IF
        END DO
      END IF

      ! Update mass stores and temperature due to release of latent heat.
      sliq(i,n)  = sliq(i,n)  - dsice
      sice(i,n)  = sice(i,n)  + dsice
      tsnow(i,n) = tsnow(i,n) + lf * dsice / csnow(i,n)

      ! Diagnostic: freezing rate (kg m-2 s-1) on this tile and point.
      IF (sf_diag%l_snice) THEN
        sf_diag%snice_freez_surft(i,surft_n) = sf_diag%snice_freez_surft(i,surft_n) &
                                                 + dsice / timestep
      END IF
    END IF

    !-----------------------------------------------------------------------
    ! 3) Capillary retention (MONARCHS-like).
    !
    ! Compute pore fraction from ice volume fraction:
    ! this is the total volume availble for refreezing (it can be liquid + air)
    !   ice volume fraction = sice / (rho_ice * ds) 
    !   pore fraction pfrac = 1 - ice volume fraction
    !
    ! Capillary retained liquid mass is a fixed fraction of pore-space water.
    ! Here: retain 5% of pore space expressed as an equivalent mass per area.
    ! Ligtenberg et al. (2011), The Cryosphere, doi:10.5194/tc-5-809-2011
    ! paper says values can be 4-13%.  
    !-----------------------------------------------------------------------
    pfrac = 1.0 - sice(i,n) / (rho_ice * ds(i,n)) !<-- Lfrac_max = 1 - cell["Sfrac"][v_lev] in calc_saturation
    pfrac = MAX(0.0, MIN(1.0, pfrac))
 
    cap_mass = 0.05 * pfrac * rho_water * ds(i,n) !<-- capillary_remain = 0.05 * (1 - cell["Sfrac"][v_lev]) 

    ! Liquid above capillary retention becomes mobile and may percolate.
    excess = MAX(0.0, sliq(i,n) - cap_mass)
    sliq(i,n) = sliq(i,n) - excess

    !-----------------------------------------------------------------------
    ! 4) Ice lens test (impermeable barrier).
    !
    ! A lens is diagnosed when the layer ice density reaches pore-closure
    ! threshold: sice/ds >= rho_firn_pore_closure.
    !
    ! Implemented as sice >= rho_close * ds to avoid division in comments,
    ! but this code uses the multiplication form already.
    !-----------------------------------------------------------------------
    ice_lens = .FALSE.
    !IF ( sice(n) >= rho_firn_pore_closure * ds(n) ) ice_lens = .TRUE.
    IF ( sice(i,n) >= 730.0 * ds(i,n) ) ice_lens = .TRUE. !<-- if cell["Sfrac"][v_lev] * cell["rho_ice"] > cell["pore_closure"]:

    !-----------------------------------------------------------------------
    ! 5) Perching and saturation-upward if blocked by a lens.
    !
    ! If the lens is present, "excess" cannot continue downward.
    ! It is redistributed upward, filling available pore-space capacity
    ! above the lens, starting from the lens layer and moving upward.
    !
    ! Any liquid that still cannot be stored above becomes outflow (win),
    ! which the caller will interpret as runoff / drainage from the column.
    !-----------------------------------------------------------------------
    IF (ice_lens .AND. excess > 0.0) THEN
      w_up = excess ! w_up is water that must go up because it cant go down

      DO m = n, 1, -1
        IF (ds(i,m) <= EPSILON(ds(i,m))) CYCLE

        ! Available pore space in layer m, expressed as a maximum liquid mass.
        pfrac_m = 1.0 - sice(i,m) / (rho_ice * ds(i,m))
        pfrac_m = MAX(0.0, MIN(1.0, pfrac_m))

        ! space is how much additional liquid can be stored before saturation.
        ! space = (maximum pore storage) − (liquid already in pores)
        space = pfrac_m * rho_water * ds(i,m) - sliq(i,m)
        space = MAX(0.0, space)

        ! Fill as much of that space as possible with perched water.
        addm = MIN(space, w_up)
        sliq(i,m) = sliq(i,m) + addm
        w_up = w_up - addm
        IF (w_up <= 0.0) EXIT
      END DO

      ! Any remaining perched water after saturating upward leaves the column.
      win = win + w_up

      ! Stop processing deeper layers because flow is blocked below the lens.
      EXIT

    ELSE
      ! No lens barrier, pass mobile liquid to the next layer down.
      win = win + excess
    END IF

  END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE percolate_monarchs
END MODULE percolate_monarchs_mod
