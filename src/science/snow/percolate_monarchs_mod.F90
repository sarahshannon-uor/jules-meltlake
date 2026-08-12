! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!-----------------------------------------------------------------------------
! Description:
!   Vertical meltwater percolation through the multilayer snow
!   column, including refreezing, capillary retention, ice lens
!   formation, and upward perching of water.
!
! Method:
!   The scheme follows a Buzzard-style percolation approach in which
!   meltwater is added at the surface and passed sequentially through
!   snow layers within a single timestep.
!
!   For each layer:
!     (1) Incoming liquid water is added.
!     (2) Refreezing occurs based on cold content (like default jules).
!     (3) If the layer is permeable:
!           - A small fraction of liquid water (5% of pore space) is
!             retained to represent capillary forces.
!           - The remaining water percolates downward to the next layer.
!     (4) If an impermeable ice lens is present:
!           - Downward percolation is blocked.
!           - Excess water is redistributed upward, filling overlying
!             layers to their pore capacity.
!      5) Excess water that remains after upward filling above an
!            impermeable ice lens is routed to lake_inflow and used to grow
!            the melt lake. Water that percolates out of the bottom of the
!            snowpack is also added to lake_inflow.
!
!
! Notes:
!   - Percolation is instantaneous within a timestep (no hydraulic delay).
!   - Full layer saturation occurs only in the presence of a lens.
!   - The subroutine loops top-down through the snow column (reverse of monarchs)
!
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
  
  SUBROUTINE percolate_monarchs( i, nsnow, land_pts, nsmax, n_wtrac_jls,     &
                               timestep, csnow, ds, tsnow,                   &
                               sice, sliq,                                   &
                               sice_wtrac, sliq_wtrac,                       &
                               win, win_wtrac,                               &
                               sf_diag, surft_n, refreeze, ice_lens_depth,   &
                               ice_lens_index, lake_inflow, has_lake)


!-----------------------------------------------------------------------------
! Inputs and supporting modules
!-----------------------------------------------------------------------------

USE jules_water_tracers_mod, ONLY: l_wtrac_jls, wtrac_calc_ratio_fn_jules
  ! l_wtrac_jls enables tracer bookkeeping.
  ! wtrac_calc_ratio_fn_jules provides tracer-to-water ratio in a store.

USE jules_snow_mod,          ONLY: rho_firn_pore_closure
  ! Density threshold at which pores are assumed closed (kg m-3).
  ! Used for the ice lens criterion.

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
USE model_time_mod, ONLY: timestep_number

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
  land_pts,                                                                    &
    ! Total number of land points.
  surft_n,                                                                     &
    ! Surface tile index for diagnostics.
  n_wtrac_jls                                                                 
    ! Number of water tracers.
  
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
    ! Model timestep (s). Used to convert masses to flux diagnostics.
  csnow(land_pts,nsmax)
    ! Areal heat capacity of each snow layer (J K-1 m-2).

LOGICAL, INTENT(IN) ::                                                        &
   has_lake(land_pts)
    ! lake depth > 10cm
   
!-----------------------------------------------------------------------------
! Scalar arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  win
    ! Liquid water mass in transit downward (kg m-2). Water out bottom at end of 
    ! subroutine

TYPE (strnewsfdiag), INTENT(IN OUT) :: sf_diag
  ! Diagnostics structure, updated here for refreezing flux.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  ds(land_pts,nsmax),                                                        &
    ! Layer thicknesses (m).
  tsnow(land_pts,nsmax),                                                     &
    ! Layer temperatures (K).
  sice(land_pts,nsmax),                                                      &
    ! Ice mass per unit area in each layer (kg m-2).
  sliq(land_pts,nsmax),                                                      &
    ! Liquid mass per unit area in each layer (kg m-2).
  sice_wtrac(land_pts,nsmax,n_wtrac_jls),                                    &
    ! Tracer mass associated with ice store (kg m-2 tracer units).
  sliq_wtrac(land_pts,nsmax,n_wtrac_jls),                                    &
    ! Tracer mass associated with liquid store.
  win_wtrac(n_wtrac_jls),                                                    &
    ! Tracer mass in the in-transit water (same units as tracer store).
  ice_lens_depth(land_pts),                                                  &
    ! Depth of ice lens
  ice_lens_index(land_pts),                                                  &
    ! Snow level of uppermost lens (as a real) 
  lake_inflow(land_pts)
    ! water mass out of passed to lake (kg m-2)
 

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                       &
     refreeze(land_pts,nsmax)
    ! Mass that is refreezing in snow layer (kgm-2). Only for plotting 

!-----------------------------------------------------------------------------
! Local indices
!-----------------------------------------------------------------------------
INTEGER ::                                                                    &
  n,                                                                          &
    ! Layer index, top-down.
  m,                                                                          &
    ! Reverse index used for upward redistribution during perching.
  i_wt, iw
    ! Water tracer index.
 
    
!-----------------------------------------------------------------------------
! Local logicals
!-----------------------------------------------------------------------------
LOGICAL :: blocked
  ! True once downward percolation has been blocked by an ice lens
  ! during the current timestep. Deeper layers no longer receive
  ! liquid water from above, although refreezing can still happen.

LOGICAL :: lens_found
  ! True if a layer satisfies the pore closure criterion
  ! during the current timestep.
  !
  ! This is used to indicate that a valid ice lens layer has been
  ! identified during the present percolation calculation. It does
  ! not itself store the persistent ice lens state.

LOGICAL :: ice_lens
  ! Persistent ice lens state for the column.
  !
  ! True if an ice lens already existed on entry to the subroutine or
  ! if a new valid ice lens is identified during this timestep.


!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  coldsnow,                                                                    &
    ! Cold content relative to melting point (J m-2).
    ! Positive means layer is below tm and can refreeze liquid.
  dsice,                                                                       &
    ! Mass of liquid refrozen into ice this timestep in a layer (kg m-2).
  cap_full,                                                                    &
    ! Maximum liquid water capacity of the current layer (kg m-2) 
  cap_full_m,                                                                  &
   ! Maximum liquid water capacity of layers above the lens (kg m-2)  
  cap_mass,                                                                    &
    ! Capillary retained liquid mass in layer (kg m-2) 5% cap_full
  excess,                                                                      &
    ! Liquid mass (kg m-2) in excess of capillary retention, available to percolate.
  w_up,                                                                        &
    ! Perched mass (kg m-2) that must be stored above a lens.
  ratio_sliq,                                                                  &
    ! Tracer to liquid ratio used for refreezing bookkeeping.
  z_top,                                                                       &
    ! depth from surface to top of current layer (m)
  z_bot,                                                                       &
   ! depth from surface to bottom of current layer (m)
  z_mid,                                                                       &
   ! depth from surface to bottom of current layer (m)
  wout,                                                                        &
   ! liquid water mass out of top after upward fill (kgm-2)
  ice_capacity,                                                                &
   !maximum additional ice mass that can fit in the layer 
  win_in
   ! water into snowpack 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PERCOLATE_MONARCHS'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

refreeze(i,:) = 0.0
lake_inflow(i) = 0.0

blocked = .FALSE.
wout    = 0.0

win_in = win

! running depth to top of current layer (m)
z_top   = 0.0

IF (has_lake(i)) THEN
   win = 0
   !print *, 'lake present bypass percolate_monarchs'
   RETURN
END IF

DO n = 1, nsnow
   
   !---------------------------------------------------------------------------
   ! keep track of the snow depth. This is needed to track the lens if the
   ! snowdepth changes from compaction/melting/sublimination/accumulation
   !---------------------------------------------------------------------------
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
  ! 2) Refreezing based on cold content (all layers, even below a lens)
  !---------------------------------------------------------------------------
  coldsnow = csnow(i,n) * (tm - tsnow(i,n) + 1.0e-3)

   IF (coldsnow > 0.0 .AND. sliq(i,n) > 0.0) THEN

    ! --- maximum additional ice mass that can fit in the layer
    ice_capacity = rho_ice * ds(i,n) - sice(i,n)
    ice_capacity = MAX(0.0, ice_capacity)

    !--- don't allow refreezing to add more ice into a layer than it can hold
    !--- otherwise bulk density (rho_snow_ml) can go > 917 and Sfrac > 1
    dsice = MIN(sliq(i,n), coldsnow / lf, ice_capacity)

   
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

    
    refreeze(i,n) = dsice   ! (kg m-2) just for output netCDF 

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

  IF (ds(i,n) <= EPSILON(0.0)) THEN
     win       = win + sliq(i,n)
     sliq(i,n) = 0.0
     z_top     = z_bot
     CYCLE
  END IF

    
  ! Maximum water capacity of layer. Limited so the bulk density does not to go > 917
  cap_full = MAX(0.0, rho_ice * ds(i,n) - sice(i,n))
  
  !---------------------------------------------------------------------------
  ! Instantaneous lens formation test (pore closure density)
  !---------------------------------------------------------------------------
  lens_found = (sice(i,n) >= rho_firn_pore_closure * ds(i,n))

  ! Update persistent lens: keep existing, or replace if new is shallower
  IF (lens_found) THEN

     ! find depth of new lens 
     IF (ice_lens_depth(i) < 0.0 .OR. z_mid < ice_lens_depth(i)) THEN
        ice_lens_depth(i) = z_mid
        
     END IF
  END IF! lens_found


! find the snow layer that contains the lens
  ice_lens = .FALSE.

  IF (ice_lens_depth(i) >= 0.0) THEN
     ice_lens = (ice_lens_depth(i) >= z_top .AND. ice_lens_depth(i) < z_bot)
  END IF


  !-----------------------------------------------------------------------
  ! Remove a previously stored ice lens if no layer satisfies the
  ! ice lens criterion after the current percolation.
  !
  ! ice_lens indicates that a lens existed on entry to the routine,
  ! while found_lens is set true if a lens is found during
  ! this timestep.
  !-----------------------------------------------------------------------
  
  IF (ice_lens .AND. .NOT. lens_found) THEN
     ice_lens_depth(i) = -1.0
     ice_lens_index(i) = -1.0
     ice_lens = .FALSE.
  END IF


  
  !---------------------------------------------------------------------------
  ! 4) If lens exists here, then upward percoalte water into layer above 
  ! up to pore space capacity
  !-------------------------------------------------------------------------
    IF (ice_lens) THEN

       ! output snow layer of ice lens
       ice_lens_index(i) = REAL(n)
       
       blocked = .TRUE.

       ! make lens completely dry
       w_up      = sliq(i,n)
       sliq(i,n) = 0.0

             
     ! upward fill 
       DO m = n-1, 1, -1

          IF (ds(i,m) <= EPSILON(0.0)) CYCLE  ! no layer so nothing to fill

          !-----------------------------------------------------------------------
          ! Original MONARCHS liquid water capacity.
          ! Fill the available pore volume with water:
          !
          !   cap_full_m = pore fraction * water density * layer thickness
          !
          ! In the JULES layer representation this can produce a total layer
          ! mass corresponding to a bulk density greater than rho_ice, so this
          ! expression is not used here.
          !-----------------------------------------------------------------------
          ! cap_full_m = pfrac_m * rho_water * ds(i,m)
          !
          ! where 
          ! pfrac_m = 1.0 - sice(i,m) / (rho_ice * ds(i,m))
          ! pfrac_m = MAX(0.0, MIN(1.0, pfrac_m))
          !-----------------------------------------------------------------------
          ! Limit the liquid water capacity so that the total mass of ice plus
          ! liquid water does not exceed rho_ice * layer thickness.
          !
          ! cap_full_m is the additional liquid water mass that can
          ! be stored before the layer reaches an equivalent bulk density of
          ! rho_ice.
          !-----------------------------------------------------------------------
          
          cap_full_m = MAX(0.0, rho_ice * ds(i,m) - sice(i,m))
          
          ! --- water to upfill
          sliq(i,m) = sliq(i,m) + w_up
          
          
          ! --- only allow layers to contain water up to pore capacity 
          IF (sliq(i,m) > cap_full_m) THEN
             w_up = sliq(i,m) - cap_full_m
             sliq(i,m) = cap_full_m
          ELSE
             w_up = 0.0 ! no more water to upfill
          END IF

                   
          IF (w_up <= 0.0) EXIT

       END DO ! m loop

     
     ! leftover after filling upward goes out of top into lake
       wout = wout + w_up
          
     IF (l_wtrac_jls) THEN
        DO i_wt = 1, n_wtrac_jls
           win_wtrac(i_wt) = 0.0
        END DO
     END IF

  ELSE 

     !-----------------------------------------------------------------------
     ! 3) No lens: downward percolation with 5% capillary retained
     !-----------------------------------------------------------------------

     cap_mass = 0.05 * cap_full
     excess    = MAX(0.0, sliq(i,n) - cap_mass)
     sliq(i,n) = sliq(i,n) - excess
     win = win + excess

     ! NOTE: problably need to add tracer stuff
  END IF

  IF (sliq(i,n) > cap_full + 1.0e-12) THEN
     print *, 'liquid exceeds full pore capacity'
     stop
  END IF
  

  z_top = z_bot

  
END DO ! nsnow

!-----------------------------------------------------------------------
! 5) water exiting the snowpack (water out bottom [win] + out top [wout])
! Check with Robin. Should water out the bottom go into lake ?
!-----------------------------------------------------------------------

lake_inflow(i) = win + wout
win = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE percolate_monarchs


END MODULE percolate_monarchs_mod
