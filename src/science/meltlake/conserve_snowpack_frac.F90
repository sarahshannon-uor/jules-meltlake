
! ****************************COPYRIGHT*******************************
!
! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]
!
! *****************************COPYRIGHT*******************************
! SUBROUTINE CONSERVE_SNOWPACK_FRAC----------------------------------
!
! Description:
!     Adjust snowpack variables following a change in lake tile
!     fraction in order to conserve gridbox snow mass.
!
! Method:
!     Apply the snow mass correction to the lowest snow layer on
!     the tile gaining fraction. The ice lens and snow layers above
!     the lowest layer are unchanged.
!
! Code Owner: s.r.shannon@reading.ac.uk
!
! Subroutine Interface:
MODULE conserve_snowpack_frac_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                      &
  ModuleName='CONSERVE_SNOWPACK_FRAC_MOD'

CONTAINS

SUBROUTINE conserve_snowpack_frac(land_pts,           &
                                  nsurft,             &
                                  surft_pts,          &
                                  surft_index,        &
                                  frac_lake_old,      &
                                  frac,               &
                                  nsnow,              &
                                  ice_lens_index,     &
                                  sice_ml,            &
                                  sliq_ml,            &
                                  tsnow_ml,           &
                                  ds_ml,              &
                                  snowdepth,          &
                                  snow_surft,         &
                                  rho_snow_ml,        &
                                  rho_snow_grnd)

USE jules_meltlake_mod, ONLY: nsmax_ml

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                       &
  land_pts,                                                                  &
    ! Number of land points.
  nsurft,                                                                    &
    ! Number of surface tiles.
  surft_pts
    ! Number of lake tile points.


!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                       &
  surft_index(land_pts)
    ! Index of lake tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
  frac_lake_old(land_pts),                                                   &
    ! Lake tile fraction before area adjustment.
  frac(land_pts,nsurft),                                                     &
    ! Fractional coverage of each surface tile after area adjustment.
  ice_lens_index(land_pts,nsurft)
    ! Snow layer containing the ice lens.


!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN OUT) ::                                                   &
  nsnow(land_pts,nsurft)
    ! Number of snow layers.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  sice_ml(land_pts,nsurft,nsmax_ml),                                         &
    ! Snow layer ice mass (kg m-2).
  sliq_ml(land_pts,nsurft,nsmax_ml),                                         &
    ! Snow layer liquid mass (kg m-2).
  tsnow_ml(land_pts,nsurft,nsmax_ml),                                        &
    ! Snow layer temperatures (K).
  ds_ml(land_pts,nsurft,nsmax_ml),                                           &
    ! Snow layer thicknesses (m).
  snowdepth(land_pts,nsurft),                                                &
    ! Snow depth (m).
  snow_surft(land_pts,nsurft),                                               &
    ! Snow mass on tiles (kg m-2).
  rho_snow_ml(land_pts,nsurft,nsmax_ml),                                     &
    ! Snow layer densities (kg m-3).
  rho_snow_grnd(land_pts,nsurft)
    ! Snowpack bulk density (kg m-3).


!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------

INTEGER ::                                                                   &
  i,                                                                         &
    ! Land-point index.
  k,                                                                         &
    ! Lake tile point index.
  ns,                                                                        &
    ! Snow layer index.
  tile_adjust
    ! Tile gaining fraction.

INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
    ! Melt-lake tile index.
  nonlake_ml = 10
    ! Non-lake elevated-ice tile index.

REAL(KIND=real_jlslsm) ::                                                    &
  delta_frac,                                                                &
    ! Change in lake tile fraction.
  snowmass_lake,                                                             &
    ! Snow mass on lake tile (kg m-2).
  snowmass_nonlake,                                                          &
    ! Snow mass on non-lake tile (kg m-2).
  delta_snowmass,                                                            &
    ! Snow mass correction on tile gaining fraction (kg m-2).
  snowmass_bottom,                                                           &
    ! Snow mass in lowest snow layer (kg m-2).
  scale_mass
    ! Factor used to adjust lowest snow layer.


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CONSERVE_SNOWPACK_FRAC'


!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_in,zhook_handle)


!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,ns,delta_frac,snowmass_lake,snowmass_nonlake,       &
!$OMP         delta_snowmass,tile_adjust,snowmass_bottom,scale_mass)

DO k = 1, surft_pts

   i = surft_index(k)

!-----------------------------------------------------------------------------
! Calculate change in lake tile fraction.
!-----------------------------------------------------------------------------

   delta_frac = frac(i,lake_ml) - frac_lake_old(i)


   IF (delta_frac /= 0.0) THEN

!-----------------------------------------------------------------------------
! Calculate snow mass on lake and non-lake tiles.
!-----------------------------------------------------------------------------

      snowmass_lake = SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml))) &
           + SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml)))

      snowmass_nonlake = SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) &
           + SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))


!-----------------------------------------------------------------------------
! Lake tile gains fraction. Non-lake tile loses fraction.
!-----------------------------------------------------------------------------

      IF (delta_frac > 0.0) THEN

         tile_adjust = lake_ml

         IF (frac(i,lake_ml) <= 0.0) THEN
            WRITE(*,*) 'ERROR: lake fraction must be positive'
            STOP
         END IF

         delta_snowmass = delta_frac / frac(i,lake_ml) &
              * (snowmass_nonlake - snowmass_lake)


!-----------------------------------------------------------------------------
! Lake tile loses fraction. Non-lake tile gains fraction.
!-----------------------------------------------------------------------------

      ELSE IF (delta_frac < 0.0) THEN

         tile_adjust = nonlake_ml

         IF (frac(i,nonlake_ml) <= 0.0) THEN
            WRITE(*,*) 'ERROR: non-lake fraction must be positive'
            STOP
         END IF

         delta_snowmass = (-delta_frac) / frac(i,nonlake_ml) &
              * (snowmass_lake - snowmass_nonlake)

      END IF


      WRITE(*,*) '----------------------------------------'
      WRITE(*,*) 'delta_frac = ', delta_frac
      WRITE(*,*) 'tile_adjust = ', tile_adjust
      WRITE(*,*) 'delta_snowmass = ', delta_snowmass

      WRITE(*,*) 'gridbox snow mass before = ',                     &
           frac_lake_old(i) * snowmass_lake                         &
           + (1.0 - frac_lake_old(i)) * snowmass_nonlake


!-----------------------------------------------------------------------------
! Apply snow mass correction to the lowest snow layer.
!-----------------------------------------------------------------------------

      ns = nsnow(i,tile_adjust)

      IF (ns <= 0) THEN
         WRITE(*,*) 'ERROR: no snow layers on tile ', tile_adjust
         STOP
      END IF

      snowmass_bottom = sice_ml(i,tile_adjust,ns) &
           + sliq_ml(i,tile_adjust,ns)


!-----------------------------------------------------------------------------
! Check that lowest layer has positive mass and sufficient mass
! remains after the correction.
!-----------------------------------------------------------------------------

      IF (snowmass_bottom <= 0.0) THEN

         WRITE(*,*) 'ERROR: bottom snow layer has zero mass'
         WRITE(*,*) 'tile_adjust = ', tile_adjust
         WRITE(*,*) 'snowmass_bottom = ', snowmass_bottom
         STOP

      END IF

      IF (snowmass_bottom + delta_snowmass <= 0.0) THEN

         WRITE(*,*) 'ERROR: insufficient mass in lowest snow layer'
         WRITE(*,*) 'tile_adjust = ', tile_adjust
         WRITE(*,*) 'snowmass_bottom = ', snowmass_bottom
         WRITE(*,*) 'delta_snowmass = ', delta_snowmass
         STOP

      END IF


!-----------------------------------------------------------------------------
! Scale ice mass, liquid mass and thickness by the same factor so that
! the density and ice/liquid fractions of the lowest layer are unchanged.
!-----------------------------------------------------------------------------

      WRITE(*,*) 'bottom layer before correction'
      WRITE(*,*) 'sice_ml = ', sice_ml(i,tile_adjust,ns)
      WRITE(*,*) 'sliq_ml = ', sliq_ml(i,tile_adjust,ns)
      WRITE(*,*) 'ds_ml   = ', ds_ml(i,tile_adjust,ns)

      scale_mass = (snowmass_bottom + delta_snowmass) &
           / snowmass_bottom

      sice_ml(i,tile_adjust,ns) = &
           sice_ml(i,tile_adjust,ns) * scale_mass

      sliq_ml(i,tile_adjust,ns) = &
           sliq_ml(i,tile_adjust,ns) * scale_mass

      ds_ml(i,tile_adjust,ns) = &
           ds_ml(i,tile_adjust,ns) * scale_mass

      WRITE(*,*) 'bottom layer after correction'
      WRITE(*,*) 'sice_ml = ', sice_ml(i,tile_adjust,ns)
      WRITE(*,*) 'sliq_ml = ', sliq_ml(i,tile_adjust,ns)
      WRITE(*,*) 'ds_ml   = ', ds_ml(i,tile_adjust,ns)


!-----------------------------------------------------------------------------
! Update total snow mass, snow depth and bulk snow density.
!-----------------------------------------------------------------------------

      snow_surft(i,tile_adjust) =                                  &
           SUM(sice_ml(i,tile_adjust,1:nsnow(i,tile_adjust)))       &
         + SUM(sliq_ml(i,tile_adjust,1:nsnow(i,tile_adjust)))

      snowdepth(i,tile_adjust) =                                   &
           SUM(ds_ml(i,tile_adjust,1:nsnow(i,tile_adjust)))

      IF (snowdepth(i,tile_adjust) > 0.0) THEN

         rho_snow_grnd(i,tile_adjust) =                            &
              snow_surft(i,tile_adjust) / snowdepth(i,tile_adjust)

      END IF


!-----------------------------------------------------------------------------
! Check gridbox snow mass conservation.
!-----------------------------------------------------------------------------

      WRITE(*,*) 'gridbox snow mass after  = ',                      &
           frac(i,lake_ml)                                           &
           * (SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml)))             &
           +  SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml))))            &
           + frac(i,nonlake_ml)                                      &
           * (SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))       &
           +  SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))))

      WRITE(*,*) 'mass difference = ',                               &
           frac(i,lake_ml)                                           &
           * (SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml)))             &
           +  SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml))))            &
           + frac(i,nonlake_ml)                                      &
           * (SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))       &
           +  SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))))      &
           - (frac_lake_old(i) * snowmass_lake                       &
           + (1.0 - frac_lake_old(i)) * snowmass_nonlake)


      IF (ABS(                                                        &
           frac(i,lake_ml)                                            &
           * (SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml)))              &
           +  SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml))))             &
           + frac(i,nonlake_ml)                                       &
           * (SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))        &
           +  SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))))       &
           - (frac_lake_old(i) * snowmass_lake                        &
           + (1.0 - frac_lake_old(i)) * snowmass_nonlake))            &
           > 0.01) THEN

         WRITE(*,*) 'ERROR: snow mass not conserved after lake frac adjustment'
         STOP

      END IF

   END IF

END DO

!$OMP END PARALLEL DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE conserve_snowpack_frac

END MODULE conserve_snowpack_frac_mod
