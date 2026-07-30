! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE apply_lake_bottom_melt-----------------------------------

! Description:
!     

! Subroutine Interface:
MODULE apply_lake_bottom_melt_mod
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='APPLY_LAKE_BOTTOM_MELT'

CONTAINS

  SUBROUTINE apply_lake_bottom_melt(land_pts,                 & ! IN
       surft_pts,                & ! IN
       surft_index,              & ! IN
       nsnow,                    & ! IN
       dhdt_lake_snow_ml,        & ! IN
       ds_ml,                    & ! IN/OUT
       sice_ml,                  & ! IN/OUT
       sliq_ml,                  & ! IN/OUT
       snowmass,                 & ! IN/OUT
       snowdepth,                & ! IN/OUT
       lake_depth_ml)              ! IN/OUT

USE water_constants_mod, ONLY:                            &
     rho_water
      ! Density of pure water (kg/m3).

USE jules_meltlake_mod, ONLY: nsmax_ml

USE um_types, ONLY: real_jlslsm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook


IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                    &
     land_pts,                                               &
     ! Total number of land points.
     surft_pts
    ! Number of tile points.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                       &
     surft_index(land_pts),                                                  &
     ! Index of tile points.
     nsnow(land_pts)
     ! Number of snow layers.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
     dhdt_lake_snow_ml(land_pts)
        ! Stefan boundary movement for this timestep (m).

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
     ds_ml(land_pts,nsmax_ml),                                               &
        ! Snow layer thicknesses (m).
     sice_ml(land_pts,nsmax_ml),                                             &
        ! Ice content of snow layers (kg/m2).
     sliq_ml(land_pts,nsmax_ml),                                             &
        ! Liquid content of snow layers (kg/m2).
     snowmass(land_pts),                                                     &
        ! Snow mass on current tile (kg/m2).
     snowdepth(land_pts),                                                    &
        ! Snow depth on current tile (m).
     lake_depth_ml(land_pts)
        ! Lake depth (m).

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                &
     i,                                                      &
        ! Land point index.
     k,                                                      &
        ! Tile point index.
     n
        ! Snow layer index.

REAL(KIND=real_jlslsm) ::                                    &
     dh_ice,                                                 &
        ! Boundary change Stefan condition (m per timestep).
     dh_remain,                                              &
        ! Remaining boundary retreat to apply (m).
     frac_melt,                                              &
        ! Fraction of current snow layer removed.
     dsice,                                                  &
        ! Ice removed from current layer (kg/m2).
     dsliq,                                                  &
        ! Liquid removed from current layer (kg/m2).
     ds_old,                                                 &
        ! Snow layer thickness before melt (m).
     dsice_tot,                                              &
        ! Total ice removed from snowpack (kg/m2).
     dsliq_tot
        ! Total liquid removed from snowpack (kg/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='APPLY_LAKE_BOTTOM_MELT'

!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Apply lake-bottom Stefan melt to snow mass and geometry.
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)               &
!$OMP PRIVATE(k,i,n,dh_ice,dh_remain,frac_melt,dsice,dsliq,    &
!$OMP         ds_old,dsice_tot,dsliq_tot)                      &
!$OMP SHARED(surft_pts,surft_index,nsnow,dhdt_lake_snow_ml,    &
!$OMP        ds_ml,sice_ml,sliq_ml,snowmass,snowdepth,         &
!$OMP        lake_depth_ml)
DO k = 1, surft_pts
   i = surft_index(k)

   
   dh_ice = dhdt_lake_snow_ml(i)

   IF (dh_ice > 0.0 .AND. nsnow(i) > 0) THEN

      dsice_tot = 0.0
      dsliq_tot = 0.0
      dh_remain = dh_ice

!-----------------------------------------------------------------------------
! Stefan retreat will normally be small compared with ds(1), but loop through
! levels in case dh_ice melts multiple snow layers. Removal is based on the
! fraction of snow layer depth that is retreating.
!-----------------------------------------------------------------------------
      DO n = 1, nsnow(i)

         IF (dh_remain <= 0.0) EXIT

         ds_old = ds_ml(i,n)

         IF (ds_old <= 0.0) CYCLE

         frac_melt = MIN(1.0, dh_remain / ds_old)

         dsice = frac_melt * sice_ml(i,n)
         dsliq = frac_melt * sliq_ml(i,n)

         dsice_tot = dsice_tot + dsice
         dsliq_tot = dsliq_tot + dsliq

!-----------------------------------------------------------------------------
! Adjust sice, sliq and ds.
!-----------------------------------------------------------------------------
         sice_ml(i,n) = sice_ml(i,n) - dsice
         sliq_ml(i,n) = sliq_ml(i,n) - dsliq
         ds_ml(i,n)   = ds_old * (1.0 - frac_melt)

         dh_remain = dh_remain - frac_melt * ds_old

      END DO

!-----------------------------------------------------------------------------
! Subtract snow mass and snow depth removed by Stefan boundary retreat.
!-----------------------------------------------------------------------------
      snowmass(i) = snowmass(i) - dsice_tot - dsliq_tot
      snowdepth(i) = snowdepth(i) - (dh_ice - dh_remain)

!-----------------------------------------------------------------------------
! Add melted ice and drained liquid to the lake depth.
!-----------------------------------------------------------------------------
      lake_depth_ml(i) = lake_depth_ml(i) +                 &
           (dsice_tot + dsliq_tot) / rho_water

   END IF

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
    
END SUBROUTINE apply_lake_bottom_melt

END MODULE apply_lake_bottom_melt_mod
