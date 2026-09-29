! *****************************COPYRIGHT*******************************

! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]

! *****************************COPYRIGHT*******************************
!  SUBROUTINE MELTLAKE_SURFACE_PRECIP-------------------------------

! Description:
!    Applies rainfall and snowfall to melt lake surface states.
!    Rain falling on a virtual or permanent lake lid is stored as
!    rain on the lid. Snowfall is added to the snow depth on the lid.
!    Precipitation falling on exposed lake water is added directly
!    to the lake as liquid water equivalent.
!
!    Snowfall on exposed water is treated as immediately melted.
!    Mass is conserved, but the latent heat required for melting is
!    not currently included.
! Code Owner: s.r.shannon@reading.ac.uk
!
! Subroutine Interface:
MODULE meltlake_surface_precip_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MELTLAKE_SURFACE_PRECIP_MOD'

CONTAINS

  SUBROUTINE meltlake_surface_precip( &
       land_pts,&
       timestep, &
       nsurft, &
       surft_pts,& 
       surft_index,& 
       ls_snow,&
       con_snow,&
       ls_rain,&
       con_rain,&                                  
       has_lid,&
       has_vlid,&
       exposed_water,&
       snow_on_lid,&                       
       lake_depth_ml,&
       lid_snow_depth_ml,&
       lid_rain_water_ml)


USE jules_snow_mod, ONLY: rho_snow_const
USE water_constants_mod, ONLY: rho_water

USE um_types, ONLY: real_jlslsm
  
    
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
     ! Total number of land points.
   nsurft
     ! Number of land tiles.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  surft_pts(nsurft),                                                           &
    ! Number of tile points.
  surft_index(land_pts,nsurft)
    ! Index of tile points.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep
   ! Timestep length (s).

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ls_snow(land_pts),                                                           &
    ! Large-scale snowfall rate (kg m-2 s-1)                              
  con_snow(land_pts),                                                          &
    ! Convective snowfall rate (kg m-2 s-1)                               
  ls_rain(land_pts),                                                           &
    ! Large-scale rainfall rate (kg m-2 s-1)                              
  con_rain(land_pts)
    ! Convective rainfall rate (kg m-2 s-1)

LOGICAL, INTENT(IN) ::                                                         &
  has_lid(land_pts,nsurft),                                                    &
    ! True where a permanent melt-lake lid is present.
  has_vlid(land_pts,nsurft),                                                   &
    ! True where a virtual melt-lake lid is present.
  exposed_water(land_pts,nsurft)
    ! True where melt-lake water is exposed at the surface.

!-----------------------------------------------------------------------------
! Array arguments with intent(inout)
!-----------------------------------------------------------------------------
LOGICAL, INTENT(IN OUT) ::                                                     &
  snow_on_lid(land_pts,nsurft)
    ! True where snow is present on a virtual or permanent lake lid.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  lake_depth_ml(land_pts,nsurft),                                              &
    ! Melt-lake water depth (m).
  lid_snow_depth_ml(land_pts,nsurft), &
    ! Snow depth on a virtual or permanent lake lid (m).
  lid_rain_water_ml(land_pts,nsurft)
    ! Liquid rainwater stored on the melt lake lid (m).
    

!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  i,                                                                           &
    ! Land point index.
  j,                                                                           &
    ! Tile point loop counter.
  n
    ! Surface tile loop counter.

REAL(KIND=real_jlslsm) ::                                                      &
  rain_add,                                                                    &
    ! Total rainfall rate at the current tile point (kg m-2 s-1).
  snow_add
    ! Total snowfall rate at the current tile point (kg m-2 s-1).

!-----------------------------------------------------------------------------

DO n = 1, nsurft

   DO j = 1, surft_pts(n)

      i = surft_index(j,n)

      rain_add = ls_rain(i) + con_rain(i)
      snow_add = ls_snow(i) + con_snow(i)

      lid_rain_water_ml(i,n) = 0.0
      
      IF (has_lid(i,n) .OR. has_vlid(i,n)) THEN

         IF (rain_add > 0.0) THEN
            lid_rain_water_ml(i,n) = rain_add * timestep / rho_water
            !lid_rain_water_ml(i,n) = lid_rain_water_ml(i,n) +                  &
             !    rain_add * timestep / rho_water
         END IF

         IF (snow_add > 0.0) THEN
            lid_snow_depth_ml(i,n) = lid_snow_depth_ml(i,n) +                 &
                 snow_add * timestep / rho_snow_const
            snow_on_lid(i,n) = .TRUE.
         END IF

      ELSE IF (exposed_water(i,n)) THEN

         lake_depth_ml(i,n) = lake_depth_ml(i,n) +                            &
              (rain_add + snow_add) * timestep / rho_water

      END IF

   END DO

END DO



END SUBROUTINE meltlake_surface_precip

END MODULE meltlake_surface_precip_mod
