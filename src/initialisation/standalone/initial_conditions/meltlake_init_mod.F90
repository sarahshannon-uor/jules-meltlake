#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE meltlake_init_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE meltlake_init(ainfo,progs)

USE jules_meltlake_mod,  ONLY: l_meltlake, nsmax_ml, dzsnow_ml,        &
                               rho_firn_efold, firn_depth_max

USE ancil_info,          ONLY: land_pts, nsurft, lice_pts
USE jules_surface_mod,   ONLY: l_elev_land_ice
USE jules_snow_mod,      ONLY: rho_snow_const, rho_snow_fresh, r0
USE water_constants_mod, ONLY: rho_ice, tm

USE layersnow_mod,   ONLY: layersnow

USE um_types, ONLY: real_jlslsm

!JULES TYPEs
USE ancil_info,    ONLY: ainfo_type
USE prognostics,   ONLY: progs_type

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the initial conditions for snowpack when melt lake model is on 
!   Note: the other melt lake variables are initialised in 
!         src/control/shared/meltlake_vars_mod.F90 where they are allocated 
!         similar to flake variables
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Arguments

!JULES TYPEs
TYPE(ainfo_type), INTENT(IN OUT) :: ainfo
TYPE(progs_type), INTENT(IN OUT) :: progs

! Work variables
INTEGER :: i, j, k, l, n, m  ! Index variables

!-----------------------------------------------------------------------------
! For idealised tests warm top --> cold bottom
!REAL(KIND=real_jlslsm), PARAMETER :: T_bot = 253.15, T_top = 272.15
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Larsen C Buzzard et al 2018, cold top (-20oC)--> warm bottom (-10oC)
! Put these in jules_meltlake.nml
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER ::                               &
 T_bot = 263.15, &
  ! Top temperature (K) 
 T_top = 253.15, &
  ! Bottom temperature (K) 
 rho_sfc = 500.0
  ! Surface density (kg/m3)
         

REAL(KIND=real_jlslsm) ::                                          &
 z_top, &
  ! Cumulative depth from the snow surface to the top of the
  ! current layer (m)
 z_mid
  ! Depth from the snow surface to the midpoint of the
  ! current layer (m)

progs%sice_surft_ml(:,:,:)     = 0.0
progs%sliq_surft_ml(:,:,:)     = 0.0
progs%rgrainl_surft_ml(:,:,:)  = r0
progs%tsnow_surft_ml(:,:,:)    = tm
progs%ds_surft_ml(:,:,:)       = 0.0
progs%rho_snow_surft_ml(:,:,:) = 0.0
progs%snow_surft_ml(:,:)       = 0.0
progs%snowdepth_surft(:,:)     = 0.0
progs%nsnow_surft(:,:)         = 0


IF (l_meltlake .AND. l_elev_land_ice) THEN

!-------------------------------------------------------------------
! Set initial snowpack depth, value from jules_meltlake.nml
!-------------------------------------------------------------------

   DO n = 1, nsurft

      IF (ainfo%l_lice_surft(n)) THEN

         DO j = 1, lice_pts
            i = ainfo%lice_index(j)

            progs%snowdepth_surft(i,n) = firn_depth_max

         END DO

      END IF

   END DO

  !-------------------------------------------------------------------
  ! Use layersnow to get number and thickness of snow layers 
  !-------------------------------------------------------------------

   DO n = 1, nsurft

      IF (ainfo%l_lice_surft(n)) THEN

         CALL layersnow(                                               &
              land_pts,                                                &
              ainfo%surft_pts(n),                                      &
              ainfo%surft_index(:,n),                                  &
              nsmax_ml,                                                &
              dzsnow_ml,                                               &
              progs%snowdepth_surft(:,n),                              &
              progs%nsnow_surft(:,n),                                  &
              progs%ds_surft_ml(:,n,:))
         
      END IF

   END DO

  !-------------------------------------------------------------------
  ! Initialise snow temperature, density and mass
  !-------------------------------------------------------------------

   DO n = 1, nsurft

      IF (ainfo%l_lice_surft(n)) THEN

         DO j = 1, lice_pts
            i = ainfo%lice_index(j)

            z_top = 0.0
            progs%snow_surft_ml(i,n) = 0.0

            DO k = 1, progs%nsnow_surft(i,n)
               
  !-------------------------------------------------------------------
  ! Midpoint depth of layer
  !-------------------------------------------------------------------             
               z_mid = z_top + 0.5 * progs%ds_surft_ml(i,n,k)

  !-------------------------------------------------------------------             
  ! Linear temperature profile, Buzzard et al 2018
  !-------------------------------------------------------------------             
               progs%tsnow_surft_ml(i,n,k) =                            &
                    T_top + (T_bot - T_top) *                           &
                    z_mid / firn_depth_max

   !-------------------------------------------------------------------
   ! Exponential firn density profile, Buzzard et al 2018
   !------------------------------------------------------------------- 
               progs%rho_snow_surft_ml(i,n,k) =                         &
                    rho_ice -                                           &
                    (rho_ice - rho_sfc) *                               &
                    EXP(-(1.9 / rho_firn_efold) * z_mid)

   !-------------------------------------------------------------------
   ! Grain size
   !-------------------------------------------------------------------
               progs%rgrainl_surft_ml(i,n,k) = r0

   
   !-------------------------------------------------------------------
   ! Ice mass is a function of density. Assume no liquid in snow 
   !-------------------------------------------------------------------
               progs%sice_surft_ml(i,n,k) =                             &
                    progs%rho_snow_surft_ml(i,n,k) *                    &
                    progs%ds_surft_ml(i,n,k)

               progs%sliq_surft_ml(i,n,k) = 0.0

   !-------------------------------------------------------------------
   ! Total column snow mass is sum of layer masses 
   !-------------------------------------------------------------------
               progs%snow_surft_ml(i,n) =                               &
                    progs%snow_surft_ml(i,n) +                          &
                    progs%sice_surft_ml(i,n,k)

   !-------------------------------------------------------------------
   ! Moving down to the top of the next layer
   !-------------------------------------------------------------------
               z_top = z_top + progs%ds_surft_ml(i,n,k)

            END DO

         END DO

      END IF

   END DO

END IF

   !-------------------------------------------------------------------
   ! write check 
   !-------------------------------------------------------------------
i = 1
n = 9

WRITE(*,*) 'MELTLAKE INITIAL CONDITIONS'
WRITE(*,'(A,I4,A,I4)') 'land point = ', i, '   surface tile = ', n
WRITE(*,'(A,I4,A,F10.4,A,F12.4)')                           &
     'nsnow = ', progs%nsnow_surft(i,n),                    &
     '   snow depth = ', progs%snowdepth_surft(i,n),         &
     '   snow mass = ', progs%snow_surft_ml(i,n)

WRITE(*,'(A)')                                               &
     ' layer    z_top      z_mid       depth       temp_C'   &
     //'      density      ice_mass    liquid_mass   grain_size'

z_top = 0.0

DO k = 1, progs%nsnow_surft(i,n)

   z_mid = z_top + 0.5 * progs%ds_surft_ml(i,n,k)

   WRITE(*,'(I6,8F12.4)')                                    &
        k,                                                    &
        z_top,                                                &
        z_mid,                                                &
        progs%ds_surft_ml(i,n,k),                             &
        progs%tsnow_surft_ml(i,n,k) - 273.15,                &
        progs%rho_snow_surft_ml(i,n,k),                       &
        progs%sice_surft_ml(i,n,k),                           &
        progs%sliq_surft_ml(i,n,k),                           &
        progs%rgrainl_surft_ml(i,n,k)

   z_top = z_top + progs%ds_surft_ml(i,n,k)

END DO


RETURN

END SUBROUTINE meltlake_init

END MODULE meltlake_init_mod
#endif

