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
USE jules_snow_mod,      ONLY: rho_snow_const, rho_snow_fresh
USE water_constants_mod, ONLY: rho_ice

!USE jules_radiation_mod, ONLY: l_snow_albedo, l_embedded_snow
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

REAL(KIND=real_jlslsm), PARAMETER :: T_min = 253.15, T_max = 272.15!263.15

REAL(KIND=real_jlslsm), PARAMETER :: rho_sfc = 500.0

REAL(KIND=real_jlslsm) :: snow_mass

REAL(KIND=real_jlslsm) :: dz_cum

REAL(KIND=real_jlslsm) :: dzsnow_cumulative(land_pts,nsurft,nsmax_ml)
              ! cumulative snow depth per layer


progs%sice_surft_ml(:,:,:)    = 0.0
progs%sliq_surft_ml(:,:,:)    = 0.0
progs%tsnow_surft_ml(:,:,:)   = 273.15
progs%ds_surft_ml(:,:,:)      = 0.0 
progs%rho_snow_surft_ml(:,:,:)= 0.0
progs%snow_surft_ml(:,:)      = 0.0
progs%snowdepth_surft(:,:)    = 0.0
  
IF (l_elev_land_ice) THEN

! get cumulative snowdepth layers from jules_meltlake.nml
dz_cum = 0.0
DO k = 1, nsmax_ml
   dz_cum = dz_cum + dzsnow_ml(k)
   dzsnow_cumulative(:,:,k) = dz_cum
END DO


DO n = 1,nsurft
    DO j = 1,lice_pts
       i = ainfo%lice_index(j)
  
       ! overwrite the snowdepth from total_snow_init_mod.F90
            IF ( ainfo%l_lice_surft(n)) THEN
 
            progs%snowdepth_surft(i,n) =  firn_depth_max
   
                DO k = 1,nsmax_ml

!-------------------------------------------------------------------------------
! For equally spaced levels
!-------------------------------------------------------------------------------
                    progs%ds_surft_ml(i,n,k) = firn_depth_max / REAL(nsmax_ml) 
                    
                    ! get cumulative snowdepth layers assuming equally spaced layers 
                    dzsnow_cumulative(i,n,k) = (REAL(k - 1) / REAL(nsmax_ml - 1)) * firn_depth_max 
!-------------------------------------------------------------------------------
! Initialise snowpack temperature. Default is a warm top, cold bottom
!-------------------------------------------------------------------------------

                    progs%tsnow_surft_ml(i,n,k) = T_max - ( (k-1) * (T_max - T_min) ) / REAL(nsmax_ml-1)
! for idealised test to control where refreezing of meltwater happens in the snowpack
            !IF (k <= 20) THEN
                ! 0–4 m: warm near-melting snow
             !   progs%tsnow_surft_ml(i,n,k) = 272.65   ! K

            !ELSE IF (k <= 25) THEN
               ! 4–5 m: cold trap, linear ramp
             !  progs%tsnow_surft_ml(i,n,k) = 271.15 + &
             !  REAL(k-21) / REAL(25-21) * (253.15 - 271.15)

            !ELSE
             ! below 3 m: cold background
            !   progs%tsnow_surft_ml(i,n,k) = 253.15   ! K

            !END IF                     
!-------------------------------------------------------------------------------
! Increase density with depth using e-folding value
!-------------------------------------------------------------------------------

               !progs%rho_snow_surft_ml(i,n,k) = 700.0
              
              ! monarchs init desnity profile   
                    progs%rho_snow_surft_ml(i,n,k) = rho_ice - &
                         (rho_ice - rho_sfc) * EXP( - (1.9 / rho_firn_efold) * dzsnow_cumulative(i,n,k) )
               print *, 'k, rho, temp', k, progs%rho_snow_surft_ml(i,n,k),progs%tsnow_surft_ml(i,n,k)-273.15    
           END DO
        END IF
    END DO
END DO

!-------------------------------------------------------------------------------
! Calculate snow layer thicknesses - testing varaible ds instead of fixed 
!-------------------------------------------------------------------------------
  !DO n = 1,nsurft
  !  CALL layersnow(land_pts, ainfo%surft_pts(n), ainfo%surft_index(:,n),       &
  !                 nsmax_ml, dzsnow_ml, progs%snowdepth_surft_ml(:,n),         &
!                   progs%nsnow_surft(:,n), progs%ds_surft_ml(:,n,:))
!  END DO

!print *, 'dzsnow_cumulative',dzsnow_cumulative(:,9,:)

!print *, 'progs%rho_snow_surft_ml(:,n)',progs%rho_snow_surft_ml(:,9,:)

!print *, 'progs%snowdepth_surft(:,n)',progs%snowdepth_surft(:,9)

!print *, 'progs%ds_surft_ml',progs%ds_surft_ml(:,9,:)

!print *, 'progs%nsnow_surft',progs%nsnow_surft(:,9)


!-------------------------------------------------------------------------------
! get snowmass on tile (kgm-2). Assume there is no liquid content 
!-------------------------------------------------------------------------------
 snow_mass = 0.0
   
DO n = 1,nsurft
   DO j = 1,lice_pts
      i = ainfo%lice_index(j)
 
      progs%snow_surft_ml(i,n) = 0.0

     ! snow layer depths are equally spaced for now
      !dz = progs%ds_surft_ml(i,n,2) - progs%ds_surft_ml(i,n,1) 
      
      DO k = 1,nsmax_ml
       !progs%sice_surft_ml(i,n,k) = progs%rho_snow_surft_ml(i,n,k) * progs%ds_surft_ml(i,n,k)
        progs%sice_surft_ml(i,n,k) = progs%rho_snow_surft_ml(i,n,k) * dzsnow_ml(k)
        progs%sliq_surft_ml(i,n,k) = 0.0
        progs%snow_surft_ml(i,n) = progs%snow_surft_ml(i,n) + progs%sice_surft_ml(i,n,k)
      END DO

   END DO
 END DO 

  



!print *, 'rho_firn_efold', rho_firn_efold
!print *, 'firn_depth_max', firn_depth_max
!print *, 'progs%ds_surft_ml', progs%ds_surft_ml(:,9,:)
!print *, 'dzsnow_ml', dzsnow_ml
!print *, 'progs%tsnow_surft_ml', progs%tsnow_surft_ml(:,9,:)-273.15
!print *, 'progs%rho_snow_surft_ml', progs%rho_snow_surft_ml(:,9,:)
!print *, 'progs%snow_surft_ml(i)',progs%snow_surft_ml
!print *, 'progs%snowdepth_surft_ml(i)',progs%snowdepth_surft_ml
!print*, 'progs%ice_mass_snow_ml(i,n,k)',progs%sice_surft_ml



END IF



RETURN

END SUBROUTINE meltlake_init
END MODULE meltlake_init_mod
#endif
