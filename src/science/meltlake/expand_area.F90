! ****************************COPYRIGHT*******************************
!
! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]
!
! *****************************COPYRIGHT*******************************
! SUBROUTINE EXPAND_AREA-----------------------------------------------
!
! Update melt-lake tile fraction using the observed volume-area
! relationship:
!
!     V = a A**b
!
! V in m3, A in m2.
! aws18_lat = -66.411999
! aws18_lon = -63.435142
!
! Code Owner: s.r.shannon@reading.ac.uk
!---------------------------------------------------------------------

MODULE expand_area_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EXPAND_AREA_MOD'

CONTAINS

  SUBROUTINE expand_area(land_pts,                    &
                         nsurft,                      &
                         surft_pts,                   &
                         surft_index,                 &
                         nsnow,                       &
                         sice_ml,                     &
                         sliq_ml,                     &
                         tsnow_ml,                    &
                         exposed_water,               &
                         frac,                        &
                         lake_depth_ml,               &
                         lake_area_ml,                &
                         sice_imbal_ml,               &
                         sliq_imbal_ml,               &
                         en_imbal_ml)
    
USE model_grid_mod,     ONLY: grid_area_ij
USE model_time_mod,     ONLY: timestep_number
USE jules_meltlake_mod, ONLY: nsmax_ml
USE water_constants_mod, ONLY: hcapi, hcapw, lf, tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                       &
  land_pts,                                                                  &
  nsurft,                                                                    &
  surft_pts

INTEGER, INTENT(IN) ::                                                       &
  surft_index(land_pts),                                                     &
  ! Index of lake tile points.
  nsnow(land_pts,nsurft)
  ! Number of snow layers.
  
REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
  sice_ml(land_pts,nsurft,nsmax_ml),                                         &
    ! Snow layer ice mass (kg m-2).
  sliq_ml(land_pts,nsurft,nsmax_ml),                                         &
    ! Snow layer liquid mass (kg m-2).
  tsnow_ml(land_pts,nsurft,nsmax_ml)
    ! Snow layer temperatures (K).
  
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  frac(land_pts,nsurft),                                                     &
   ! Fractional coverage of each surface tile.
  lake_depth_ml(land_pts),                                                   &
   ! Mean lake depth over the lake tile (m).
  lake_area_ml(land_pts),                                                    &
   ! lake area in (m2)
  sice_imbal_ml(land_pts),                                                    &
   ! Cumulative compensating snow ice mass adjustment (kg m-2 gridbox).
  sliq_imbal_ml(land_pts),                                                   &
   ! Cumulative compensating snow liquid water mass adjustment (kg m-2 gridbox).
  en_imbal_ml(land_pts)
   ! Cumulative compensating snowpack enthalpy adjustment (J m-2 gridbox).

LOGICAL, INTENT(IN) ::                                                       &
  exposed_water(land_pts)
     
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------

INTEGER ::                                                                   &
  i,                                                                         &
  k, & 
  n

INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
  nonlake_ml = 10

REAL(KIND=real_jlslsm) ::                                                    &
  lake_volume,                                                               &
! Lake volume in gridbox (m3).
  new_lake_frac, &
! New lake tile fraction.
  lake_area, &
  delta_frac, &
  en_lake,                                                                   &
  en_nonlake

REAL(KIND=real_jlslsm), PARAMETER ::                                         &
  ! most obs   
  !a_va         = 0.401,                                                      &
  !b_va         = 1.018,                                                      &
  ! slow
  !a_va         = 5.0e-8,                                                     &
  !b_va         = 2.0,                                                        &
  a_va         = 1.14e-4,                                                     &
  b_va         = 1.5,                                                        &
  ! synthetic
  ! a_va         = 0.01,                                                     &
  ! b_va         = 1.20,                                                     &
  gridbox_area = 3.0814E8
  ! ERA5 0.25 degree gridbox area at AWS18 (m2).
  


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND_AREA'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,lake_volume,new_lake_frac,n,en_lake,en_nonlake)

DO k = 1, surft_pts

   i = surft_index(k)

!-----------------------------------------------------------------------------
! Only evolve frac if there is exposed water. Don't change if there is a lid
!-----------------------------------------------------------------------------
   
   IF (exposed_water(i) .AND. &
        lake_depth_ml(i) > 0.0) THEN
       
!-----------------------------------------------------------------------------
! Calculate current lake volume
!-----------------------------------------------------------------------------

      lake_volume = lake_depth_ml(i)                                         &
                  * frac(i,lake_ml)                                          &
                  * gridbox_area

!-----------------------------------------------------------------------------
! Calculate lake area from V = a A^b
!-----------------------------------------------------------------------------

      lake_area_ml(i) = (lake_volume / a_va)**(1.0 / b_va)

      new_lake_frac = lake_area_ml(i) / gridbox_area

      new_lake_frac = MIN(1.0, new_lake_frac)

!-----------------------------------------------------------------------------
! Overwrite lake depth with depth consistant with volume and area
!-----------------------------------------------------------------------------

      lake_area_ml(i) = &
         new_lake_frac * gridbox_area

      IF (lake_area_ml(i) <= 0.0) THEN
         PRINT *, 'EXPAND_AREA: zero lake area'
         PRINT *, 'Timestep = ', timestep_number
         PRINT *, 'Lake depth = ', lake_depth_ml(i)
         PRINT *, 'Lake fraction = ', frac(i,lake_ml)
         PRINT *, 'Lake volume = ', lake_volume
         PRINT *, 'Lake area = ', lake_area_ml(i)
         PRINT *, 'Exposed water = ', exposed_water(i)
         STOP 'EXPAND_AREA: zero lake area'
      END IF

      lake_depth_ml(i) = &
         lake_volume / lake_area_ml(i)

!-----------------------------------------------------------------------------
! The compensating adjustment to masses from changing frac
!-----------------------------------------------------------------------------

      delta_frac = new_lake_frac - frac(i,lake_ml)

      IF (delta_frac /= 0.0) THEN

         sice_imbal_ml(i) = sice_imbal_ml(i) + delta_frac * ( &
              SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) &
              - SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml))) )

         sliq_imbal_ml(i) = sliq_imbal_ml(i) + delta_frac * ( &
              SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) &
              - SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml))) )

!-----------------------------------------------------------------------------
! Calculate snowpack enthalpy for each tile (J m-2 tile).
! 
! 
! H_snow = SUM_j [ (sice_j*hcapi + sliq_j*hcapw) * (tsnow_j - tm)
!                  + lf*sliq_j ]
!
! Sensible energy: energy associated with snow temperature.
! Latent energy:   energy associated with liquid water rather than ice.
!
! Reference state: ice at 0 degrees C has zero enthalpy.
!-----------------------------------------------------------------------------
! ** Cold ice has negative enthalpy.
! ** Ice at 0°C has zero enthalpy.
! ** Liquid water at 0°C has positive enthalpy.
!-----------------------------------------------------------------------------

      en_lake = 0.0

      DO n = 1, nsnow(i,lake_ml)

         en_lake = en_lake + &
              (sice_ml(i,lake_ml,n) * hcapi + &
              sliq_ml(i,lake_ml,n) * hcapw) * &
              (tsnow_ml(i,lake_ml,n) - tm) !+ &
              !lf * sliq_ml(i,lake_ml,n)

      END DO

      en_nonlake = 0.0

      DO n = 1, nsnow(i,nonlake_ml)
         
         en_nonlake = en_nonlake + &
              (sice_ml(i,nonlake_ml,n) * hcapi + &
              sliq_ml(i,nonlake_ml,n) * hcapw) * &
              (tsnow_ml(i,nonlake_ml,n) - tm) !+ &
              !lf * sliq_ml(i,nonlake_ml,n)

      END DO

!-----------------------------------------------------------------------------
! Accumulate compensating gridbox enthalpy adjustment (J m-2 gridbox).
!-----------------------------------------------------------------------------

      en_imbal_ml(i) = en_imbal_ml(i) + &
           delta_frac * (en_nonlake - en_lake)
         
      END IF

      
      frac(i,lake_ml)    = new_lake_frac
      frac(i,nonlake_ml) = 1.0 - new_lake_frac

      
      !WRITE(*,*) 'EXPAND_AREA: timestep = ', timestep_number

      !WRITE(*,*) 'Old lake frac = ', frac(i,lake_ml)
      !WRITE(*,*) 'New lake frac = ', new_lake_frac
      !WRITE(*,*) 'Delta frac    = ', delta_frac

      !WRITE(*,*) 'Lake volume (m3) = ', lake_volume
      !WRITE(*,*) 'Lake area (m2)   = ', lake_area_ml(i)
      !WRITE(*,*) 'Lake depth (m)   = ', lake_depth_ml(i)

      !WRITE(*,*) 'Lake snow ice mass (kg m-2) = ', &
      !     SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml)))

      !WRITE(*,*) 'Non-lake snow ice mass (kg m-2) = ', &
      !     SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml)))

      !WRITE(*,*) 'sice_imbal_ml (kg m-2 gridbox) = ', sice_imbal_ml(i)
      !WRITE(*,*) 'sliq_imbal_ml (kg m-2 gridbox) = ', sliq_imbal_ml(i)
      !WRITE(*,*) 'en_imbal_ml (J m-2 gridbox) = ', en_imbal_ml(i)

      !WRITE(*,*) 'Lake snow enthalpy = ', en_lake
      !WRITE(*,*) 'Non-lake snow enthalpy = ', en_nonlake

      !WRITE(*,*) 'Enthalpy adjustment this timestep = ', &
      !     delta_frac * (en_nonlake - en_lake)

      !WRITE(*,*) 'Cumulative latent adjustment = ', &
      !     lf * sliq_imbal_ml(i)

      !WRITE(*,*) 'Cumulative sensible adjustment = ', &
      !     en_imbal_ml(i) - lf * sliq_imbal_ml(i)

   END IF

END DO


!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE expand_area

END MODULE expand_area_mod

