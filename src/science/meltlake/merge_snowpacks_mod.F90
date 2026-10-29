
! ****************************COPYRIGHT*******************************
!
! (c) [University of Reading]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC237]
!
! *****************************COPYRIGHT*******************************
! SUBROUTINE MERGE_SNOWPACKS------------------------------------
!
! Description:
!   Merge lake and non-lake snowpacks after complete lake freeze-over.
!
! Method:
!   Area-weight layer thickness, ice mass, liquid mass and sensible
!   energy. Apply cumulative adjustments from lake area expansion
!   and check conservation.
!
! Code Owner: s.r.shannon@reading.ac.uk
!
!---------------------------------------------------------------------

MODULE merge_snowpacks_mod

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                      &
  ModuleName='MERGE_SNOWPACKS_MOD'

CONTAINS

SUBROUTINE merge_snowpacks(land_pts,          &
                           nsurft,            &
                           surft_pts,         &
                           surft_index,       &
                           did_insert_lid,    &
                           frac,              &
                           nsnow,             &
                           sice_ml,           &
                           sliq_ml,           &
                           tsnow_ml,          &
                           ds_ml,             &
                           snowdepth,         &
                           snow_surft,        &
                          ! ice_lens_index,    &
                           sice_imbal_ml,     &
                           sliq_imbal_ml,     &
                           en_imbal_ml)

USE jules_meltlake_mod, ONLY: nsmax_ml
USE model_time_mod,     ONLY: timestep_number
USE water_constants_mod, ONLY: hcapi, hcapw, tm

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

LOGICAL, INTENT(IN) ::                                                       &
  did_insert_lid(land_pts)
    ! True where the lake lid has been inserted into the snowpack.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                        &
  frac(land_pts,nsurft),                                                     &
    ! Fractional coverage of each surface tile.
  !ice_lens_index(land_pts,nsurft),                                           &
    ! Snow layer containing the impermeable ice lens.
  sice_imbal_ml(land_pts),                                                   &
    ! Cumulative compensating ice mass adjustment (kg m-2 gridbox).
  sliq_imbal_ml(land_pts),                                                   &
    ! Cumulative compensating liquid mass adjustment (kg m-2 gridbox).
  en_imbal_ml(land_pts)
    ! Cumulative compensating sensible energy adjustment (J m-2 gridbox).

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
    ! Total snowpack depth (m).
  snow_surft(land_pts,nsurft)
    ! Snow mass on tiles (kg m-2).

INTEGER, PARAMETER :: merge_unit = 99

CHARACTER(LEN=*), PARAMETER :: &
  merge_file = 'snowpack_merge_t6134.csv'

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
  nmerge
    ! Number of layers in the preliminary merged snowpack.

REAL(KIND=real_jlslsm) ::                                                    &
  sice_merge(nsmax_ml),                                                      &
  sliq_merge(nsmax_ml),                                                      &
  ds_merge(nsmax_ml),                                                        &
  en_merge(nsmax_ml)

INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
    ! Melt-lake tile index.
  nonlake_ml = 10
    ! Non-lake elevated-ice tile index.

REAL(KIND=real_jlslsm) ::                                                    &
  sice_gridbox_total,                                                        &
    ! Ice mass conservation target (kg m-2 gridbox).
  sliq_gridbox_total,                                                        &
    ! Liquid mass conservation target (kg m-2 gridbox).
  en_gridbox_total
    ! Sensible energy conservation target (J m-2 gridbox).

REAL(KIND=real_jlslsm) ::                                                    &
  tsnow_merge(nsmax_ml),                                                     &
  en_lake, en_nonlake,                                                       &
  delta_t,                                                                   &
  sice_error, sliq_error, en_error

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MERGE_SNOWPACKS'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,ns,nmerge,sice_merge,sliq_merge,ds_merge,en_merge, &
!$OMP         tsnow_merge,en_lake,en_nonlake, &
!$OMP         sice_gridbox_total,sliq_gridbox_total,en_gridbox_total, &
!$OMP         delta_t,sice_error,sliq_error,en_error)

DO k = 1, surft_pts

   i = surft_index(k)

   IF (did_insert_lid(i)) THEN

      print *, 'did_insert_lid at timestep_number', timestep_number

!-----------------------------------------------------------------------------
! Write original snowpacks after lid insertion, before merging.
!-----------------------------------------------------------------------------

      !IF (timestep_number == 6134 .AND. k == 1) THEN

        ! OPEN(UNIT=merge_unit,FILE=merge_file, &
        !      STATUS='REPLACE',ACTION='WRITE')

        ! WRITE(merge_unit,'(A)') &
        !      'column,layer,ds_m,sice_kgm2,sliq_kgm2,T_C,rho_kgm3,snowdepth_m'

       !  DO ns = 1, nsnow(i,lake_ml)

         !   WRITE(merge_unit,'(A,",",I0,6(",",ES17.9))') &
        !         'lake',ns, &
        !         ds_ml(i,lake_ml,ns), &
        !         sice_ml(i,lake_ml,ns), &
        !         sliq_ml(i,lake_ml,ns), &
        !         tsnow_ml(i,lake_ml,ns)-tm, &
        !         (sice_ml(i,lake_ml,ns)+sliq_ml(i,lake_ml,ns)) / &
        !         ds_ml(i,lake_ml,ns), &
        !         snowdepth(i,lake_ml)

        ! END DO

        ! DO ns = 1, nsnow(i,nonlake_ml)

         !   WRITE(merge_unit,'(A,",",I0,6(",",ES17.9))') &
         !        'nonlake',ns, &
         !        ds_ml(i,nonlake_ml,ns), &
         !        sice_ml(i,nonlake_ml,ns), &
         !        sliq_ml(i,nonlake_ml,ns), &
         !        tsnow_ml(i,nonlake_ml,ns)-tm, &
         !        (sice_ml(i,nonlake_ml,ns)+sliq_ml(i,nonlake_ml,ns)) / &
         !        ds_ml(i,nonlake_ml,ns), &
         !        snowdepth(i,nonlake_ml)

         !END DO

      !END IF

!-----------------------------------------------------------------------------
! The lid has been inserted into the snopack changing sice , sliq, tsnow
! so recalcualte enthalpy from sensible energy
!-----------------------------------------------------------------------------

      en_lake = 0.0

      DO ns = 1, nsnow(i,lake_ml)

         en_lake = en_lake + &
              (sice_ml(i,lake_ml,ns)*hcapi + &
               sliq_ml(i,lake_ml,ns)*hcapw) * &
              (tsnow_ml(i,lake_ml,ns)-tm)

      END DO

      en_nonlake = 0.0

      DO ns = 1, nsnow(i,nonlake_ml)

         en_nonlake = en_nonlake + &
              (sice_ml(i,nonlake_ml,ns)*hcapi + &
               sliq_ml(i,nonlake_ml,ns)*hcapw) * &
              (tsnow_ml(i,nonlake_ml,ns)-tm)

      END DO

!-----------------------------------------------------------------------------
! Calculate conservation targets for the merged snowpack.
! These targets are used to verify conservation after merging.
!-----------------------------------------------------------------------------

      sice_gridbox_total = frac(i,lake_ml) * &
                       SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml))) + &
                       frac(i,nonlake_ml) * &
                       SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) + &
                       sice_imbal_ml(i)

      sliq_gridbox_total = frac(i,lake_ml) * &
                       SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml))) + &
                       frac(i,nonlake_ml) * &
                       SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) + &
                       sliq_imbal_ml(i)

      en_gridbox_total = frac(i,lake_ml)*en_lake + &
                     frac(i,nonlake_ml)*en_nonlake + &
                     en_imbal_ml(i)

!-----------------------------------------------------------------------------
! Merge snowpack layers for the lake and nonlake tiles
! Each snowpack maney have different number of nsnow
!
!-----------------------------------------------------------------------------

      nmerge = MAX(nsnow(i,lake_ml),nsnow(i,nonlake_ml))

      sice_merge = 0.0
      sliq_merge = 0.0
      ds_merge = 0.0
      en_merge = 0.0

      DO ns = 1, nmerge

         sice_merge(ns) = frac(i,lake_ml)*sice_ml(i,lake_ml,ns) + &
                          frac(i,nonlake_ml)*sice_ml(i,nonlake_ml,ns)

         sliq_merge(ns) = frac(i,lake_ml)*sliq_ml(i,lake_ml,ns) + &
                          frac(i,nonlake_ml)*sliq_ml(i,nonlake_ml,ns)

         en_merge(ns) = frac(i,lake_ml) * &
              (sice_ml(i,lake_ml,ns)*hcapi + &
               sliq_ml(i,lake_ml,ns)*hcapw) * &
              (tsnow_ml(i,lake_ml,ns)-tm) + &
              frac(i,nonlake_ml) * &
              (sice_ml(i,nonlake_ml,ns)*hcapi + &
               sliq_ml(i,nonlake_ml,ns)*hcapw) * &
              (tsnow_ml(i,nonlake_ml,ns)-tm)

         IF (ns <= nsnow(i,lake_ml)) THEN
            ds_merge(ns) = ds_merge(ns) + &
                           frac(i,lake_ml)*ds_ml(i,lake_ml,ns)
         END IF

         IF (ns <= nsnow(i,nonlake_ml)) THEN
            ds_merge(ns) = ds_merge(ns) + &
                           frac(i,nonlake_ml)*ds_ml(i,nonlake_ml,ns)
         END IF

      END DO

!-----------------------------------------------------------------------------
! Distribute accumulated mass and sensible energy imbalances across
! the merged snowpack.
!-----------------------------------------------------------------------------

! Calculate merged layer temperature from merged enthalpy
      tsnow_merge(1:nmerge) = tm + en_merge(1:nmerge) / &
           (sice_merge(1:nmerge)*hcapi + sliq_merge(1:nmerge)*hcapw)

! Distribute ice imbalance proportionally
      sice_merge(1:nmerge) = sice_merge(1:nmerge) * &
           (sice_gridbox_total/SUM(sice_merge(1:nmerge)))

! Distribute liquid imbalance proportionally
      IF (SUM(sliq_merge(1:nmerge)) > 0.0) THEN
         sliq_merge(1:nmerge) = sliq_merge(1:nmerge) * &
              (sliq_gridbox_total/SUM(sliq_merge(1:nmerge)))
      END IF

! Apply uniform temperature adjustment to conserve sensible energy
      delta_t = (en_gridbox_total - &
           SUM((sice_merge(1:nmerge)*hcapi + sliq_merge(1:nmerge)*hcapw) * &
           (tsnow_merge(1:nmerge)-tm))) / &
           SUM(sice_merge(1:nmerge)*hcapi + sliq_merge(1:nmerge)*hcapw)

      tsnow_merge(1:nmerge) = tsnow_merge(1:nmerge) + delta_t

!-----------------------------------------------------------------------------
! Check conservation of ice mass, liquid mass and sensible energy.
!-----------------------------------------------------------------------------

      sice_error = SUM(sice_merge(1:nmerge)) - sice_gridbox_total

      sliq_error = SUM(sliq_merge(1:nmerge)) - sliq_gridbox_total

      en_error = SUM((sice_merge(1:nmerge)*hcapi + &
                      sliq_merge(1:nmerge)*hcapw) * &
                     (tsnow_merge(1:nmerge)-tm)) - en_gridbox_total

      WRITE(*,*) 'Ice mass error (kg m-2) = ', sice_error
      WRITE(*,*) 'Liquid mass error (kg m-2) = ', sliq_error
      WRITE(*,*) 'Sensible energy error (J m-2) = ', en_error

!-----------------------------------------------------------------------------
! Write merged snowpack after applying imbalances.
!-----------------------------------------------------------------------------

     ! IF (timestep_number == 6134 .AND. k == 1) THEN

      !   DO ns = 1, nmerge

       !     WRITE(merge_unit,'(A,",",I0,6(",",ES17.9))') &
       !          'merged',ns, &
       !          ds_merge(ns), &
       !          sice_merge(ns), &
       !          sliq_merge(ns), &
       !          tsnow_merge(ns)-tm, &
       !          (sice_merge(ns)+sliq_merge(ns))/ds_merge(ns), &
       !          SUM(ds_merge(1:nmerge))

       !  END DO

       !  CLOSE(merge_unit)

      !END IF

      
      IF (ABS(sice_error) > MAX(0.001,1.0E-6*ABS(sice_gridbox_total)) .OR. &
           ABS(sliq_error) > MAX(0.001,1.0E-6*ABS(sliq_gridbox_total)) .OR. &
           ABS(en_error) > MAX(1000.0,1.0E-6*ABS(en_gridbox_total))) THEN

         STOP 'MERGE_SNOWPACKS: Mass or energy not conserved'

      END IF




!-----------------------------------------------------------------------------
! Assign corrected merged snowpack to the nonlake tile
!-----------------------------------------------------------------------------

      nsnow(i,nonlake_ml) = nmerge

      sice_ml(i,nonlake_ml,1:nmerge) = sice_merge(1:nmerge)

      sliq_ml(i,nonlake_ml,1:nmerge) = sliq_merge(1:nmerge)

      tsnow_ml(i,nonlake_ml,1:nmerge) = tsnow_merge(1:nmerge)

      ds_ml(i,nonlake_ml,1:nmerge) = ds_merge(1:nmerge)

      snowdepth(i,nonlake_ml) = SUM(ds_merge(1:nmerge))

      snow_surft(i,nonlake_ml) = SUM(sice_merge(1:nmerge)) + &
                                SUM(sliq_merge(1:nmerge))


      
   END IF ! did_insert_lid

END DO

!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE merge_snowpacks

END MODULE merge_snowpacks_mod
