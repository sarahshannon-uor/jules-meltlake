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
                         ds_ml,                       &
                         snowdepth,                   &
                         snow_surft,                  &
                         rho_snow_ml,                 &
                         rho_snow_grnd,               &
                         exposed_water,               &
                         frac,                        &
                         lake_depth_ml,               &
                         lake_area_ml,                &
                         sice_imbal_ml,               &
                         sliq_imbal_ml,               &
                         en_imbal_ml,                 &
                         correct_imbalance,           &
                         correct_frac_immediately,    &
                         correct_lake_tile,           &
                         immediate_merge)
         
USE model_time_mod,     ONLY: timestep_number
USE water_constants_mod, ONLY: hcapi, hcapw, tm
USE jules_meltlake_mod, ONLY: nsmax_ml

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
  surft_index(land_pts)
  ! Index of lake tile points.

INTEGER, INTENT(IN OUT) :: nsnow(land_pts,nsurft)
  ! Number of snow layers.
  
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                    &
  sice_ml(land_pts,nsurft,nsmax_ml),                                         &
  sliq_ml(land_pts,nsurft,nsmax_ml),                                         &
  tsnow_ml(land_pts,nsurft,nsmax_ml),                                        &
  ds_ml(land_pts,nsurft,nsmax_ml),                                           &
  snowdepth(land_pts,nsurft),                                                &
  snow_surft(land_pts,nsurft),                                               &
  rho_snow_ml(land_pts,nsurft,nsmax_ml),                                      &
  rho_snow_grnd(land_pts,nsurft),                                            &
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

LOGICAL, INTENT(IN) :: correct_frac_immediately
LOGICAL, INTENT(IN) :: correct_imbalance
LOGICAL, INTENT(IN) :: correct_lake_tile
LOGICAL, INTENT(IN) :: immediate_merge

LOGICAL, INTENT(IN) ::                                                       &
  exposed_water(land_pts)
     
!-----------------------------------------------------------------------------
! Local scalars
!-----------------------------------------------------------------------------

INTEGER ::                                                                   &
  i,                                                                         &
  k,                                                                         &
  n

INTEGER, PARAMETER ::                                                        &
  lake_ml    = 9,                                                            &
  nonlake_ml = 10

REAL(KIND=real_jlslsm) ::                                                    &
  lake_volume,                                                               &
! Lake volume in gridbox (m3).
  new_lake_frac, &
! New lake tile fraction.
  delta_frac,                                                                 &
  en_lake,                                                                   &
  en_nonlake


!Larger va_b: slower relative increase in lake area with volume.
!Smaller va_b: faster relative increase in lake area with volume.
!Larger va_a: Deeper lakes for a given area, or smaller areas for a given volume.
!Smaller va_a: Shallower lakes for a given area, or larger areas for a given volume.
! Reference curve: 10 m depth at approximately 20% coverage

REAL(KIND=real_jlslsm), PARAMETER ::                                         &
 a_va = 0.045,                                                           &
 b_va = 1.2,                                                                 &
 gridbox_area = 3.0814E8
  ! ERA5 0.25 degree gridbox area at AWS18 (m2).
  
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXPAND_AREA'

!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(k,i,n,lake_volume,new_lake_frac,delta_frac,en_lake,en_nonlake)

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

      !lake_volume = lake_depth_ml(i)                                         &
      !            * frac(i,lake_ml)                                          &
      !            * gridbox_area

!-----------------------------------------------------------------------------
! Calculate lake area from V = a A^b
!-----------------------------------------------------------------------------

      !lake_area_ml(i) = (lake_volume / a_va)**(1.0 / b_va)

      IF (frac(i,lake_ml) > 0.0) THEN

         ! Existing lake
         lake_volume = lake_depth_ml(i) * frac(i,lake_ml) * gridbox_area

         lake_area_ml(i) = (lake_volume / a_va)**(1.0 / b_va)

      ELSE

         ! New lake after complete freeze-over
         ! From lake_depth = a_va * lake_area**(b_va - 1)
         print *, '---------- growing from nothing ----------- ' 
         lake_area_ml(i) = &
              (lake_depth_ml(i) / a_va)**(1.0 / (b_va - 1.0))

         lake_volume = lake_depth_ml(i) * lake_area_ml(i)

      END IF

      

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

         IF (correct_imbalance .AND. correct_frac_immediately) THEN
            
            
            CALL correct_snowpack_frac(i, land_pts, nsurft, nsnow,          &
                 sice_ml, sliq_ml, tsnow_ml, ds_ml, snowdepth, snow_surft, &
                 rho_snow_ml, rho_snow_grnd, frac(i,lake_ml), new_lake_frac, &
                 correct_lake_tile, immediate_merge)

         ELSE

            ! Accumulate compensating gridbox snow ice and liquid mass
            ! adjustments for application when snowpacks merge.
            sice_imbal_ml(i) = sice_imbal_ml(i) + delta_frac * ( &
                 SUM(sice_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) &
                 - SUM(sice_ml(i,lake_ml,1:nsnow(i,lake_ml))) )

            sliq_imbal_ml(i) = sliq_imbal_ml(i) + delta_frac * ( &
                 SUM(sliq_ml(i,nonlake_ml,1:nsnow(i,nonlake_ml))) &
                 - SUM(sliq_ml(i,lake_ml,1:nsnow(i,lake_ml))) )

            ! Sensible snowpack enthalpy relative to ice at tm.
            ! Keep the convention used by merge_snowpacks.
            en_lake = 0.0
            DO n = 1, nsnow(i,lake_ml)
               en_lake = en_lake + &
                    (sice_ml(i,lake_ml,n)*hcapi + &
                     sliq_ml(i,lake_ml,n)*hcapw) * &
                    (tsnow_ml(i,lake_ml,n) - tm)
            END DO

            en_nonlake = 0.0
            DO n = 1, nsnow(i,nonlake_ml)
               en_nonlake = en_nonlake + &
                    (sice_ml(i,nonlake_ml,n)*hcapi + &
                     sliq_ml(i,nonlake_ml,n)*hcapw) * &
                    (tsnow_ml(i,nonlake_ml,n) - tm)
            END DO

            en_imbal_ml(i) = en_imbal_ml(i) + &
                 delta_frac * (en_nonlake - en_lake)

         END IF

      END IF

      frac(i,lake_ml)    = new_lake_frac
      frac(i,nonlake_ml) = 1.0 - new_lake_frac

     
   END IF

END DO


!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                         &
                        zhook_out,zhook_handle)

RETURN

END SUBROUTINE expand_area

!-----------------------------------------------------------------------------
! Conserve gridbox snow ice, liquid water and sensible energy when frac changes.
! Apply the correction to the selected target tile.
!-----------------------------------------------------------------------------
SUBROUTINE correct_snowpack_frac(i, land_pts, nsurft, nsnow,             &
     sice_ml, sliq_ml, tsnow_ml, ds_ml, snowdepth, snow_surft,           &
     rho_snow_ml, rho_snow_grnd, old_frac, new_frac,                 &
     correct_lake_tile, immediate_merge)

USE jules_meltlake_mod, ONLY: nsmax_ml
USE water_constants_mod, ONLY: hcapi, hcapw, tm
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER, INTENT(IN) :: i, land_pts, nsurft
INTEGER, INTENT(IN OUT) :: nsnow(land_pts,nsurft)
REAL(KIND=real_jlslsm), INTENT(IN) :: old_frac, new_frac
LOGICAL, INTENT(IN) :: correct_lake_tile, immediate_merge
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                            &
     sice_ml(land_pts,nsurft,nsmax_ml),                            &
     sliq_ml(land_pts,nsurft,nsmax_ml),                            &
     tsnow_ml(land_pts,nsurft,nsmax_ml),                           &
     ds_ml(land_pts,nsurft,nsmax_ml),                              &
     snowdepth(land_pts,nsurft),                                   &
     snow_surft(land_pts,nsurft),                                  &
     rho_snow_ml(land_pts,nsurft,nsmax_ml),                         &
     rho_snow_grnd(land_pts,nsurft)

INTEGER, PARAMETER :: lake_ml=9, nonlake_ml=10
INTEGER :: n, nl, nn, target, nmerge, source, destination
REAL(KIND=real_jlslsm) :: df, ft, dmi, dmw, de, el, en
REAL(KIND=real_jlslsm) :: mi, mw, new_mi, new_mw, fi, fw
REAL(KIND=real_jlslsm) :: old_mass, new_mass, cap, total_cap, dt
REAL(KIND=real_jlslsm) :: layer_energy(nsmax_ml)
REAL(KIND=real_jlslsm) :: si_l, sl_l, d_l, u_l
REAL(KIND=real_jlslsm) :: si_n, sl_n, d_n, u_n
REAL(KIND=real_jlslsm) :: u_j, heatcap, w_old, w_added, w_new

df = new_frac - old_frac
IF (df == 0.0) RETURN

!-----------------------------------------------------------------------------
! number of active snow layers in lake, non lake snowpacks 
!-----------------------------------------------------------------------------
nl = nsnow(i,lake_ml)
nn = nsnow(i,nonlake_ml)

!-----------------------------------------------------------------------------
! Mix into the growing tile. The shrinking tile is left unchanged.
! During expansion nonlake snow joins the lake snowpack.
! During contraction, including freeze-over, lake snow joins nonlake.
!-----------------------------------------------------------------------------

IF (immediate_merge .AND. correct_lake_tile) THEN

   IF (df > 0.0) THEN
      source      = nonlake_ml
      destination = lake_ml
      w_old       = old_frac
      w_added     = df
      w_new       = new_frac
   ELSE
      source      = lake_ml
      destination = nonlake_ml
      w_old       = 1.0 - old_frac
      w_added     = -df
      w_new       = 1.0 - new_frac
   END IF

   IF (w_new <= 0.0) &
      STOP 'CORRECT_SNOWPACK_FRAC: zero destination fraction'

   nl = nsnow(i,destination)
   nn = nsnow(i,source)
   nmerge = MAX(nl,nn)

   DO n = 1, nmerge

      si_l = 0.0
      sl_l = 0.0
      d_l  = 0.0
      u_l  = 0.0

      si_n = 0.0
      sl_n = 0.0
      d_n  = 0.0
      u_n  = 0.0

      IF (n <= nl) THEN
         si_l = sice_ml(i,destination,n)
         sl_l = sliq_ml(i,destination,n)
         d_l  = ds_ml(i,destination,n)
         u_l = (si_l*hcapi + sl_l*hcapw) * &
               (tsnow_ml(i,destination,n)-tm)
      END IF

      IF (n <= nn) THEN
         si_n = sice_ml(i,source,n)
         sl_n = sliq_ml(i,source,n)
         d_n  = ds_ml(i,source,n)
         u_n = (si_n*hcapi + sl_n*hcapw) * &
               (tsnow_ml(i,source,n)-tm)
      END IF

      u_j = (w_old*u_l + w_added*u_n)/w_new

      sice_ml(i,destination,n) = &
           (w_old*si_l + w_added*si_n)/w_new
      sliq_ml(i,destination,n) = &
           (w_old*sl_l + w_added*sl_n)/w_new
      ds_ml(i,destination,n) = &
           (w_old*d_l + w_added*d_n)/w_new

      heatcap = sice_ml(i,destination,n)*hcapi + &
                sliq_ml(i,destination,n)*hcapw

      IF (heatcap > 0.0) THEN
         tsnow_ml(i,destination,n) = tm + u_j/heatcap
      ELSE
         tsnow_ml(i,destination,n) = tm
      END IF

      IF (ds_ml(i,destination,n) > 0.0) THEN
         rho_snow_ml(i,destination,n) = &
              (sice_ml(i,destination,n) + &
               sliq_ml(i,destination,n)) / ds_ml(i,destination,n)
      ELSE
         rho_snow_ml(i,destination,n) = 0.0
      END IF

   END DO

   nsnow(i,destination) = nmerge
   snowdepth(i,destination) = SUM(ds_ml(i,destination,1:nmerge))
   snow_surft(i,destination) = &
        SUM(sice_ml(i,destination,1:nmerge)) + &
        SUM(sliq_ml(i,destination,1:nmerge))

   IF (snowdepth(i,destination) > 0.0) THEN
      rho_snow_grnd(i,destination) = &
           snow_surft(i,destination)/snowdepth(i,destination)
   END IF

   RETURN
END IF

dmi = df * (SUM(sice_ml(i,nonlake_ml,1:nn)) - &
            SUM(sice_ml(i,lake_ml,1:nl)))
dmw = df * (SUM(sliq_ml(i,nonlake_ml,1:nn)) - &
            SUM(sliq_ml(i,lake_ml,1:nl)))

el = SUM((hcapi*sice_ml(i,lake_ml,1:nl) + &
          hcapw*sliq_ml(i,lake_ml,1:nl)) * &
         (tsnow_ml(i,lake_ml,1:nl) - tm))
en = SUM((hcapi*sice_ml(i,nonlake_ml,1:nn) + &
          hcapw*sliq_ml(i,nonlake_ml,1:nn)) * &
         (tsnow_ml(i,nonlake_ml,1:nn) - tm))
de = df * (en - el)

! Select tile receiving the correction; imbalance calculations above
! remain independent of the target tile.
target = nonlake_ml
ft = 1.0 - new_frac
IF (correct_lake_tile .AND. df > 0.0) THEN
   target = lake_ml
   ft = new_frac
END IF

IF (ft <= 0.0 .OR. nsnow(i,target) <= 0) THEN
   PRINT *, 'CORRECT_SNOWPACK_FRAC: invalid target tile'
   PRINT *, 'Grid point = ', i
   PRINT *, 'Target tile = ', target
   PRINT *, 'Old lake fraction = ', old_frac
   PRINT *, 'New lake fraction = ', new_frac
   PRINT *, 'Target fraction = ', ft
   PRINT *, 'Target nsnow = ', nsnow(i,target)
   STOP 'CORRECT_SNOWPACK_FRAC: invalid target tile'
END IF

mi = SUM(sice_ml(i,target,1:nsnow(i,target)))
mw = SUM(sliq_ml(i,target,1:nsnow(i,target)))
new_mi = mi + dmi/ft
new_mw = mw + dmw/ft

IF (new_mi < 0.0 .OR. new_mw < 0.0) THEN
   STOP 'CORRECT_SNOWPACK_FRAC: negative snow mass'
END IF

fi = 0.0
fw = 0.0
IF (mi > 0.0) fi = new_mi/mi
IF (mw > 0.0) fw = new_mw/mw

! Save layer energies, then change layer masses and thicknesses.
DO n = 1, nsnow(i,target)
   old_mass = sice_ml(i,target,n) + sliq_ml(i,target,n)
   layer_energy(n) = (hcapi*sice_ml(i,target,n) + &
                      hcapw*sliq_ml(i,target,n)) * &
                     (tsnow_ml(i,target,n) - tm)

   sice_ml(i,target,n) = sice_ml(i,target,n)*fi
   sliq_ml(i,target,n) = sliq_ml(i,target,n)*fw
   IF (n == 1) THEN
      IF (mi == 0.0) sice_ml(i,target,n) = new_mi
      IF (mw == 0.0) sliq_ml(i,target,n) = new_mw
   END IF

   new_mass = sice_ml(i,target,n) + sliq_ml(i,target,n)

!-----------------------------------------------------------------------------
! scale layer thickness by *new_mass/old_mass to preserve bulk density profile
!-----------------------------------------------------------------------------     

   IF (old_mass > 0.0) THEN
      ds_ml(i,target,n) = ds_ml(i,target,n)*new_mass/old_mass
   ELSE IF (new_mass > 0.0) THEN
      STOP 'CORRECT_SNOWPACK_FRAC: mass added to empty layer'
   END IF

!-----------------------------------------------------------------------------
! recalcaulte bulk density using scaled ds 
!-----------------------------------------------------------------------------    
   IF (ds_ml(i,target,n) > 0.0) &
      rho_snow_ml(i,target,n) = new_mass/ds_ml(i,target,n)
END DO

! Restore original layer energies and add the gridbox energy correction.
total_cap = SUM(hcapi*sice_ml(i,target,1:nsnow(i,target)) + &
                hcapw*sliq_ml(i,target,1:nsnow(i,target)))
IF (total_cap <= 0.0) STOP 'CORRECT_SNOWPACK_FRAC: zero heat capacity'
dt = (de/ft)/total_cap

DO n = 1, nsnow(i,target)
   cap = hcapi*sice_ml(i,target,n) + hcapw*sliq_ml(i,target,n)
   IF (cap > 0.0) THEN
      tsnow_ml(i,target,n) = tm + layer_energy(n)/cap + dt
   ELSE IF (layer_energy(n) /= 0.0) THEN
      STOP 'CORRECT_SNOWPACK_FRAC: lost layer energy'
   END IF
END DO

snowdepth(i,target) = SUM(ds_ml(i,target,1:nsnow(i,target)))
snow_surft(i,target) = SUM(sice_ml(i,target,1:nsnow(i,target))) + &
                       SUM(sliq_ml(i,target,1:nsnow(i,target)))
IF (snowdepth(i,target) > 0.0) &
   rho_snow_grnd(i,target) = snow_surft(i,target)/snowdepth(i,target)

END SUBROUTINE correct_snowpack_frac

END MODULE expand_area_mod
