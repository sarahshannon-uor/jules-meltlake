! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_meltlake_mod

USE max_dimensions,    ONLY: snow_layers_max_ml
USE jules_surface_mod, ONLY: l_elev_land_ice
USE jules_science_fixes_mod, ONLY: l_fix_neg_snow
!-----------------------------------------------------------------------------
! Description:
!   Contains snow options and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE missing_data_mod, ONLY: rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER :: nsmax_ml = 200        !  Maximum number of snow layers

REAL(KIND=real_jlslsm) :: rho_firn_efold = 37

REAL(KIND=real_jlslsm) :: firn_depth_max = 35.0

!REAL(KIND=real_jlslsm) :: dzsnow_ml(snow_layers_max_ml) = rmdi
                ! Prescribed thickness of snow layers (m)
                ! This is the thickness of each snow layer when it is not
                ! the bottom layer (note that dzSnow(nsMax) is not used
                ! because that is always the bottom layer)
REAL(KIND=real_jlslsm) :: dzsnow_ml(200) = rmdi
!-----------------------------------------------------------------------------
! Switches
!-----------------------------------------------------------------------------
LOGICAL :: l_meltlake = .FALSE.                                               
               
!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_meltlake/                                                     &
! Switches
     l_meltlake, nsmax_ml, dzsnow_ml, firn_depth_max, rho_firn_efold 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_MELTLAKE_MOD'

CONTAINS

#if !defined(RIVERS_ONLY)
SUBROUTINE check_jules_meltlake()

USE ereport_mod, ONLY: ereport

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_MELTLAKE namelist for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER :: errorstatus
CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_JULES_MELTLAKE'

!-----------------------------------------------------------------------------

! Check that we have a suitable number of snow layers
IF ( nsmax_ml > snow_layers_max_ml ) THEN
  errorstatus = 101
  CALL ereport(RoutineName, errorstatus,                                       &
  "Too many snow layers specified - increase snow_layers_max and recompile")
END IF

! Elevated ice tile is needed for the melt lake model 
IF ( .NOT. l_elev_land_ice) THEN
  errorstatus = 101
  CALL ereport(RoutineName, errorstatus,                                       &
  "Melt lake model requires l_elev_land_ice=.true")
END IF

! fix for negative snow is needed for the melt lake model 
IF ( .NOT. l_fix_neg_snow) THEN
  errorstatus = 101
  CALL ereport(RoutineName, errorstatus,                                       &
  "Melt lake model requires l_fix_neg_snow=.true.")
END IF

END SUBROUTINE check_jules_meltlake
#endif

SUBROUTINE print_nlist_jules_meltlake()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer


!-----------------------------------------------------------------------------


CALL jules_print('jules_meltlake',                                                 &
                 'Contents of namelist jules_meltlake')

WRITE(lineBuffer, *) '  l_meltlake = ', l_meltlake
CALL jules_print('jules_mellake', lineBuffer)

WRITE(lineBuffer, *) '  nsmax_ml = ', nsmax_ml
CALL jules_print('jules_meltlake', lineBuffer)

WRITE(lineBuffer, *) '  dzsnow_ml = ', dzsnow_ml
CALL jules_print('jules_meltlake', lineBuffer)

WRITE(lineBuffer, *) '  firn_depth_max = ', firn_depth_max
CALL jules_print('jules_meltlake', lineBuffer)

WRITE(lineBuffer, *) '  rho_firn_efold = ', rho_firn_efold
CALL jules_print('jules_meltlake', lineBuffer)


CALL jules_print('jules_meltlake',                                            &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_meltlake

#if defined(UM_JULES) && !defined(LFRIC)
SUBROUTINE read_nml_jules_meltlake (unitnumber)

! Description:
!  Read the JULES_MELTLAKE namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype
USE parkind1,         ONLY: jprb,  jpim
USE yomhook,          ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber

INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_SNOW'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 7
INTEGER, PARAMETER :: n_real = 25 + snow_layers_max + 5 * npft_max
INTEGER, PARAMETER :: n_log = 5 + npft_max

TYPE :: my_namelist
  SEQUENCE
  LOGICAL :: l_meltlake
  INTEGER :: nsmax_ml
  REAL(KIND=real_jlslsm) :: firn_depth_max
  REAL(KIND=real_jlslsm) :: rho_firn_efold
  REAL(KIND=real_jlslsm) :: dzsnow_ml(snow_layers_max_ml)
  
 
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,               &
                    n_real_in = n_real, n_log_in = n_log)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_meltlake, IOSTAT = errorstatus,         &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_meltlake", iomessage)

  my_nml % l_meltlake      = l_meltlake
  my_nml % nsmax_ml        = nsmax_ml
  my_nml % dzsnow_ml       = dzsnow_ml
  my_nml % firn_depth_max  = firn_depth_max
  my_nml % rho_firn_efold  = rho_firn_efold
  
  
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  
  l_meltlake             = my_nml % l_meltlake
  nsmax_ml               = my_nml % nsmax_ml
  dzsnow_ml              = my_nml % dzsnow_ml
  firn_depth_max         = my_nml % firn_depth_max
  rho_firn_efold         = my_nml % rho_firn_efold

  
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_meltlake
#endif

END MODULE jules_meltlake_mod
