#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE init_meltlake_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_meltlake(nml_dir)
!-----------------------------------------------------------------------------
! Description:
!   Reads in the meltlake namelist items and checks them for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
USE io_constants, ONLY: namelist_unit

USE string_utils_mod, ONLY: to_string

USE jules_meltlake_mod, ONLY: jules_meltlake, check_jules_meltlake, nsmax_ml,  &
   print_nlist_jules_meltlake

USE logging_mod, ONLY: log_info, log_fatal

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
CHARACTER(LEN=*), INTENT(IN) :: nml_dir  ! The directory containing the
                                         ! namelists
! Work variables
INTEGER :: ERROR  ! Error indicator
CHARACTER(LEN=errormessagelength) :: iomessage

!-----------------------------------------------------------------------------
! First, read the meltlake namelist
!-----------------------------------------------------------------------------
CALL log_info("init_meltlake", "Reading JULES_MELTLAKE namelist...")

OPEN(namelist_unit, FILE=(TRIM(nml_dir) // '/' // 'jules_meltlake.nml'),       &
     STATUS='old', POSITION='rewind', ACTION='read', IOSTAT = ERROR,           &
     IOMSG = iomessage)
IF ( ERROR /= 0 )                                                              &
  CALL log_fatal("init_meltlake",                                              &
                 "Error opening namelist file jules_meltlake.nml " //          &
                 "(IOSTAT=" // TRIM(to_string(ERROR)) // " IOMSG=" //          &
                 TRIM(iomessage) // ")")

READ(namelist_unit, NML = jules_meltlake, IOSTAT = ERROR, IOMSG = iomessage)
IF ( ERROR /= 0 )                                                              &
  CALL log_fatal("init_meltlake",                                              &
                 "Error reading namelist JULES_MELTLAKE " //                   &
                 "(IOSTAT=" // TRIM(to_string(ERROR)) // " IOMSG=" //          &
                 TRIM(iomessage) // ")")

CLOSE(namelist_unit, IOSTAT = ERROR, IOMSG = iomessage)
IF ( ERROR /= 0 )                                                              &
  CALL log_fatal("init_meltlake",                                              &
                 "Error closing namelist file jules_meltlake.nml " //          &
                 "(IOSTAT=" // TRIM(to_string(ERROR)) // " IOMSG=" //          &
                 TRIM(iomessage) // ")")

CALL print_nlist_jules_meltlake()
CALL check_jules_meltlake()



RETURN
END SUBROUTINE init_meltlake

END MODULE init_meltlake_mod
#endif
