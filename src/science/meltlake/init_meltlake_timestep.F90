! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
!
! Code Description:
!   Language: FORTRAN 90
!   Initialise varaible every timestep.
!   Note: This subroutine is doing nothing. Keep in case it's useful later 
! Code Owner: Please refer to ModuleLeaders.txt
!

MODULE init_meltlake_timestep_mod

USE parkind1, ONLY: jpim

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = "INIT_MELTLAKE_TIMESTEP_MOD"

CONTAINS

SUBROUTINE init_meltlake_timestep(meltlake_vars)
  
USE meltlake_vars_mod, ONLY: meltlake_vars_type
USE yomhook, ONLY: lhook, dr_hook  
USE parkind1, ONLY: jprb
USE model_time_mod, ONLY: timestep_number

IMPLICIT NONE
!
! Description:
!   Initialise meltlake varaibles every model timestep 
!   
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!

TYPE(meltlake_vars_type), INTENT(IN OUT) :: meltlake_vars

REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*),  PARAMETER :: RoutineName = "INIT_MELTLAKE_TIMESTEP"

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------
! Initialisation
!-------------------------------------------------
!meltlake_vars%lake_inflow(:,:)  = 0.0


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_meltlake_timestep
END MODULE init_meltlake_timestep_mod
