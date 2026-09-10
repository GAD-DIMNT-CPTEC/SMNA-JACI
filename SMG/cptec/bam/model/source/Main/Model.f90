
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.28 $
!
PROGRAM Main
 USE AtmosModelMod, Only : atmos_model_init,atmos_model_run,atmos_model_finalize
 IMPLICIT NONE
 CALL atmos_model_init()
 CALL atmos_model_run()
 CALL atmos_model_finalize()
 STOP
END PROGRAM Main
