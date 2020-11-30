#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

      subroutine smooth_driver_3D (
     U     fld_in,
     I     smooth_options,smooth3DMask,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_driver_3d
C     | o determine which smoothing algorithm to employ
C     | o read in the operator
C     | o be cool, be smooth
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      character*(20) smooth_options
      _RL smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      integer smoothOpNb
      integer myThid

C --- Local Parameters
      logical doNormalize

      doNormalize=.TRUE.
c read smoothing [i.e diffusion] operator and normalization:
      CALL smooth_read_operator_3d(smoothOpNb,myThid)

      if (smooth_options(1:4).eq.'WC01')
     &    CALL smooth_correl3d(fld_in,smoothOpNb,myThid)

      if (smooth_options(1:6).eq.'smooth')
     &    CALL smooth_diff3d(fld_in,smooth3Dnbt(smoothOpNb),myThid)

      if (smooth_options(1:10).eq.'matern_pre')
     &    CALL smooth_inverse_matern_3d(fld_in,
     &        smooth3DMask, smoothOpNb, doNormalize, myThid)

      CALL EXCH_XYZ_RL ( fld_in , myThid )

      end
