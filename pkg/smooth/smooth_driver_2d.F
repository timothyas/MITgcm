#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

      subroutine smooth_driver_2D (
     U     fld_in,
     I     smooth_options,smooth2DMask,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_driver_2d
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
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      character*(20) smooth_options
      _RL smooth2DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      integer myThid

C --- Local Parameters
      logical doNormalize
      logical doWC01, doSmooth, doMaternPre, doMaternReg

      doNormalize=.TRUE.
c read smoothing [i.e diffusion] operator and normalization:
      CALL smooth_read_operator_2d(smoothOpNb,myThid)

      if (smooth_options(1:4).eq.'WC01')
     &      CALL smooth_correl2d(fld_in,smooth2DMask,smoothOpNb,myThid)

      if (smooth_options(1:6).eq.'smooth')
     &      CALL smooth_diff2d(fld_in,smooth2DMask,
     &          smooth2Dnbt(smoothOpNb),myThid)

      if (smooth_options(1:10).eq.'matern_pre')
     &      CALL smooth_inverse_matern_2d(fld_in,
     &          smooth2DMask, smoothOpNb, doNormalize, myThid)

      CALL EXCH_XY_RL ( fld_in , myThid )

      end
