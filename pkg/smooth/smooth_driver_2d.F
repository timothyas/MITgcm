#include "SMOOTH_OPTIONS.h"

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
c#include "tamc.h"
#include "SMOOTH.h"
C this seems like bad form, but
C 1. this routine is called by the control package
C 2. there should not be more than one preprocessing option
#include "CTRL_SIZE.h"


C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smooth2DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      character*(20) smooth_options(maxCtrlProc)
      integer smoothOpNb
      integer myThid

C --- Local Parameters
      character*( 80) fnamegeneric
      integer k2
      logical doNormalize

c read smoothing [i.e diffusion] operator and normalization:
      CALL smooth_read_operator_2d(smoothOpNb,myThid)

      doNormalize=.TRUE.
      do k2 = 1, maxCtrlProc
        if (smooth_options(k2).eq.'WC01')
     &      CALL smooth_correl2d(fld_in,smooth2DMask,smoothOpNb,myThid)
        if (smooth_options(k2).eq.'smooth')
     &      CALL smooth_diff2d(fld_in,smooth2DMask,
     &          smooth2Dnbt(smoothOpNb),myThid)
        if (smooth_options(k2).eq.'matern_precond')
     &      CALL smooth_inverse_matern_2d(fld_in,
     &          smooth2DMask, smoothOpNb, doNormalize, myThid)
      enddo

      CALL EXCH_XY_RL ( fld_in , myThid )

      end
