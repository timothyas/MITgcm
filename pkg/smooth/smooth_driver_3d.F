#include "SMOOTH_OPTIONS.h"

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
c#include "tamc.h"
#include "SMOOTH.h"
C this seems like bad form, but
C 1. this routine is called by the control package
C 2. there should not be more than one preprocessing option
#ifdef ALLOW_CTRL
# include "CTRL_SIZE.h"
#else
      integer  maxCtrlProc
      parameter ( maxCtrlProc = 1 )
#endif


C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      _RL smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      character*(20) smooth_options(maxCtrlProc)
      integer smoothOpNb
      integer myThid

C --- Local Parameters
      character*( 80) fnamegeneric
      integer k2
      logical doNormalize

c read smoothing [i.e diffusion] operator and normalization:
      CALL smooth_read_operator_3d(smoothOpNb,myThid)

      doNormalize=.TRUE.
      do k2 = 1, maxCtrlProc
        if (smooth_options(k2).eq.'WC01')
     &      CALL smooth_correl3d(fld_in,smoothOpNb,myThid)
        if (smooth_options(k2).eq.'smooth')
     &      CALL smooth_diff3d(fld_in,smooth3Dnbt(smoothOpNb),myThid)
        if (smooth_options(k2).eq.'matern_precond')
     &      CALL smooth_inverse_matern_3d(fld_in,
     &          smooth3DMask, smoothOpNb, doNormalize, myThid)
      enddo

      CALL EXCH_XYZ_RL ( fld_in , myThid )

      end
