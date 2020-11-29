#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

      subroutine smooth_driver_2D (
     U     fld_in, fld_reg,
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
#ifdef ALLOW_CTRL
# include "CTRL_SIZE.h"
#else
      integer  maxCtrlProc
      parameter ( maxCtrlProc = 1 )
#endif


C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fld_reg(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      character*(20) smooth_options(maxCtrlProc)
      _RL smooth2DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      integer myThid

C --- Local Parameters
      character*( 80) fnamegeneric
      integer k2
      logical doNormalize
      logical doWC01, doSmooth, doMaternPre, doMaternReg

c read smoothing [i.e diffusion] operator and normalization:
      CALL smooth_read_operator_2d(smoothOpNb,myThid)

c determine algorithm to use
      doNormalize=.TRUE.
      doWC01=.FALSE.
      doSmooth=.FALSE.
      doMaternPre=.FALSE.
      doMaternReg=.FALSE.
      do k2 = 1, maxCtrlProc
        if (smooth_options(k2)(1:4).eq.'WC01') doWC01 = .TRUE.
        if (smooth_options(k2)(1:6).eq.'smooth') doSmooth = .TRUE.
        if (smooth_options(k2)(1:10).eq.'matern_pre') doMaternPre=.TRUE.
        if (smooth_options(k2)(1:10).eq.'matern_reg') doMaternReg=.TRUE.
      enddo

      if (doWC01)
     &      CALL smooth_correl2d(fld_in,smooth2DMask,smoothOpNb,myThid)
      if (doSmooth)
     &      CALL smooth_diff2d(fld_in,smooth2DMask,
     &          smooth2Dnbt(smoothOpNb),myThid)
      if (doMaternPre)
     &      CALL smooth_inverse_matern_2d(fld_in,
     &          smooth2DMask, smoothOpNb, doNormalize, myThid)
      if (doMaternReg)
     &      CALL smooth_apply_matern_2d(fld_in, fld_reg,
     &          smooth2DMask, smoothOpNb, doNormalize, myThid)

      CALL EXCH_XY_RL ( fld_in , myThid )

      end
