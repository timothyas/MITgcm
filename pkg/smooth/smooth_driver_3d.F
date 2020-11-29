#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

      subroutine smooth_driver_3D (
     U     fld_in, fld_reg,
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
#include "tamc.h"
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
      _RL fld_reg(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      character*(20) smooth_options(maxCtrlProc)
      _RL smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      integer smoothOpNb
      integer myThid

C --- Local Parameters
      character*( 80) fnamegeneric
      integer k2
      logical doNormalize
      logical doWC01, doSmooth, doMaternPre, doMaternReg

c read smoothing [i.e diffusion] operator and normalization:
      CALL smooth_read_operator_3d(smoothOpNb,myThid)

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
     &    CALL smooth_correl3d(fld_in,smoothOpNb,myThid)
      if (doSmooth)
     &    CALL smooth_diff3d(fld_in,smooth3Dnbt(smoothOpNb),myThid)
      if (doMaternPre)
     &    CALL smooth_inverse_matern_3d(fld_in,
     &        smooth3DMask, smoothOpNb, doNormalize, myThid)
      if (doMaternReg)
     &    CALL smooth_apply_matern_3d(fld_in, fld_reg,
     &        smooth3DMask, smoothOpNb, doNormalize, myThid)

      CALL EXCH_XYZ_RL ( fld_in , myThid )

      end
