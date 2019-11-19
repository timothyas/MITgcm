#include "AUTODIFF_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_DIAGNOSTICS
# include "DIAG_OPTIONS.h"
#endif
#include "AD_CONFIG.h"

C     !ROUTINE: addummy_for_hfac
C     !INTERFACE:
      subroutine addummy_for_hfac( myname, myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE addummy_for_hfac
C     *==========================================================*
C     Extract adjoint variables associated with r_low sensitivites
C     Including:
C
C       myname  variable(s)
C       ------  -----------
C       'C'     adhFacC
C       'W'     adhFacW
C       'S'     adhFacS
C       'R'     sensitivities to hFac reciprocals:
C                   adrecip_hFacC, adrecip_hFacW, adrecip_hFacS
C       'A'     A matrix components in CG2D solve:
C                   adAC2d, adAW2d, adAS2d
C      
C       'O'     write filled diagnostics out to file
C
C     Modified from verification/bottom_ctrl_5x5/code_ad/dummy_in_hfac
C     by Tim Smith
C     *==========================================================*
C     | SUBROUTINE addummy_for_hfac
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "AUTODIFF_PARAMS.h"
#ifdef ALLOW_AUTODIFF_MONITOR
# include "AUTODIFF.h"
# include "adcommon.h"
# ifdef ALLOW_DIAGNOSTICS
#  include "DIAGNOSTICS_SIZE.h"
#  include "DIAGNOSTICS.h"
# endif
#endif

      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
      INTEGER  IO_ERRCOUNT
      EXTERNAL IO_ERRCOUNT

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid - Thread number for this instance of the routine.
      CHARACTER*(1) myname
      _RL myTime
      integer myIter
      integer myThid

#ifdef ALLOW_AUTODIFF_MONITOR
C     !LOCAL VARIABLES:
c     == local variables ==
C     suff - Hold suffix part of a filename
C     beginIOErrCount - Begin and end IO error counts
C     endIOErrCount
C     msgBuf - Error message buffer
      CHARACTER*(MAX_LEN_FNAM) suff
      INTEGER beginIOErrCount
      INTEGER endIOErrCount
      CHARACTER*(MAX_LEN_MBUF) msgBuf
# ifdef ALLOW_DIAGNOSTICS
      logical modelEnd, modelStart
# endif

c     == end of interface ==
CEOP

#ifdef ALLOW_DEPTH_CONTROL

      if (
     & different_multiple(adjDumpFreq, myTime, deltaTClock)) then

#ifdef ALLOW_DIAGNOSTICS
        if ( useDiag4AdjOutp ) then
          modelEnd = myTime.eq.endTime .or. myIter.eq.nEndIter
          modelStart = myTime.eq.startTime .or. myIter.eq.nIter0
          call diagnostics_switch_onoff( myTime, myIter, myThid )
        endif
#endif

        CALL TIMER_START('I/O (WRITE)        [ADJOINT LOOP]', myThid )

C--     Set suffix for this set of data files.
        WRITE(suff,'(I10.10)') myIter
        writeBinaryPrec = writeStatePrec

C--     Read IO error counter
        beginIOErrCount = IO_ERRCOUNT(myThid)

        IF ( myname(1:1) .eq. 'C' ) THEN
         CALL WRITE_FLD_XYZ_RL('ADJhFacC.',suff,adhfacc,myIter,myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL(adhFacC,'ADJhFacC', 0,Nr,0,1,1,myThid)
#endif

        ELSE IF ( myname(1:1) .eq. 'W' ) THEN
         CALL WRITE_FLD_XYZ_RL('ADJhFacW.',suff,adhfacw,myIter,myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL(adhFacW,'ADJhFacW',0,Nr,0,1,1,myThid)
#endif

        ELSE IF ( myname(1:1) .eq. 'S' ) THEN
         CALL WRITE_FLD_XYZ_RL ('ADJhFacS.',suff,adhfacs,myIter,myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL(adhFacS,'ADJhFacS',0,Nr,0,1,1,myThid)
#endif
        else if (myname(1:1) .eq. 'R' ) then
         CALL WRITE_FLD_XYZ_RL('ADJrhFaC.', suff, adrecip_hFacC,
     &                         myIter, myThid)
         CALL WRITE_FLD_XYZ_RL('ADJrhFaW.', suff, adrecip_hFacW,
     &                         myIter, myThid)
         CALL WRITE_FLD_XYZ_RL('ADJrhFaS.', suff, adrecip_hFacS,
     &                         myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL(adrecip_hFacC,'ADJrhFaC',
     &                         0,Nr,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL( adrecip_hFacW, 'ADJrhFaW',
     &                         0,Nr,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL( adrecip_hFacS, 'ADJrhFaS',
     &                         0,Nr,0,1,1,myThid)
#endif
        else if (myname(1:1) .eq. 'A' ) then
         CALL WRITE_FLD_XY_RL('ADJaC2d.', suff, adac2d,myIter, myThid)
         CALL WRITE_FLD_XY_RL('ADJaW2d.', suff, adaw2d,myIter, myThid)
         CALL WRITE_FLD_XY_RL('ADJaS2d.', suff, adas2d,myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL(adac2d,'ADJaC2d',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(adaw2d,'ADJaW2d',0,1,0,1,1,myThid)
         CALL DIAGNOSTICS_FILL(adas2d,'ADJaS2d',0,1,0,1,1,myThid)
#endif

#ifdef ALLOW_DIAGNOSTICS
C--
C   For now, this is a hack to print these variables through
C   the diagnostics package. I don't want the write statement
C   to conflict with addummy_in_stepping.
C   This only happens for iter0 (the final adjoint step)
C   and at this stage I want to print everything that's been filled.
       else if (myname(1:1) .eq. 'O') then
        if ( modelStart )
     &   call diagnostics_write_adj(modelStart, myTime, myIter, myThid)
#endif
       ELSE
         WRITE(msgBuf,'(A,I10)') 'addummy_for_hfac: ',
     &      'no valid input for myname specified'
         CALL PRINT_ERROR( msgBuf, 1 )
       END IF
C--     Reread IO error counter
        endIOErrCount = IO_ERRCOUNT(myThid)

C--     Check for IO errors
        IF ( endIOErrCount .NE. beginIOErrCount ) THEN
         WRITE(msgBuf,'(A)')  'S/R WRITE_STATE'
         CALL PRINT_ERROR( msgBuf, 1 )
         WRITE(msgBuf,'(A)')  'Error writing out model state'
         CALL PRINT_ERROR( msgBuf, 1 )
         WRITE(msgBuf,'(A,I10)') 'Timestep ',myIter
         CALL PRINT_ERROR( msgBuf, 1 )
        ENDIF
        CALL TIMER_STOP( 'I/O (WRITE)        [ADJOINT LOOP]', myThid )
      endif ! different_multiple

#endif /* ALLOW_DEPTH_CONTROL */
#endif /* ALLOW_AUTODIFF_MONITOR */

      RETURN
      END
