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
C     Extract adjoint variable from TAMC/TAF-generated
C     adjoint common blocks, contained in adcommon.h
C     and write fields to file;
C     Make sure common blocks in adcommon.h are up-to-date
C     w.r.t. current adjoint code.
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
      CHARACTER*(30) myname
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
      CHARACTER*(30) myfullname
# ifdef ALLOW_DIAGNOSTICS
      logical modelEnd, modelStart
# endif

c     == end of interface ==
CEOP

#ifdef ALLOW_DEPTH_CONTROL
        WRITE(msgBuf,'(2A)')  'addummy_for_hfac, myname ',
     &  myname
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )

      if (
     & different_multiple(adjDumpFreq, myTime, deltaTClock)) then

#ifdef ALLOW_DIAGNOSTICS
        if ( useDiag4AdjOutp ) then
          modelEnd = myTime.eq.endTime .or. myIter.eq.nEndIter
          modelStart = myTime.eq.startTime .or. myIter.eq.nIter0
          call diagnostics_switch_onoff( myTime, myIter, myThid )
        endif
#endif

        WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac '
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )

        CALL TIMER_START('I/O (WRITE)        [ADJOINT LOOP]', myThid )

C--     Set suffix for this set of data files.
        WRITE(suff,'(I10.10)') myIter
        writeBinaryPrec = writeStatePrec

C--     Read IO error counter
        beginIOErrCount = IO_ERRCOUNT(myThid)
        WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac, if statement 1 '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )

        IF ( myname(1:1) .eq. 'C' ) THEN
         WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac, if hFacC!! '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )
           myfullname = 'hFacC'
           CALL WRITE_FLD_XYZ_RL ( 'ADJhFacC.', suff, adhfacc,
     &          myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( adhFacC, 'ADJhFacC', 0,Nr,0,1,1,myThid )
#endif

        ELSE IF ( myname(1:1) .eq. 'W' ) THEN
         WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac, if hFacW!! '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )
           myfullname = 'hFacW'
           CALL WRITE_FLD_XYZ_RL ( 'ADJhFacW.', suff, adhfacw,
     &          myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( adhFacW, 'ADJhFacW', 0,Nr,0,1,1,myThid )
#endif
        ELSE IF ( myname(1:1) .eq. 'S' ) THEN
         WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac, if hFacS!! '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )
           myfullname = 'hFacS'
           CALL WRITE_FLD_XYZ_RL ( 'ADJhFacS.', suff, adhfacs,
     &          myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( adhFacS, 'ADJhFacS', 0,Nr,0,1,1,myThid )
#endif
        else if (myname(1:1) .eq. 'R' ) then
         WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac, reciprocals!! '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )

           CALL WRITE_FLD_XYZ_RL ( 'ADJrhFaC.', suff, adrecip_hFacC,
     &          myIter, myThid)
           CALL WRITE_FLD_XYZ_RL ( 'ADJrhFaW.', suff, adrecip_hFacW,
     &          myIter, myThid)
           CALL WRITE_FLD_XYZ_RL ( 'ADJrhFaS.', suff, adrecip_hFacS,
     &          myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( adrecip_hFacC, 'ADJrhFaC', 0,Nr,0,1,1,
     &                        myThid )
         CALL DIAGNOSTICS_FILL( adrecip_hFacW, 'ADJrhFaW', 0,Nr,0,1,1,
     &                        myThid )
         CALL DIAGNOSTICS_FILL( adrecip_hFacS, 'ADJrhFaS', 0,Nr,0,1,1,
     &                        myThid )
#endif
        else if (myname(1:1) .eq. 'A' ) then
         WRITE(msgBuf,'(A)')  ' inside addummy_for_hfac, A matrix!! '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, myThid )

           CALL WRITE_FLD_XY_RL ( 'ADJaC2d.', suff, adac2d,
     &          myIter, myThid)
           CALL WRITE_FLD_XY_RL ( 'ADJaW2d.', suff, adaw2d,
     &          myIter, myThid)
           CALL WRITE_FLD_XY_RL ( 'ADJaS2d.', suff, adas2d,
     &          myIter, myThid)
#ifdef ALLOW_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( adac2d, 'ADJaC2d', 0,1,0,1,1,
     &                        myThid )
         CALL DIAGNOSTICS_FILL( adaw2d, 'ADJaW2d', 0,1,0,1,1,
     &                        myThid )
         CALL DIAGNOSTICS_FILL( adas2d, 'ADJaS2d', 0,1,0,1,1,
     &                        myThid )
#endif

#ifdef ALLOW_DIAGNOSTICS
C--
C   For now, this is a hack to print these variables through
C   the diagnostics package. I don't want the write statement
C   to conflict with addummy_in_stepping.
C   This only happens for iter0 (the final adjoint step)
C   and at this stage I want to print everything that's been filled.
       else if (myname(1:17) .eq. 'writeAdjHfacDiags') then
        if ( modelStart )
     &   call diagnostics_write_adj(modelStart, myTime, myIter, myThid)
#endif
       ELSE
          write(*,*) 'addummy_for_hfac: no valid myname specified'
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
        ELSE
         WRITE(msgBuf,'(A,I10)')
     &    '// ad'//myfullname//' written, timestep', myIter
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, 1 )
         WRITE(msgBuf,'(A)')  ' '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &    SQUEEZE_RIGHT, 1 )
        ENDIF

        CALL TIMER_STOP( 'I/O (WRITE)        [ADJOINT LOOP]', myThid )
      endif ! different_multiple

#endif /* ALLOW_DEPTH_CONTROL */
#endif /* ALLOW_AUTODIFF_MONITOR */

      RETURN
      END
