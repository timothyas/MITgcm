#include "AUTODIFF_OPTIONS.h"

      subroutine g_dummy_for_hfac( myname, myTime, myIter, myThid )

C     *==========================================================*
C     | SUBROUTINE g_dummy_for_hfac
C     *==========================================================*

      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     == Routine arguments ==
C     myThid - Thread number for this instance of the routine.
      CHARACTER*(1) myname
      _RL myTime
      INTEGER myIter
      INTEGER myThid

      RETURN
      END
