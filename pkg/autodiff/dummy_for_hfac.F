#include "AUTODIFF_OPTIONS.h"

      subroutine dummy_for_hfac( myname, mytime, myIter, myThid )

C     *==========================================================*
C     | SUBROUTINE dummy_for_hfac
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
