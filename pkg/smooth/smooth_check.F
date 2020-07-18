#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_SHELFICE
# include "SHELFICE_OPTIONS.h"
#endif

      subroutine smooth_check( myThid )

C     *==========================================================*
C     | SUBROUTINE smooth_check
C     | o Routine that checks for inconsistencies in the smooth package
C     *==========================================================*

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "PARAMS.h"

#include "SMOOTH.h"

c     == routine arguments ==

      integer myThid

c     == local variables ==
C     msgBuf      - Informational/error message buffer
C     iUnit       - Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      integer smoothOpNb

c     == end of interface ==

C     TODO:
C --- Check to make sure operator file exists in this case
C       when smooth2Dtype .ne.0 and smooth2Dsize.eq.1
C       when smooth3Dtype .ne.0 and smooth3Dsize.eq.1 or 2

C --- Check to make sure length scale file exists in this case
C       when smooth2Dtype .ne.0 and smooth2Dsize.eq.2
C       when smooth3Dtype .ne.0 and smooth3Dsize.eq.3

C --- Check to make sure norm exists if 
C       smooth2dfilter.eq.1
C       smooth3dfilter.eq.1
      

C --- Check for legal mask names
      DO smoothOpNb=1,smoothOpNbMax
        IF ( (smooth2DmaskName(smoothOpNb).ne.'maskC').and.
     &       (smooth2DmaskName(smoothOpNb).ne.'maskW').and.
     &       (smooth2DmaskName(smoothOpNb).ne.'maskS').and.
     &       (smooth2DmaskName(smoothOpNb).ne.'maskI') ) then
          WRITE(msgBuf,'(2A,I2,2A)') 'SMOOTH_CHECK: ',
     &      'Unrecognized option for smooth2DmaskName(', smoothOpNb, 
     &      '): ', smooth2DmaskName(smoothOpNb)
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'SMOOTH_CHECK: ',
     &      'must be: maskC, maskW, maskS, or maskI'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_CHECK'
        ENDIF

C --- only maskC is allowed for smoothing 3D fields... for now I guess.
        IF ( (smooth3DmaskName(smoothOpNb).ne.'maskC') ) then
          WRITE(msgBuf,'(2A)') 'SMOOTH_CHECK: ',
     &      'Only maskC implemented for 3D smoothing'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_CHECK'
        ENDIF

        IF ( (smooth3DmaskName(smoothOpNb).ne.'maskC').and.
     &       (smooth3DmaskName(smoothOpNb).ne.'maskW').and.
     &       (smooth3DmaskName(smoothOpNb).ne.'maskS')) then
          WRITE(msgBuf,'(2A,I2,2A)') 'SMOOTH_CHECK: ',
     &      'Unrecognized option for smooth3DmaskName(', smoothOpNb, 
     &      '): ', smooth3DmaskName(smoothOpNb)
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'SMOOTH_CHECK: ',
     &      'must be: maskC, maskW, or maskS'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_CHECK'
        ENDIF



C --- make sure proper macros are defined for maskI (iceshelf domain)
#if !( defined(ALLOW_SHI2D_CONTROL) || defined(ALLOW_SHIFWFLX_CONTROL) )
        IF ( (smooth2DmaskName(smoothOpNb).eq.'maskI') .or. 
     &       (smooth3DmaskName(smoothOpNb).eq.'maskI') )then
          WRITE(msgBuf,'(3A)') 'SMOOTH_CHECK: ',
     &      'Must compile with option ALLOW_SHI2D_CONTROL or ',
     &      'ALLOW_SHIFWFLX_CONTROL'
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(3A)') 'SMOOTH_CHECK: ',
     &      ' in SHELFICE_OPTIONS.h ',
     &      '(and CTRL_OPTIONS.h) to define maskI'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_CHECK'
        ENDIF
#endif
      ENDDO

      RETURN
      END
