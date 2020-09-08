#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_SHELFICE
# include "SHELFICE_OPTIONS.h"
#endif
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
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

C --- Make sure 2D Dims are 'xy', 'xz', or 'yz'
        IF ((smooth2DDims(smoothOpNb)(1:2).NE.'xy').AND.
     &      (smooth2DDims(smoothOpNb)(1:2).NE.'xz').AND.
     &      (smooth2DDims(smoothOpNb)(1:2).NE.'yz')) THEN
          WRITE(msgBuf,'(2A,I2,2A)') 'SMOOTH_CHECK: ',
     &      'Unrecognized option for smooth2DDims(', smoothOpNb,
     &      '): ', smooth2DDims(smoothOpNb)
          CALL PRINT_ERROR( msgBuf, myThid )
          WRITE(msgBuf,'(2A)') 'SMOOTH_CHECK: ',
     &      'must be: xy , xzn, xzs, yze, or yzw (lower case)'
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_CHECK'
        ENDIF

C --- Make sure correct OBCS flags turned on for smoothing xz or yz
C     Also, not implemented for WC01 algorithm ... yet?
        IF ((smooth2DDims(smoothOpNb)(1:2).EQ.'xz').OR.
     &      (smooth2DDims(smoothOpNb)(1:2).EQ.'yz')) THEN
#ifndef ALLOW_OBCS
          WRITE(msgBuf,'(A)') 'SMOOTH_CHECK: ',
     &      'Must compile with ALLOW_OBCS to smooth xz/yz slice'
          CALL PRINT_ERROR( msgBuf, myThid )
#endif
          IF (smooth2DAlgorithm(smoothOpNb)(1:6).NE.'matern') THEN
            WRITE(msgBuf,'(3A)') 'SMOOTH_CHECK: ',
     &        'smoothing of xz/yz slices only implemented for matern ',
     &        'algorithm'
            CALL PRINT_ERROR( msgBuf, myThid )
          ENDIF
        ENDIF


        IF (smooth2DDims(smoothOpNb)(1:2).EQ.'xz') THEN
#if !( defined(ALLOW_OBCS_NORTH) || defined(ALLOW_OBCS_SOUTH) )
          WRITE(msgBuf,'(3A)') 'SMOOTH_CHECK: ',
     &      'Must compile with ALLOW_OBCS_NORTH or ALLOW_OBCS_SOUTH',
     &      ' to smooth xz slice'
          CALL PRINT_ERROR( msgBuf, myThid )
#endif
          IF ((smooth2DDims(smoothOpNb)(3:3).NE.'n').AND.
     &        (smooth2DDims(smoothOpNb)(3:3).NE.'s')) THEN
          WRITE(msgBuf,'(A,A,I,A)') 'SMOOTH_CHECK: ',
     &      'Must specify north or south for xz slice as',
     &      'smooth2DDims(',smoothOpNb,') = xzn or xzs'
          CALL PRINT_ERROR( msgBuf, myThid )
          ENDIF
        ENDIF

        IF (smooth2DDims(smoothOpNb)(1:2).EQ.'yz') THEN
#if !( defined(ALLOW_OBCS_EAST) || defined(ALLOW_OBCS_WEST) )
          WRITE(msgBuf,'(3A)') 'SMOOTH_CHECK: ',
     &      'Must compile with ALLOW_OBCS_WEST or ALLOW_OBCS_WEST',
     &      ' to smooth yz slice'
          CALL PRINT_ERROR( msgBuf, myThid )
#endif
          IF ((smooth2DDims(smoothOpNb)(3:3).NE.'e').AND.
     &        (smooth2DDims(smoothOpNb)(3:3).NE.'w')) THEN
          WRITE(msgBuf,'(A,A,I,A)') 'SMOOTH_CHECK: ',
     &      'Must specify east or west for yz slice as',
     &      'smooth2DDims(',smoothOpNb,') = yze or yzw'
          CALL PRINT_ERROR( msgBuf, myThid )
          ENDIF
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
