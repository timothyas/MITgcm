#include "CTRL_OPTIONS.h"

CBOP
C     !ROUTINE: ctrl_depth_ini
C     !INTERFACE:
      subroutine ctrl_depth_ini( mythid )

C     !DESCRIPTION: \bv
c     *=================================================================
c     | SUBROUTINE ctrl_depth_ini
c     | Add the depth part of the control vector to the model state
c     | and update the tile halos.
c     | The control vector is defined in the header file "ctrl.h".
c     *=================================================================
C     \ev

C     !USES:
      implicit none

c     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
#include "GRID.h"
#include "CTRL_SIZE.h"
#include "ctrl.h"
#include "ctrl_dummy.h"
#include "optim.h"
#ifdef ALLOW_GENARR2D_CONTROL
# include "CTRL_GENARR.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
c     == routine arguments ==
      integer mythid

#ifdef ALLOW_DEPTH_CONTROL
C     !LOCAL VARIABLES:
c     == local variables ==

#ifndef ALLOW_GENARR2D_CONTROL

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer il

      logical equal
      logical doglobalread
      logical ladinit

      character*( 80)   fnamedepth
      character*(max_len_mbuf) msgbuf

      _RL     fac

c     == external ==
      integer  ilnblnk
      external ilnblnk

#else

      integer igen_depth, iarr

#endif


c     == end of interface ==
CEOP

#ifdef ALLOW_GENARR2D_CONTROL

Cts ---
C ctrl_map_genarr2d takes in field: r_low, and adds update to it:
C xx_r_low. This is different from old setup: xx_r_low includes
C contribution from r_low and ctrl file.

      ! Need to set igen index for ctrl_map_genarr2d
      igen_depth=0
      do iarr = 1, maxCtrlArr2D
        if (xx_genarr2d_weight(iarr).NE.' ') then
          if (xx_genarr2d_file(iarr)(1:8).EQ.'xx_depth') then
            igen_depth=iarr
          endif
        endif   
      enddo

      ! r_low = r_low + xx_depth.*.data
      if (igen_depth.gt.0)
     & call ctrl_map_genarr2d(r_low, igen_depth, myThid)

cts ---
#else

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

      doglobalread = .false.
      ladinit      = .false.

      equal = .true.

      if ( equal ) then
        fac = 1. _d 0
      else
        fac = 0. _d 0
      endif

      write(standardmessageunit,'(21x,a)')
     &     'ctrl_depth_ini: ctrl update R_low,'
      write(standardmessageunit,'(21x,a)')
     &     '                adding the control vector.'

C     Re-initialize hFacC, so that TAMC/TAF can see it
C     Once hFacC is the control variable, and not its anomaly
C     this will be no longer necessary
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do j = 1-oly,sny+oly
         do i = 1-olx,snx+olx
          tmpfld2d(i,j,bi,bj) = 0. _d 0
         enddo
        enddo
       enddo
      enddo

C--

      do bj = jtlo,jthi
        do bi = itlo,ithi
          do j = 1-oly,sny+oly
            do i = 1-olx,snx+olx
              xx_r_low(i,j,bi,bj) = 0. _d 0
            enddo
          enddo
        enddo
      enddo

      il=ilnblnk( xx_depth_file )
      write(fnamedepth(1:80),'(2a,i10.10)')
     &     xx_depth_file(1:il),'.',optimcycle
      call active_read_xy( fnamedepth, tmpfld2d, 1,
     &                     doglobalread, ladinit, optimcycle,
     &                     mythid, xx_depth_dummy )
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do j = jmin,jmax
         do i = imin,imax
          xx_r_low(i,j,bi,bj) = r_low(i,j,bi,bj)
     &                        + fac*tmpfld2d(i,j,bi,bj)
         enddo
        enddo
       enddo
      enddo

c--   Update the tile edges.

      _EXCH_XY_RL( xx_r_low, myThid )

#endif /* ALLOW_GENARR2D_CONTROL */

#endif /* ALLOW_DEPTH_CONTROL */

      return
      end

