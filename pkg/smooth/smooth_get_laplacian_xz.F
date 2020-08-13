#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_get_laplacian_xz (
     U     center, east, west, up, down,
     I     smoothOpNb, smooth2DMask, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_get_laplacian_xz
C     | o get the matrix to the discretized 2D 
C     |   laplacian operator: div( kappa * grad( ) )
C     | o only the diagonal elements of the diffusion tensor:
C     |
C     |             Kux Kuz    Kux  0 
C     |    kappa =  Kwx Kwz ->  0  Kwz
C     |
C     |   are implemented currently
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"
#ifdef ALLOW_OBCS
# include "OBCS_GRID.h"
#endif

C --- Inputs
C   center, east, west, east, west, up, down :
C       fields that define the laplacian stencil at each point
C       for a point (i,k)
C           center = (i,k)
C           east   = (i+1, k  ); west   = (i-1, k  )
C           up     = (i  , k+1); down   = (i  , k-1)
      _RL center(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,Nr,nSx,nSy)

      integer smoothOpNb
      _RL smooth2Dmask(1-OLx:sNx+OLx,Nr,nSx,nSy)
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_NORTH) || (defined ALLOW_OBCS_SOUTH)

C --- Local parameters
      _RL xA(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL zA(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL recip_vol(1-OLx:sNx+OLx,Nr,nSx,nSy)
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Initialize output
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO k = 1,Nr
         DO i = 1-OLx,sNx+OLx
          center(i,k,bi,bj)=0. _d 0
          east  (i,k,bi,bj)=0. _d 0
          west  (i,k,bi,bj)=0. _d 0
          up    (i,k,bi,bj)=0. _d 0
          down  (i,k,bi,bj)=0. _d 0
          xA    (i,k,bi,bj)=0. _d 0
          zA    (i,k,bi,bj)=0. _d 0
          recip_vol(i,k,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Get discretized Laplacian
C     Note: these are the (swapped sign of) the diagonal
C           elements in applied to input field in smooth_rhs
      DO bj = jtlo,jthi
       DO bi = itlo,ithi

C ---   Northside
        IF ((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &      (tileHasOBN(bi,bj))) then
         DO i = 1-OLx,sNx+OLx
          j = OB_Jn(i,bi,bj)
          DO k = 1,Nr

            xA(i,k,bi,bj) = _dxG(i,j,bi,bj) * drF(k) *
     &          _hFacS(i,j,k,bi,bj)  * _recip_dyC(i,j,bi,bj) *
     &          smoothXZ_Kux(i,k,bi,bj)

            zA(i,k,bi,bj) = _rA(i,j,bi,bj) * 
     &          _hFacC(i,j,k,bi,bj)  * recip_drC(k) *
     &          smoothXZ_Kwz(i,k,bi,bj)

            recip_vol(i,k,bi,bj) = _recip_hFacC(i,j,k,bi,bj)*
     &          recip_rA(i,j,bi,bj) * recip_drF(k)

          ENDDO
         ENDDO

C ---   Southside
        ELSEIF ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &          (tileHasOBS(bi,bj))) then
         DO i = 1-OLx,sNx+OLx
          j = OB_Js(i,bi,bj)
          DO k = 1,Nr

            xA(i,k,bi,bj) = _dxG(i,j,bi,bj) * drF(k) *
     &          _hFacS(i,j,k,bi,bj)  * _recip_dyC(i,j,bi,bj) *
     &          smoothXZ_Kux(i,k,bi,bj)

            zA(i,k,bi,bj) = _rA(i,j,bi,bj) * 
     &          _hFacC(i,j,k,bi,bj)  * recip_drC(k) *
     &          smoothXZ_Kwz(i,k,bi,bj)

            recip_vol(i,k,bi,bj) = _recip_hFacC(i,j,k,bi,bj)*
     &          recip_rA(i,j,bi,bj) * recip_drF(k)

          ENDDO
         ENDDO
        ENDIF
         

        IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &       (tileHasOBN(bi,bj))).OR.
     &      ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &       (tileHasOBS(bi,bj)))) then
         DO k = 1,Nr
          DO i = 1-OLx+1,sNx+OLx-1

           east(i,k,bi,bj)=east(i,k,bi,bj) + xA(i+1,k,bi,bj)*
     &      smooth2Dmask(i,k,bi,bj)*smooth2Dmask(i+1,k,bi,bj)*
     &      recip_vol(i,k,bi,bj)

           west(i,k,bi,bj)=west(i,k,bi,bj) + xA(i,k,bi,bj)*
     &      smooth2Dmask(i-1,k,bi,bj)*smooth2Dmask(i,k,bi,bj)*
     &      recip_vol(i,k,bi,bj)

           if (k.eq.1) then
            up(i,k,bi,bj) = 0. _d 0
           else
            up(i,k,bi,bj)=up(i,k,bi,bj) + zA(i,k,bi,bj)*
     &       smooth2Dmask(i,k-1,bi,bj)*smooth2Dmask(i,k,bi,bj)*
     &       recip_vol(i,k,bi,bj)
           endif

           if (k.eq.Nr) then
            down(i,k,bi,bj) = 0. _d 0
           else
            down(i,k,bi,bj)=down(i,k,bi,bj) + zA(i,k+1,bi,bj)*
     &       smooth2Dmask(i,k,bi,bj)*smooth2Dmask(i,k+1,bi,bj)*
     &       recip_vol(i,k,bi,bj)
           endif

           center(i,k,bi,bj) = center(i,k,bi,bj) -
     &       east(i,k,bi,bj) - west(i,k,bi,bj) -
     &       up  (i,k,bi,bj) - down(i,k,bi,bj)

          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

      CALL SMOOTH_EXCH_XZ_RL( center, smoothOpNb, myThid )
      CALL SMOOTH_EXCH_XZ_RL( east, smoothOpNb, myThid )
      CALL SMOOTH_EXCH_XZ_RL( west, smoothOpNb, myThid )
      CALL SMOOTH_EXCH_XZ_RL( up, smoothOpNb, myThid )
      CALL SMOOTH_EXCH_XZ_RL( down, smoothOpNb, myThid )

# endif
#endif

      end
