#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_SHELFICE
# include "SHELFICE_OPTIONS.h"
#endif


      subroutine smooth_filtervar2D (smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_filtervar2D
C     | o Routine that computes the filter variance
C     |   field associated with a diffusion operator, as part
C     |   a 2D spatial correlation operator (smooth_correld2D.F)
C     |   See Weaver and Courtier 01 for details.
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SMOOTH.h"
#if (defined ALLOW_SHI2D_CONTROL) || (defined ALLOW_SHIFWFLX_CONTROL)
# include "SHELFICE.h"
#endif

      integer smoothOpNb, myThid

      Real*8   port_rand, port_rand_norm
      EXTERNAL PORT_RAND, PORT_RAND_NORM

      Real*8 rand_seed
      integer timeVals(8)

      integer i,j,k, bi, bj, ii, jj, kk
      integer itlo,ithi,jtlo,jthi
      integer diLoc,djLoc,dkLoc
      integer nbt_in
      character*( 80) fnamegeneric
      logical doNormalize
      _RL smoothTmpFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smoothTmpVar (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smoothTmpMean(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smooth2Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      doNormalize=.FALSE.

c if smooth2Dfilter(smoothOpNb)=0: the filter variance field
c has been computed earlier and is already in the run directory
c so this routine does not do anything

      IF (smooth2Dfilter(smoothOpNb).NE.0) then

       nbt_in=smooth2Dnbt(smoothOpNb)/2

c --- allow a different mask other than maskC
c     Note: this is essentially a copy of ctrl_get_mask
c     but is repeated here for package independence
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
          if (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskC') then
           smooth2Dmask(i,j,bi,bj) = maskC(i,j,1,bi,bj)
          elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskW') then
             smooth2Dmask(i,j,bi,bj) = maskW(i,j,1,bi,bj)
          elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskS') then
             smooth2Dmask(i,j,bi,bj) = maskS(i,j,1,bi,bj)
#if (defined ALLOW_SHI2D_CONTROL) || (defined ALLOW_SHIFWFLX_CONTROL)
          elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskI') then
             smooth2Dmask(i,j,bi,bj) = maskSHI(i,j,1,bi,bj)
#endif
          endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Read smoothing [i.e. diffusion] operator:
      CALL smooth_read_operator_2d(smoothOpNb,myThid)

      IF (smooth2Dfilter(smoothOpNb).EQ.2) then
c compute the normalization matrix using the approximate method
c
c This method can be quite expensive -- so that the approximate
c method (see below) is usually the prefered one.
c The exact method can be used to check the accuracy
c of the approximate method results (that can be predicted).
c
c note: the exact method requires the adjoint of smooth_diff2D.F (see below)

C --- Check which algorithm is used
      if (smooth2DAlgorithm(smoothOpNb)(1:6).eq.'matern') then
        WRITE(errorMessageUnit,'(A,/,A)' )
     &   "Exact method for matern prior covariance not "
     &   "implemented yet... but would be pretty expensive."
        CALL ALL_PROC_DIE( myThid )
        STOP 'ABNORMAL END: S/R smooth_filtervar2D'

       else

        diLoc=15 !int(5*smooth_L/smooth_dx)
        djLoc=20 !int(5*smooth_L/smooth_dx)

        DO ii=1,diLoc
        DO jj=1,djLoc

        DO bj=jtlo,jthi
         DO bi=itlo,ithi
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             smoothTmpFld(i,j,bi,bj)=0.
            ENDDO
           ENDDO

           DO j=jj,sNy,djLoc
            DO i=ii,sNx,diLoc
             smoothTmpFld(i,j,bi,bj)=1.
            ENDDO
           ENDDO
         ENDDO
        ENDDO

c note: as we go to adjoint part, we need to have 0 in overlaps
c       so we must NOT have done an exchange for smoothTmpFld

c adjoint:
       WRITE(errorMessageUnit,'(A,/,A)' )
     &  "you need to have adsmooth_diff2D compiled and then:",
     &  "uncomment the line below and comment the stop"
       CALL ALL_PROC_DIE( myThid )
       STOP 'ABNORMAL END: S/R smooth_filtervar2D'
c       call adsmooth_diff2D(smoothTmpFld,smooth2dmask,nbt_in,mythid)

c division by sqrt(area)*sqrt(area) [1 to end adj, 1 to begin fwd]
       DO bj = jtlo,jthi
        DO bi = itlo,ithi
         DO j = 1,sNy
          DO i = 1,sNx
c division by ~volume:
       smoothTmpFld(i,j,bi,bj)=smoothTmpFld(i,j,bi,bj)
     &  *recip_rA(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

c coming out of adjoint part: overlaps are 0
c going in fwd part: we need to fill them up
       CALL EXCH_XY_RL ( smoothTmpFld,myThid )

c fwd:
       CALL smooth_diff2D(smoothTmpFld,smooth2dmask,nbt_in,mythid)

c convert variance to normalization factor:
       DO bj=jtlo,jthi
        DO bi=itlo,ithi
          DO j=jj,sNy,djLoc
           DO i=ii,sNx,diLoc
            if (smooth2dmask(i,j,bi,bj).NE.0) then
               smooth2Dnorm(i,j,bi,bj)=
     &         1/sqrt(smoothTmpFld(i,j,bi,bj))
            endif
           ENDDO
          ENDDO
        ENDDO
       ENDDO

       ENDDO      !DO ii=1,diLoc
       ENDDO      !DO jj=1,djLoc

      endif ! end algorithm check

      ELSEIF (smooth2Dfilter(smoothOpNb).EQ.1) then
c compute the normalization matrix using the approximate method

C --- 1. Initialize mean and variance
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           smoothTmpMean(i,j,bi,bj) = 0. _d 0
           smoothTmpVar(i,j,bi,bj)  = 0. _d 0
          ENDDO
         ENDDO
       ENDDO
      ENDDO

c --- 2. Initialize random number generator
c     seed is 1023 * pid * milliseconds from date/time
      rand_seed = 1023. _d 0
      call date_and_time(VALUES=timeVals)
      rand_seed = rand_seed * REAL((1+myPid)*timeVals(8),8)
      smoothTmpFld(1,1,1,1)=port_rand(rand_seed)

      DO ii=1,smooth2DNbRand(smoothOpNb)
c --- 3. set points wehre mask!=0 to random number
c --- WC01: divide by sqrt(area)
       DO bj=jtlo,jthi
        DO bi=itlo,ithi
         DO j=1,sNy
          DO i=1,sNx

           smoothTmpFld(i,j,bi,bj)   = 0. _d 0
           if (smooth2dmask(i,j,bi,bj).NE.0) then

             smoothTmpFld(i,j,bi,bj)=port_rand_norm(-rand_seed)

             if (smooth2DAlgorithm(smoothOpNb)(1:4).eq.'WC01') then
               smoothTmpFld(i,j,bi,bj)=
     &         smoothTmpFld(i,j,bi,bj)*sqrt(recip_rA(i,j,bi,bj))
             elseif (smooth2DAlgorithm(smoothOpNb)(1:6).eq.'matern') then
               smoothTmpFld(i,j,bi,bj)=smoothTmpFld(i,j,bi,bj)*
     &         smooth2DRandNorm(i,j,bi,bj)
            endif
           endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       CALL EXCH_XY_RL ( smoothTmpFld, myThid )

c --- 4. Smooth random field with appropriate algorithm
       if (smooth2DAlgorithm(smoothOpNb)(1:4).eq.'WC01')
     &  call smooth_diff2D(smoothTmpFld,smooth2dmask,nbt_in,mythid)

       if (smooth2DAlgorithm(smoothOpNb)(1:6).eq.'matern')
     &  call smooth_inverse_matern_2d(
     &          smoothTmpFld,smooth2DMask,smoothOpNb,doNormalize,myThid)

c --- 5. accumulate statistics (to compute the variance later)
       DO bj=jtlo,jthi
        DO bi=itlo,ithi
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           smoothTmpVar(i,j,bi,bj)=smoothTmpVar(i,j,bi,bj)
     &      +smoothTmpFld(i,j,bi,bj)*smoothTmpFld(i,j,bi,bj)/
     &       smooth2DNbRand(smoothOpNb)
           smoothTmpMean(i,j,bi,bj)=smoothTmpMean(i,j,bi,bj)
     &      +smoothTmpFld(i,j,bi,bj)/smooth2DNbRand(smoothOpNb)
          ENDDO
         ENDDO
        ENDDO
       ENDDO

C --- temporary, write out all samples
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dfld',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      1, smoothTmpFld, ii, 1, myThid )

       WRITE(standardMessageUnit,'(A,I4,A,I4)')
     & 'smooth_filtervar2D: ',ii,' members done out of ',
     &  smooth2DNbRand(smoothOpNb)

      ENDDO ! end ii->nbRand

c --- 6. compute variance and convert it to normalization factor:
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          if (smooth2dmask(i,j,bi,bj).NE.0) then
           smooth2Dnorm(i,j,bi,bj)=
     &      1/sqrt ( smooth2DNbRand(smoothOpNb)/
     &      (smooth2DNbRand(smoothOpNb)-1)*
     &      (smoothTmpVar(i,j,bi,bj) -
     &          smoothTmpMean(i,j,bi,bj)*smoothTmpMean(i,j,bi,bj)) )
           endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO

      ENDIF ! end smooth2Dfilter == 1 or 2

c write smooth2Dnorm to file:
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dnorm',smoothOpNb
      CALL WRITE_REC_3D_RL(fnamegeneric,smoothprec,
     &            1,smooth2Dnorm,1,1,mythid)
      CALL EXCH_XY_RL ( smooth2Dnorm,  myThid )

c write out mean too
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dmean',smoothOpNb
      CALL WRITE_REC_3D_RL(fnamegeneric,smoothprec,
     &            1,smoothTmpMean,1,1,mythid)

      ENDIF ! end smooth2Dtype != 0

      END

