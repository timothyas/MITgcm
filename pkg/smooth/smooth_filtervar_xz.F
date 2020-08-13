#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif 

      subroutine smooth_filtervar_xz (smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_filtervar2D_xz
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
#ifdef ALLOW_OBCS
# include "OBCS_GRID.h"
#endif

      integer smoothOpNb, myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_NORTH) || (defined ALLOW_OBCS_SOUTH)

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
      _RL smoothTmpFld (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL smoothTmpVar (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL smoothTmpMean(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL smooth2Dmask (1-OLx:sNx+OLx,Nr,nSx,nSy)

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

C ---   Northside
        IF ((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &      (tileHasOBN(bi,bj))) then
         DO i = 1-OLx,sNx+OLx
          j = OB_Jn(i,bi,bj)
          DO k = 1,Nr
           if (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskC') then
            smooth2Dmask(i,k,bi,bj) = maskC(i,j,k,bi,bj)
           elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskW') then
            smooth2Dmask(i,k,bi,bj) = maskW(i,j,k,bi,bj)
           elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskS') then
            smooth2Dmask(i,k,bi,bj) = maskS(i,j,k,bi,bj)
           endif
          ENDDO
         ENDDO

C ---   Southside
        ELSEIF ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &          (tileHasOBS(bi,bj))) then
         DO i = 1-OLx,sNx+OLx
          j = OB_Js(i,bi,bj)
          DO k = 1,Nr
           if (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskC') then
            smooth2Dmask(i,k,bi,bj) = maskC(i,j,k,bi,bj)
           elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskW') then
            smooth2Dmask(i,k,bi,bj) = maskW(i,j,k,bi,bj)
           elseif (smooth2DmaskName(smoothOpNb)(1:5).eq.'maskS') then
            smooth2Dmask(i,k,bi,bj) = maskS(i,j,k,bi,bj)
           endif
          ENDDO
         ENDDO
        ELSE
         DO k=1,Nr
          DO i=1-OLx,sNx+OLx
           smooth2DMask(i,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

C --- Read smoothing [i.e. diffusion] operator:
      CALL smooth_read_operator_xz(smoothOpNb,myThid)

      IF (smooth2Dfilter(smoothOpNb).EQ.2) then

        WRITE(errorMessageUnit,'(A,/,A)' )
     &   "Exact method for filter variance of xz slice not implemented"
        CALL ALL_PROC_DIE( myThid )
        STOP 'ABNORMAL END: S/R smooth_filtervar2D_xz'

      ELSEIF (smooth2Dfilter(smoothOpNb).EQ.1) then
c compute the normalization matrix using the approximate method

C --- 1. Initialize mean and variance
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k = 1,Nr
         DO i = 1-OLx,sNx+OLx
          smoothTmpMean(i,k,bi,bj) = 0. _d 0
          smoothTmpVar (i,k,bi,bj)  = 0. _d 0
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
       DO bj=jtlo,jthi
        DO bi=itlo,ithi

         IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &        (tileHasOBN(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &        (tileHasOBS(bi,bj)))) then
          DO k = 1,Nr
           DO i = 1-OLx,sNx+OLx
            if (smooth2DMask(i,k,bi,bj).NE.0) then

             smoothTmpFld(i,k,bi,bj)=port_rand_norm(-rand_seed)

             if (smooth2DAlgorithm(smoothOpNb).eq.'matern') then
               smoothTmpFld(i,k,bi,bj)=smoothTmpFld(i,k,bi,bj)*
     &         smoothXZRandNorm(i,k,bi,bj)
             endif
            else
             smoothTmpFld(i,k,bi,bj)   = 0. _d 0
            endif
           ENDDO
          ENDDO
         ELSE
          DO k = 1,Nr
           DO i = 1-OLx,sNx+OLx
            smoothTmpFld(i,k,bi,bj)   = 0. _d 0
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO
       CALL SMOOTH_EXCH_XZ_RL ( smoothTmpFld, smoothOpNb, myThid )

c --- 4. Smooth random field with appropriate algorithm
       if (smooth2DAlgorithm(smoothOpNb).eq.'matern')
     &  call smooth_inverse_matern_xz(
     &          smoothTmpFld,smooth2DMask,smoothOpNb,doNormalize,myThid)

c --- 5. accumulate statistics (to compute the variance later)
       DO bj=jtlo,jthi
        DO bi=itlo,ithi

         IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &        (tileHasOBN(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &        (tileHasOBS(bi,bj)))) then
          DO k = 1,Nr
           DO i = 1-OLx,sNx+OLx
            smoothTmpVar(i,k,bi,bj)=smoothTmpVar(i,k,bi,bj)
     &       +smoothTmpFld(i,k,bi,bj)*smoothTmpFld(i,k,bi,bj)/
     &        smooth2DNbRand(smoothOpNb)
            smoothTmpMean(i,k,bi,bj)=smoothTmpMean(i,k,bi,bj)
     &       +smoothTmpFld(i,k,bi,bj)/smooth2DNbRand(smoothOpNb)
           ENDDO
          ENDDO
         ENDIF
        ENDDO
       ENDDO

       WRITE(standardMessageUnit,'(A,I4,A,I4)')
     & 'smooth_filtervar_xz: ',ii,' members done out of ',
     &  smooth2DNbRand(smoothOpNb)

      ENDDO ! end ii->nbRand

c --- 6. compute variance and convert it to normalization factor:
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &       (tileHasOBN(bi,bj))).OR.
     &      ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &       (tileHasOBS(bi,bj)))) then
         DO k = 1,Nr
          DO i = 1-OLx,sNx+OLx
           if (smooth2dmask(i,k,bi,bj).NE.0) then
            smoothXZnorm(i,k,bi,bj)=
     &       1/sqrt ( smooth2DNbRand(smoothOpNb)/
     &       (smooth2DNbRand(smoothOpNb)-1)*
     &       (smoothTmpVar(i,k,bi,bj) -
     &           smoothTmpMean(i,k,bi,bj)*smoothTmpMean(i,k,bi,bj)) )
           endif
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

      ENDIF ! end smooth2Dfilter == 1 or 2

c write smooth2Dnorm to file:
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dnorm',smoothOpNb
      CALL WRITE_REC_XZ_RL(fnamegeneric,smoothprec,
     &            Nr,smoothXZnorm,1,1,mythid)
      CALL SMOOTH_EXCH_XZ_RL ( smoothXZnorm,  smoothOpNb, myThid )

c write out mean too
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dmean',smoothOpNb
      CALL WRITE_REC_XZ_RL(fnamegeneric,smoothprec,
     &            Nr,smoothTmpMean,1,1,mythid)

      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dfld',smoothOpNb
      CALL WRITE_REC_XZ_RL( fnamegeneric, smoothprec,
     &                      Nr, smoothTmpFld, 1, 1, myThid )

      ENDIF ! end smooth2Dtype != 0

# endif /* ALLOW_OBCS_NORTH or ALLOW_OBCS_SOUTH*/
#endif /* ALLOW_OBCS*/

      END

