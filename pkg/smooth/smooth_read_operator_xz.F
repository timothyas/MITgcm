#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_read_operator_xz (smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_read_operator_xz
C     | o Read in the tensor Kappa:
C     |
C     |             Kux Kuz
C     |     Kappa = Kwx Kwz
C     |
C     |   and also the field and delta (matern)
C     |
C     | o for matern algortihm, read three files for three
C     |   entries to the Kappa tensor
C     |     smoothXZKux<smoothOpNb>
C     |     smoothXZKwz<smoothOpNb>
C     |
C     |   as well as the file smoothXZDelta<smoothOpNb>
C     |
C     |   which form the operator: A = (delta - div(Kappa grad())
C     |
C     | o for matern, all entries except Kux, Kwz are zero
C     |   and if those two files are not found, they are set
C     |   to 1's. If delta is not found, it is set to zero.
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      integer smoothOpNb
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_NORTH) || (defined ALLOW_OBCS_SOUTH)

C --- Local parameters
      character*( 80) fnamegeneric
      character*(MAX_LEN_MBUF) msgBuf
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      logical fexist

C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)


      if (smooth2DAlgorithm(smoothOpNb).eq.'matern') then

c --- Initialize tensor with default values
        do bj = jtlo,jthi
         do bi = itlo,ithi
          do k = 1,Nr
           do i = 1-OLx,sNx+OLx
            smoothXZ_Kux(i,k,bi,bj)  = 1. _d 0
            smoothXZ_Kwz(i,k,bi,bj)  = 1. _d 0
            smoothXZDelta(i,k,bi,bj) = 0. _d 0
            smoothXZRandNorm(i,k,bi,bj) = 1. _d 0
            smoothXZNorm(i,k,bi,bj) = 0. _d 0
           enddo
          enddo
         enddo
        enddo

C ---   Kux
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth2DKux',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_XZ_RL(fnamegeneric, smoothprec, Nr,
     &          smoothXZ_Kux, 1, 1, myThid )
          CALL SMOOTH_EXCH_XZ_RL( smoothXZ_Kux, smoothOpNb, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_XZ: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kux to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   Kwz
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth2DKwz',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_XZ_RL(fnamegeneric, smoothprec, Nr,
     &          smoothXZ_Kwz, 1, 1, myThid )
          CALL SMOOTH_EXCH_XZ_RL( smoothXZ_Kwz, smoothOpNb, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_XZ: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kwz to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   delta:
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth2DDelta',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        IF (fexist) THEN
          CALL READ_REC_XZ_RL(fnamegeneric,smoothprec,
     &               Nr, smoothXZDelta,1, 1, myThid )
          CALL SMOOTH_EXCH_XZ_RL( smoothXZDelta, smoothOpNb, myThid )
        ELSE
          WRITE(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_XZ: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to zeros'
          CALL PRINT_ERROR(msgBuf, myThid)
        ENDIF

C ---   Look for normalization of white noise process
C       i.e. 1/sqrt(det(F))
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth2DRandNorm',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_XZ_RL(fnamegeneric,smoothprec,
     &               Nr, smoothXZRandNorm,1, 1, myThid )
          CALL SMOOTH_EXCH_XZ_RL( smoothXZRandNorm, smoothOpNb, myThid )
        else
          write(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_XZ: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to ones'
          call print_message(msgBuf, standardMessageUnit,
     &        SQUEEZE_RIGHT, myThid)
        endif

      ENDIF ! Check algorithm 'matern'

c ---  read normalization field [i.e. 1/sqrt(var(filter))]:
c     ctrl_check ensures this file exists or is computed, if needed
      fexist = .FALSE.
      write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &    'smooth2Dnorm',smoothOpNb,'.data'
      inquire(file=fnamegeneric,exist=fexist)
      IF (fexist) THEN
        CALL READ_REC_XZ_RL(fnamegeneric,smoothprec,
     &             Nr, smoothXZnorm,1,1,mythid)
        CALL SMOOTH_EXCH_XZ_RL( smoothXZNorm, smoothOpNb, myThid )
      ENDIF


# endif
#endif

      RETURN
      END
