#include "SMOOTH_OPTIONS.h"

      subroutine smooth_read_operator_3d (smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_read_operator_3d
C     | o Read in the tensor Kappa:
C     |
C     |             Kux Kuy Kuz
C     |     Kappa = Kvx Kvy Kvz
C     |             Kwx Kwy Kwz
C     |
C     |   and also the fields kappaR (WC01), and delta (matern)
C     |
C     | o for WC01, read single file: smooth3Doperator<smoothOpNb>
C     |   with 10 fields in order:
C     |
C     |     Kwx, Kwy, Kwz, Kux, Kvy, Kuz, Kvz, Kuy, Kvx, kappaR
C     |
C     | o for matern algortihm, read three files for three
C     |   entries to the Kappa tensor
C     |     smooth3DKux<smoothOpNb>
C     |     smooth3DKvy<smoothOpNb>
C     |     smooth3DKwz<smoothOpNb>
C     |
C     |   as well as the file smooth3DDelta<smoothOpNb>
C     |
C     |   which form the operator: A = (delta - div(Kappa grad())
C     |
C     | o for matern, all entries except Kux, Kvy, Kwz are zero
C     |   and if those three files are not found, they are set
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


C --- WC01: look for single file, smooth3Doperator
      if (smooth3DAlgorithm(smoothOpNb).eq.'WC01') then
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth3Doperator',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kwx, 1, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kwy, 2, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kwz, 3, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kux, 4, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kvy, 5, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kuz, 6, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kvz, 7, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kuy, 8, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kvx, 9, 1, myThid )
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_kappaR, 10, 1, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_3D: ',
     &      'cannot find file ', fnamegeneric
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_READ_OPERATOR_3D' 
        endif ! look for smooth3Doperator file

      elseif (smooth3DAlgorithm(smoothOpNb).eq.'matern') then

c --- Initialize tensor with default values
        do bj = jtlo,jthi
         do bi = itlo,ithi
          do k = 1,sNy
           do j = 1-OLy,sNy+OLy
            do i = 1-OLx,sNx+OLx
             smooth3D_Kux(i,j,k) = 1. _d 0
             smooth3D_Kvy(i,j,k) = 1. _d 0
             smooth3D_Kwz(i,j,k) = 1. _d 0
             smooth3DDelta(i,j,k) = 0. _d 0
             smooth3DRandNorm(i,j,k) = 1. _d 0

             smooth3D_Kuy(i,j,k) = 0. _d 0
             smooth3D_Kuz(i,j,k) = 0. _d 0
             smooth3D_Kvx(i,j,k) = 0. _d 0
             smooth3D_Kvz(i,j,k) = 0. _d 0
             smooth3D_Kwx(i,j,k) = 0. _d 0
             smooth3D_Kwy(i,j,k) = 0. _d 0
             smooth3D_kappaR(i,j,k) = 0. _d 0
            enddo
           enddo
          enddo
         enddo
        enddo

C ---   Kux
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth3DKux',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kux, 1, 1, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_3D: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kux to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   Kvy
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth3DKvy',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kvy, 1, 1, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_3D: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kvy to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   Kwz
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth3DKwz',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, Nr,
     &          smooth3D_Kwz, 1, 1, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_3D: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kwz to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   delta:
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth3DDelta',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        IF (fexist) THEN
          CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &               Nr, smooth3DDelta,1, 1, myThid )
          CALL EXCH_XYZ_RL ( smooth3DDelta, myThid )
        ELSE
          WRITE(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_3D: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to zeros'
          CALL PRINT_ERROR(msgBuf, myThid)
        ENDIF

C ---   Look for normalization of white noise process
C       i.e. 1/sqrt(det(F))
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth3DRandNorm',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &               Nr, smooth3DRandNorm,1, 1, myThid )
          CALL EXCH_XYZ_RL ( smooth3DRandNorm, myThid )
        else
          write(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_3D: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to ones'
          call print_message(msgBuf, standardMessageUnit,
     &        SQUEEZE_RIGHT, myThid)
        endif

      ENDIF ! Check algorithm 'WC01' or 'matern'

c ---  read normalization field [i.e. 1/sqrt(var(filter))]:
c     ctrl_check ensures this file exists or is computed, if needed
      fexist = .FALSE.
      write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &    'smooth3Dnorm',smoothOpNb,'.data'
      inquire(file=fnamegeneric,exist=fexist)
      IF (fexist) THEN
        CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &             Nr, smooth3Dnorm,1,1,mythid)
        CALL EXCH_XYZ_RL ( smooth3Dnorm, myThid )
      ENDIF

      CALL EXCH_XYZ_RL ( smooth3D_Kwx, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kwy, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kwz, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kux, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kvy, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kuz, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kvz, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kuy, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kvx, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_kappaR, myThid )
      RETURN
      END
