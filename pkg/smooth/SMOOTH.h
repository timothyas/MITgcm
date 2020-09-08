C pkg/smooth constants

      integer     smoothprec
      parameter ( smoothprec = 32 )

      logical smooth3DdoImpldiff
      parameter ( smooth3DdoImpldiff = .TRUE. )

      integer smoothOpNbMax
      PARAMETER ( smoothOpNbMax = 10 )

      _RL smooth2DdelTime,smooth3DdelTime
      PARAMETER ( smooth2DdelTime = 1. _d 0 )
      PARAMETER ( smooth3DdelTime = 1. _d 0 )

c fields:

      _RS
     & smooth_recip_hFacC(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smooth_hFacW(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smooth_hFacS(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy)

C --- 3D Fields
      integer smooth3Dnbt(smoothOpNbMax),
     & smooth3DtypeZ(smoothOpNbMax),smooth3DsizeZ(smoothOpNbMax),
     & smooth3DtypeH(smoothOpNbMax),smooth3DsizeH(smoothOpNbMax),
     & smooth3Dfilter(smoothOpNbMax),smooth3DNbRand(smoothOpNbMax),
     & smooth3DJacobiMaxIters(smoothOpNbMax)

      _RL smooth3DtotTime,
     & smooth3D_Lx0(smoothOpNbMax),
     & smooth3D_Ly0(smoothOpNbMax), smooth3D_Lz0(smoothOpNbMax),
     & smooth3D_Lx(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smooth3D_Ly(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smooth3D_Lz(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smooth3Dnorm (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy)

      _RL 
     & smooth3D_kappaR (1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smooth3D_Kwx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kwy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kwz(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kvy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kuz(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kvz(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kuy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3D_Kvx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3DDelta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth3DRandNorm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      character*(5) smooth3DmaskName(smoothOpNbMax)
      character*(11) smooth3DAlgorithm(smoothOpNbMax)

C --- 2D Fields
      integer smooth2Dnbt(smoothOpNbMax),
     & smooth2Dtype(smoothOpNbMax),smooth2Dsize(smoothOpNbMax),
     & smooth2Dfilter(smoothOpNbMax),smooth2DNbRand(smoothOpNbMax),
     & smooth2DJacobiMaxIters(smoothOpNbMax)

      _RL smooth2DtotTime,
     & smooth2D_Lx0(smoothOpNbMax),smooth2D_Ly0(smoothOpNbMax),
     & smooth2D_Lx(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy),
     & smooth2D_Ly(1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy),
     & smooth2Dnorm (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smoothXZnorm (1-OLx:sNx+OLx,Nr,nSx,nSy),
     & smoothYZnorm (1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smooth2DSOROmega(smoothOpNbMax)

      character*(5) smooth2DmaskName(smoothOpNbMax)
      character*(11) smooth2DAlgorithm(smoothOpNbMax)
      character*(3) smooth2DDims(smoothOpNbMax)

      _RL 
     & smooth2D_Kux (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy),
     & smooth2D_Kvy (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy),
     & smooth2DDelta(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smooth2DRandNorm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy),
     & smoothXZ_Kux (1-Olx:sNx+Olx,Nr,nSx,nSy),
     & smoothXZ_Kwz (1-Olx:sNx+Olx,Nr,nSx,nSy),
     & smoothXZDelta(1-OLx:sNx+OLx,Nr,nSx,nSy),
     & smoothXZRandNorm(1-OLx:sNx+OLx,Nr,nSx,nSy),
     & smoothYZ_Kvy (1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smoothYZ_Kwz (1-Oly:sNy+Oly,Nr,nSx,nSy),
     & smoothYZDelta(1-OLy:sNy+OLy,Nr,nSx,nSy),
     & smoothYZRandNorm(1-OLy:sNy+OLy,Nr,nSx,nSy)

C --- Common Blocks
      COMMON /smooth_flds_rs/
     & smooth_recip_hFacC, smooth_hFacW, smooth_hFacS

      COMMON /smooth_flds_rl/
     & smooth3Dnorm,smooth3DtotTime,
     & smooth2Dnorm,smooth2DtotTime,
     & smooth3D_Lx0,smooth3D_Ly0,smooth3D_Lz0,
     & smooth3D_Lx,smooth3D_Ly,smooth3D_Lz,
     & smooth2D_Lx0,smooth2D_Ly0,
     & smooth2D_Lx,smooth2D_Ly,
     & smoothXZNorm, smoothYZNorm,
     & smooth2DSOROmega


      COMMON /smooth_flds_c/
     & smooth3DAlgorithm, smooth2DAlgorithm,
     & smooth3DmaskName, smooth2DmaskName,
     & smooth2DDims

      COMMON /smooth_operators_i/
     & smooth3Dnbt, smooth2Dnbt,
     & smooth2Dtype, smooth2Dsize,
     & smooth3DtypeZ, smooth3DsizeZ,
     & smooth3DtypeH, smooth3DsizeH,
     & smooth2Dfilter, smooth3Dfilter,
     & smooth2DNbRand, smooth3DNbRand,
     & smooth2DJacobiMaxIters, smooth3DJacobiMaxIters

      COMMON /smooth_operators_r/
     & smooth3D_kappaR,smooth3D_Kwx,smooth3D_Kwy,smooth3D_Kwz,
     & smooth3D_Kux,smooth3D_Kvy,smooth3D_Kuz,smooth3D_Kvz,
     & smooth3D_Kuy,smooth3D_Kvx,
     & smooth3DDelta,smooth3DRandNorm,
     & smooth2D_Kux,smooth2D_Kvy,
     & smooth2DDelta,smooth2DRandNorm,
     & smoothXZ_Kux,smoothXZ_Kwz,
     & smoothXZDelta,smoothXZRandNorm,
     & smoothYZ_Kvy,smoothYZ_Kwz,
     & smoothYZDelta,smoothYZRandNorm


 
