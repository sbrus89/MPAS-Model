!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! \file mpas_ocn_turbulence.F
!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_turbulence
!
!>  \brief MPAS ocean turbulence arrays for ADC mixing structure
!!  \author Luke Van Roekel
!!  \date   14 Jan 2020
!!  \details
!!  This module creates and maintains a primary ocean mesh structure
!!  for the ADC mixing scheme.  Eventually it will send arrays over to GPU
!!  but for now just loads on to cpu
!
!-------------------------------------------------------------------------------

module ocn_turbulence

   ! module dependencies
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   use ocn_config

   implicit none
   private

   !----------------------------------------------------------------------------
   !
   ! Public parameters
   !
   !----------------------------------------------------------------------------

   !type, public :: adcMixingType

   real(kind=RKIND),public :: epsilon, &! small value below which w2, v2, u2 cannot go
                       sigmat, Ko, gamma1, beta5, c1, c2, c8, c10, c11, &
                       alpha1, alpha2, alpha3, B1, Kt, grav, c_1, c_2, c_mom, &
                       c_therm, c_mom_w3, c_pt, c_ps, c_pv, &
                       kappa_FL, kappa_w3, kappa_VAR, Cww_D, Cww_E

   real(kind=RKIND),public, dimension(:,:), allocatable :: KspsU, KspsD, ze, zm, &
                       eps, lendn, lenup, length, lenspsD, lenspsU, KhU, KhD, KmU,  &
                       KmD, wt_spsU, wt_spsD, ws_spsU, ws_spsD, uw2,  &
                       vw2, u2w, v2w, w2t, w2s, wts, uvw, uwt, vwt,   &
                       uws, vws, ws2, wt2, areaFraction, Entrainment,        &
                       Detrainment, tumd, sumd, wumd, Mc

   real(kind=RKIND), public, dimension(:,:), allocatable :: w2tend1, w2tend2,     &
                       w2tend3, w2tend4, w2tend5, wttend1, wttend2,   &
                       wttend3, wttend4, wttend5, w3tend1, w3tend2,   &
                       w3tend3, w3tend5, w3tend4, wstend1, wstend2,   &
                       wstend3, wstend4, wstend5, uwtend1, uwtend2,   &
                       uwtend3, uwtend4, uwtend5, vwtend1, vwtend2,   &
                       vwtend3, vwtend4, vwtend5, u2tend1, u2tend2,   &
                       u2tend3, u2tend4, u2tend5, v2tend1, v2tend2,   &
                       v2tend3, v2tend4, v2tend5, u2cliptend,         &
                       v2cliptend, w2cliptend

   real(kind=RKIND), public, dimension(:,:,:), allocatable :: u2, v2, w2, t2, s2, &
                       uw, vw, wt, ws, w3, uv, ut, vt, us, vs, ts

!   end type adcMixingType

   !ntype (adcMixingType), public :: adcMixing

   !----------------------------------------------------------------------------
   !
   ! Public member functions
   !
   !----------------------------------------------------------------------------

   public ::                      &
      ocn_turbulenceCreate,       &
      ocn_turbulenceUpdateFields, &
      ocn_turbulenceDestroy

   contains

   subroutine ocn_turbulenceCreate(domain)

   type(domain_type) :: &
         domain                    !< [in] MPAS type to describe domain

   integer :: &
         blockCount               ! counter for number of blocks

   type(block_type), pointer :: &
         block                    ! variables in current subblock

   type (mpas_pool_type), pointer :: &
         adcDiagnosticsPool, adcPrognosticsPool, adcTendPool, meshPool  ! mesh variables in MPAS pool structure

   real(kind=RKIND), dimension(:,:), pointer :: KspsUTmp, KspsDTmp, zeTmp, zmTmp, &
                       epsTmp, lengthTmp, lenspsDTmp, lenspsUTmp, KhUTmp, KhDTmp, KmUTmp,  &
                       KmDTmp, wt_spsUTmp, wt_spsDTmp, ws_spsUTmp, ws_spsDTmp, uw2Tmp,  &
                       vw2Tmp, u2wTmp, v2wTmp, w2tTmp, w2sTmp, wtsTmp, uvwTmp, uwtTmp, vwtTMP,   &
                       uwsTmp, vwsTmp, ws2Tmp, wt2Tmp, areaFractionTmp, EntrainmentTmp,        &
                       DetrainmentTmp, tumdTmp, sumdTmp, wumdTmp, McTmp, lendnTmp, lenupTmp

   real(kind=RKIND), dimension(:,:), pointer :: w2tend1Tmp, w2tend2Tmp,     &
                       w2tend3Tmp, w2tend4Tmp, w2tend5Tmp, wttend1Tmp, wttend2Tmp,   &
                       wttend3Tmp, wttend4Tmp, wttend5Tmp, w3tend1Tmp, w3tend2Tmp,   &
                       w3tend3Tmp, w3tend5Tmp, w3tend4Tmp, wstend1Tmp, wstend2Tmp,   &
                       wstend3Tmp, wstend4Tmp, wstend5Tmp, uwtend1Tmp, uwtend2Tmp,   &
                       uwtend3Tmp, uwtend4Tmp, uwtend5Tmp, vwtend1Tmp, vwtend2Tmp,   &
                       vwtend3Tmp, vwtend4Tmp, vwtend5Tmp, u2tend1Tmp, u2tend2Tmp,   &
                       u2tend3Tmp, u2tend4Tmp, u2tend5Tmp, v2tend1Tmp, v2tend2Tmp,   &
                       v2tend3Tmp, v2tend4Tmp, v2tend5Tmp, u2cliptendTmp,         &
                       v2cliptendTmp, w2cliptendTmp

   real(kind=RKIND), dimension(:,:), pointer :: u2Tmp, v2Tmp, w2Tmp, t2Tmp, s2Tmp, &
                       uwTmp, vwTmp, wtTmp, wsTmp, w3Tmp, uvTmp, utTmp, vtTmp, usTmp,&
                       vsTmp, tsTmp

   integer, pointer :: nVertLevels, nVertLevelsP1
   integer :: nCells
   integer, dimension(:), pointer :: nCellsArray

   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      nCells = nCellsArray( size(nCellsArray) )
      call mpas_pool_get_subpool(block % structs, 'adcDiagnosticArrays', adcDiagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'adcTendArrays', adcTendPool)
      call mpas_pool_get_subpool(block % structs, 'adcPrognosticArrays', adcPrognosticsPool)
      sigmat = config_adc_sigmat
      epsilon = config_adc_epsilon
      Ko = config_adc_Ko
      gamma1 = config_adc_gamma1
      beta5 = config_adc_beta5
      c1 = 2.88_RKIND
      c2 = 3.84_RKIND
      c8 = 3.0_RKIND
      c10 = 5.0_RKIND
      c11 = 0.1_RKIND
      alpha1 = 0.984_RKIND
      alpha2 = 0.568_RKIND
      alpha3 = 0.29_RKIND
      B1 = 16.6_RKIND
      Kt = 0.4_RKIND
      c_1 = config_adc_c1
      c_2 = config_adc_c2
      c_mom = config_adc_Cmom
      c_therm = config_adc_Ctherm
      c_mom_w3 = config_adc_Cmom_w3
      c_pt = config_adc_cpt
      c_ps = config_adc_cps
      c_pv = config_adc_cpv
      Cww_D = config_adc_CwwD
      Cww_E = config_adc_CwwE
      kappa_FL = config_adc_kappaFL
      kappa_VAR = config_adc_kappaVAR
      kappa_w3 = config_adc_kappaW3

      call mpas_pool_get_array(adcDiagnosticsPool, 'KspsD', KspsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KspsU', KspsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ze', zeTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'zm', zmTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'eps', epsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'length', lengthTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenup', lenupTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lendn', lendnTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenspsD', lenspsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'lenspsU', lenspsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KhU', KhUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KhD', KhDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KmU', KmUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'KmD', KmDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wt_spsU', wt_spsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wt_spsD', wt_spsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ws_spsU', ws_spsUTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ws_spsD', ws_spsDTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uw2', uw2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'vw2', vw2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'u2w', u2wTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'v2w', v2wTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'w2t', w2tTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'w2s', w2sTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wts', wtsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uvw', uvwTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uwt', uwtTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'vwt', vwtTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'uws', uwsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'vws', vwsTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'ws2', ws2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wt2', wt2Tmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'areaFraction', areaFractionTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Entrainment', EntrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Detrainment', DetrainmentTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'tumd', tumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'sumd', sumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'wumd', wumdTmp)
      call mpas_pool_get_array(adcDiagnosticsPool, 'Mc', McTmp)

      allocate( KspsD( nVertLevelsP1, nCells ), KspsU( nVertLevelsP1, nCells) )
      allocate( ze( nVertLevelsP1, nCells ), zm( nVertLevels, nCells ))
      allocate( eps( nVertLevelsP1, nCells ), length( nVertLevelsP1, nCells ))
      allocate( lendn( nVertLevelsP1, nCells ), lenup( nVertLevelsP1, nCells ))
      allocate( lenspsU( nVertLevelsP1, nCells), lenspsD( nVertLevelsP1, nCells))
      allocate( KhU( nVertLevelsP1, nCells), KhD( nVertLevelsP1, nCells))
      allocate( KmU( nVertLevelsP1, nCells), KmD( nVertLevelsP1, nCells))
      allocate( wt_spsU( nVertLevelsP1, nCells), wt_spsD( nVertLevelsP1, nCells))
      allocate( ws_spsU( nVertLevelsP1, nCells), ws_spsD( nVertLevelsP1, nCells))
      allocate( uw2( nVertLevels, nCells ), vw2( nVertLevels, nCells))
      allocate( u2w( nVertLevels, nCells ), v2w( nVertLevels, nCells ))
      allocate( w2t( nVertLevels, nCells ), w2s( nVertLevels, nCells ))
      allocate( wts( nVertLevels, nCells ), uvw( nVertLevels, nCells ))
      allocate( uwt( nVertLevels, nCells ), vwt( nVertLevels, nCells ))
      allocate( uws( nVertLevels, nCells ), vws( nVertLevels, nCells ))
      allocate( ws2( nVertLevels, nCells ), wt2( nVertLevels, nCells ))
      allocate( areaFraction( nVertLevelsP1, nCells ))
      allocate( Entrainment( nVertLevelsP1, nCells ), Detrainment( nVertLevelsP1, nCells ))
      allocate( tumd( nVertLevelsP1, nCells ), sumd( nVertLevelsP1, nCells ))
      allocate( wumd( nVertLevelsP1, nCells ), Mc( nVertLevelsP1, nCells ))

      KspsD = KspsDTmp
      KspsU = KspsUTmp
      ze = zeTmp
      zm = zmTmp
      eps = epsTmp
      length = lengthTmp
      lenup = lenupTmp
      lendn = lendnTmp
      lenspsD = lenspsDTmp
      lenspsU = lenspsUTmp
      KhU = KhUTmp
      KhD = KhDTmp
      KmU = KmUTmp
      KmD = KmDTmp
      wt_spsU = wt_spsUTmp
      wt_spsD = wt_spsDTmp
      ws_spsU = ws_spsUTmp
      ws_spsD = ws_spsDTmp
      uw2 = uw2Tmp
      vw2 = vw2Tmp
      u2w = u2wTmp
      v2w = v2wTmp
      w2t = w2tTmp
      w2s = w2sTmp
      wts = wtsTmp
      uvw = uvwTmp
      uwt = uwtTmp
      vwt = vwtTmp
      uws = uwsTmp
      vws = vwsTmp
      ws2 = ws2Tmp
      wt2 = wt2Tmp
      areaFraction = areaFractionTmp
      Entrainment = EntrainmentTmp
      Detrainment = DetrainmentTmp
      tumd = tumdTmp
      sumd = sumdTmp
      wumd = wumdTmp

      call mpas_pool_get_array(adcTendPool, 'w2tend1', w2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend2', w2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend3', w2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend4', w2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w2tend5', w2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend1', w3tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend2', w3tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend3', w3tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend4', w3tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'w3tend5', w3tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend1', wttend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend2', wttend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend3', wttend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend4', wttend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wttend5', wttend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend1', wstend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend2', wstend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend3', wstend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend4', wstend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'wstend5', wstend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend1', uwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend2', uwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend3', uwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend4', uwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'uwtend5', uwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend1', vwtend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend2', vwtend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend3', vwtend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend4', vwtend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'vwtend5', vwtend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend1', u2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend2', u2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend3', u2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend4', u2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2tend5', u2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend1', v2tend1Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend2', v2tend2Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend3', v2tend3Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend4', v2tend4Tmp)
      call mpas_pool_get_array(adcTendPool, 'v2tend5', v2tend5Tmp)
      call mpas_pool_get_array(adcTendPool, 'u2cliptend', u2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'v2cliptend', v2cliptendTmp)
      call mpas_pool_get_array(adcTendPool, 'w2cliptend', w2cliptendTmp)

      allocate( w2tend1( nVertLevelsP1, nCells ))
      allocate( w2tend2( nVertLevelsP1, nCells ))
      allocate( w2tend3( nVertLevelsP1, nCells ))
      allocate( w2tend4( nVertLevelsP1, nCells ))
      allocate( w2tend5( nVertLevelsP1, nCells ))
      allocate( wttend1( nVertLevelsP1, nCells ))
      allocate( wttend2( nVertLevelsP1, nCells ))
      allocate( wttend3( nVertLevelsP1, nCells ))
      allocate( wttend4( nVertLevelsP1, nCells ))
      allocate( wttend5( nVertLevelsP1, nCells ))
      allocate( wstend1( nVertLevelsP1, nCells ))
      allocate( wstend2( nVertLevelsP1, nCells ))
      allocate( wstend3( nVertLevelsP1, nCells ))
      allocate( wstend4( nVertLevelsP1, nCells ))
      allocate( wstend5( nVertLevelsP1, nCells ))
      allocate( uwtend1( nVertLevelsP1, nCells ))
      allocate( uwtend2( nVertLevelsP1, nCells ))
      allocate( uwtend3( nVertLevelsP1, nCells ))
      allocate( uwtend4( nVertLevelsP1, nCells ))
      allocate( uwtend5( nVertLevelsP1, nCells ))
      allocate( vwtend1( nVertLevelsP1, nCells ))
      allocate( vwtend2( nVertLevelsP1, nCells ))
      allocate( vwtend3( nVertLevelsP1, nCells ))
      allocate( vwtend4( nVertLevelsP1, nCells ))
      allocate( vwtend5( nVertLevelsP1, nCells ))
      allocate( v2tend1( nVertLevelsP1, nCells ))
      allocate( v2tend2( nVertLevelsP1, nCells ))
      allocate( v2tend3( nVertLevelsP1, nCells ))
      allocate( v2tend4( nVertLevelsP1, nCells ))
      allocate( v2tend5( nVertLevelsP1, nCells ))
      allocate( u2tend1( nVertLevelsP1, nCells ))
      allocate( u2tend2( nVertLevelsP1, nCells ))
      allocate( u2tend3( nVertLevelsP1, nCells ))
      allocate( u2tend4( nVertLevelsP1, nCells ))
      allocate( u2tend5( nVertLevelsP1, nCells ))
      allocate( w3tend1( nVertLevels, nCells ))
      allocate( w3tend2( nVertLevels, nCells ))
      allocate( w3tend3( nVertLevels, nCells ))
      allocate( w3tend4( nVertLevels, nCells ))
      allocate( w3tend5( nVertLevels, nCells ))
      allocate( u2cliptend( nVertLevelsP1, nCells ))
      allocate( v2cliptend( nVertLevelsP1, nCells ))
      allocate( w2cliptend( nVertLevelsP1, nCells ))

      w2tend1 = w2tend1Tmp
      w2tend2 = w2tend2Tmp
      w2tend3 = w2tend3Tmp
      w2tend4 = w2tend4Tmp
      w2tend5 = w2tend5Tmp
      w3tend1 = w3tend1Tmp
      w3tend2 = w3tend2Tmp
      w3tend3 = w3tend3Tmp
      w3tend4 = w3tend4Tmp
      w3tend5 = w3tend5Tmp
      wttend1 = wttend1Tmp
      wttend2 = wttend2Tmp
      wttend3 = wttend3Tmp
      wttend4 = wttend4Tmp
      wttend5 = wttend5Tmp
      wstend1 = wstend1Tmp
      wstend2 = wstend2Tmp
      wstend3 = wstend3Tmp
      wstend4 = wstend4Tmp
      wstend5 = wstend5Tmp
      uwtend1 = uwtend1Tmp
      uwtend2 = uwtend2Tmp
      uwtend3 = uwtend3Tmp
      uwtend4 = uwtend4Tmp
      uwtend5 = uwtend5Tmp
      vwtend1 = vwtend1Tmp
      vwtend2 = vwtend2Tmp
      vwtend3 = vwtend3Tmp
      vwtend4 = vwtend4Tmp
      vwtend5 = vwtend5Tmp
      u2tend1 = u2tend1Tmp
      u2tend2 = u2tend2Tmp
      u2tend3 = u2tend3Tmp
      u2tend4 = u2tend4Tmp
      u2tend5 = u2tend5Tmp
      v2tend1 = v2tend1Tmp
      v2tend2 = v2tend2Tmp
      v2tend3 = v2tend3Tmp
      v2tend4 = v2tend4Tmp
      v2tend5 = v2tend5Tmp
      u2cliptend = u2cliptendTmp
      v2cliptend = v2cliptendTmp
      w2cliptend = w2cliptendTmp

      call mpas_pool_get_array(adcPrognosticsPool, 'u2', u2Tmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'v2', v2Tmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'w2', w2Tmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'wt', wtTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'ws', wsTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'uw', uwTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'vw', vwTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'w3', w3Tmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'uv', uvTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 't2', t2Tmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'ut', utTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'vt', vtTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'us', usTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'vs', vsTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 'ts', tsTmp, 1)
      call mpas_pool_get_array(adcPrognosticsPool, 's2', s2Tmp, 1)

      allocate( u2( 2, nVertLevelsP1, nCells ) )
      allocate( v2( 2, nVertLevelsP1, nCells ) )
      allocate( w2( 2, nVertLevelsP1, nCells ) )
      allocate( wt( 2, nVertLevelsP1, nCells ) )
      allocate( ws( 2, nVertLevelsP1, nCells ) )
      allocate( uw( 2, nVertLevelsP1, nCells ) )
      allocate( vw( 2, nVertLevelsP1, nCells ) )
      allocate( w3( 2, nVertLevels, nCells ) )
      allocate( uv( 2, nVertLevelsP1, nCells ) )
      allocate( t2( 2, nVertLevelsP1, nCells ) )
      allocate( ut( 2, nVertLevelsP1, nCells ) )
      allocate( vt( 2, nVertLevelsP1, nCells ) )
      allocate( us( 2, nVertLevelsP1, nCells ) )
      allocate( vs( 2, nVertLevelsP1, nCells ) )
      allocate( ts( 2, nVertLevelsP1, nCells ) )
      allocate( s2( 2, nVertLevelsP1, nCells ) )

      u2(1,:,:) = u2Tmp(:,:nCells)
      v2(1,:,:) = v2Tmp(:,:nCells)
      w2(1,:,:) = w2Tmp(:,:nCells)
      wt(1,:,:) = wtTmp(:,:nCells)
      ws(1,:,:) = wsTmp(:,:nCells)
      uw(1,:,:) = uwTmp(:,:nCells)
      vw(1,:,:) = vwTmp(:,:nCells)
      w3(1,:,:) = w3Tmp(:,:nCells)
      uv(1,:,:) = uvTmp(:,:nCells)
      t2(1,:,:) = t2Tmp(:,:nCells)
      ut(1,:,:) = utTmp(:,:nCells)
      vt(1,:,:) = vtTmp(:,:nCells)
      us(1,:,:) = usTmp(:,:nCells)
      vs(1,:,:) = vsTmp(:,:nCells)
      ts(1,:,:) = tsTmp(:,:nCells)
      s2(1,:,:) = s2Tmp(:,:nCells)

      call mpas_pool_get_array(adcPrognosticsPool, 'u2', u2Tmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'v2', v2Tmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'w2', w2Tmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'wt', wtTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'ws', wsTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'uw', uwTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'vw', vwTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'w3', w3Tmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'uv', uvTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 't2', t2Tmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'ut', utTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'vt', vtTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'us', usTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'vs', vsTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 'ts', tsTmp, 2)
      call mpas_pool_get_array(adcPrognosticsPool, 's2', s2Tmp, 2)

      u2(2,:,:) = u2Tmp(:,:nCells)
      v2(2,:,:) = v2Tmp(:,:nCells)
      w2(2,:,:) = w2Tmp(:,:nCells)
      wt(2,:,:) = wtTmp(:,:nCells)
      ws(2,:,:) = wsTmp(:,:nCells)
      uw(2,:,:) = uwTmp(:,:nCells)
      vw(2,:,:) = vwTmp(:,:nCells)
      w3(2,:,:) = w3Tmp(:,:nCells)
      uv(2,:,:) = uvTmp(:,:nCells)
      t2(2,:,:) = t2Tmp(:,:nCells)
      ut(2,:,:) = utTmp(:,:nCells)
      vt(2,:,:) = vtTmp(:,:nCells)
      us(2,:,:) = usTmp(:,:nCells)
      vs(2,:,:) = vsTmp(:,:nCells)
      ts(2,:,:) = tsTmp(:,:nCells)
      s2(2,:,:) = s2Tmp(:,:nCells)
      block => block % next
   end do
!   at end could do a copy to accelerator

   end subroutine ocn_turbulenceCreate

   subroutine ocn_turbulenceUpdateFields

!   would only be needed to do device updates

   end subroutine ocn_turbulenceUpdateFields

   subroutine ocn_turbulenceDestroy(ierr)

     integer, intent(inout) :: ierr
!   deallocate stuff, set integers and constants to zero, nullify stuff on gpu

      sigmat = 0.0_RKIND
      epsilon = 0.0_RKIND
      Ko = 0.0_RKIND
      gamma1 = 0.0_RKIND
      beta5 = 0.0_RKIND
      c1 = 0.0_RKIND
      c2 = 0.0_RKIND
      c8 = 0.0_RKIND
      c10 = 0.0_RKIND
      c11 = 0.0_RKIND
      alpha1 = 0.0_RKIND
      alpha2 = 0.0_RKIND
      alpha3 = 0.0_RKIND
      B1 = 0.0_RKIND
      Kt = 0.0_RKIND
      c_1 = 0.0_RKIND
      c_2 = 0.0_RKIND
      c_mom = 0.0_RKIND
      c_therm = 0.0_RKIND
      c_mom_w3 = 0.0_RKIND
      c_pt = 0.0_RKIND
      c_ps = 0.0_RKIND
      c_pv = 0.0_RKIND
      Cww_D = 0.0_RKIND
      Cww_E = 0.0_RKIND
      kappa_FL = 0.0_RKIND
      kappa_VAR = 0.0_RKIND
      kappa_w3 = 0.0_RKIND

      deallocate(KspsD,    &
                 KspsU,    &
                 ze,       &
                 zm,       &
                 eps,      &
                 length,   &
                 lenup,    &
                 lendn,    &
                 lenspsD,  &
                 lenspsU,  &
                 KhU,      &
                 KhD,      &
                 KmU,      &
                 KmD,      &
                 wt_spsU,  &
                 wt_spsD,  &
                 ws_spsU,  &
                 ws_spsD,  &
                 uw2,      &
                 vw2,      &
                 u2w,      &
                 v2w,      &
                 w2t,      &
                 w2s,      &
                 wts,      &
                 uvw,      &
                 uwt,      &
                 vwt,      &
                 uws,      &
                 vws,      &
                 ws2,      &
                 wt2,      &
                 areaFraction,    &
                 Entrainment, &
                 Detrainment, &
                 tumd,        &
                 sumd,        &
                 wumd,        &
                 w2tend1,     &
                 w2tend2,     &
                 w2tend3,     &
                 w2tend4,     &
                 w2tend5,     &
                 w3tend1,     &
                 w3tend2,     &
                 w3tend3,     &
                 w3tend4,     &
                 w3tend5,     &
                 wttend1,     &
                 wttend2,     &
                 wttend3,     &
                 wttend4,     &
                 wttend5,     &
                 wstend1,     &
                 wstend2,     &
                 wstend3,     &
                 wstend4,     &
                 wstend5,     &
                 uwtend1,     &
                 uwtend2,     &
                 uwtend3,     &
                 uwtend4,     &
                 uwtend5,     &
                 vwtend1,     &
                 vwtend2,     &
                 vwtend3,     &
                 vwtend4,     &
                 vwtend5,     &
                 u2tend1,     &
                 u2tend2,     &
                 u2tend3,     &
                 u2tend4,     &
                 u2tend5,     &
                 v2tend1,     &
                 v2tend2,     &
                 v2tend3,     &
                 v2tend4,     &
                 v2tend5,     &
                 u2cliptend,  &
                 v2cliptend,  &
                 w2cliptend,  &
                 u2,          &
                 v2,          &
                 w2,          &
                 wt,          &
                 ws,          &
                 uw,          &
                 vw,          &
                 w3,          &
                 uv,          &
                 t2,          &
                 ut,          &
                 vt,          &
                 us,          &
                 vs,          &
                 ts,          &
                 s2)

   end subroutine ocn_turbulenceDestroy

   end module ocn_turbulence
