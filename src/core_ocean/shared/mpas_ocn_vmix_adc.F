!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_adc
!
!> \brief MPAS ocean vertical mixing interface to new ADC closure 
!> \author Luke Van Roekel
!> \date   10 April 2020
!> \details
!>  This module contains the routines for calls into the new ADC library
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_adc

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log
   use mpas_vector_reconstruction
   use mpas_timekeeping

   use ocn_adc_mixing
   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes
   use ocn_turbulence
   use ocn_constants
   use ocn_config

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vmix_adc_mix, &
             ocn_vmix_adc_init, &
             ocn_vmix_adc_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: i1, i2
   logical :: adcOn
   real (kind=RKIND) :: backgroundVisc, backgroundDiff
   integer, parameter :: num_water_types = 5
   real (kind=RKIND), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_RKIND, 0.62_RKIND, 0.67_RKIND, 0.77_RKIND, 0.78_RKIND /), &
      depth1 = (/ 0.35_RKIND, 0.60_RKIND, 1.00_RKIND, 1.50_RKIND, 1.40_RKIND /), &
      depth2 = (/ 23.0_RKIND, 20.0_RKIND, 17.0_RKIND, 14.0_RKIND, 7.90_RKIND /)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_adc_mix
!
!> \brief   Computes mixing tendency using the assumed distribution closure
!> \author  Luke Van Roekel
!> \date    10 April 2020
!> \details
!>  This routine computes the vertical mixing tendencies for momentum
!>  and tracers by calling ADC mixing library
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_adc_mix(meshPool, statePool, forcingPool, diagnosticsPool, adcMixingPool, dt, &
       err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in), optional :: timeLevelIn !< Input: time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool         !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool   !< Input/Output: diagnostic information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

      type (mpas_pool_type), intent(inout) :: &
         adcMixingPool !< INput/Output: turbulence arrays

      real(KIND=RKIND), intent(in) :: dt
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer :: &
        maxLevelCell, nEdgesOnCell

      real (kind=RKIND), dimension(:), pointer :: &
        stressX, stressY, stressZ, stressZonal, stressMeridional, fCell, &
        surfaceStress, surfaceStressMagnitude, bottomDepth,  &
        penetrativeTemperatureFlux, angleEdge, surfaceThicknessFlux, &
        surfaceThicknessFluxRunoff, dThreshMLD

      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThickness, velocityZonal, velocityMeridional, &
        tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
        tracerGroupSurfaceFluxRemoved, BruntVaisalaFreqTop,  &
        inSituThermalExpansionCoeff, inSituSalinityContractionCoeff

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, pointer :: indexTemperature, indexSalinity

      integer :: k, i, iCell, jCell, iNeighbor, iter, timeLevel, kIndexOBL, kav, iEdge, nCells
      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray

     real (kind=RKIND), dimension(:), allocatable :: zonalAverage, meridionalAverage, uw_sfc, &
                           vw_sfc, latVal, ws_sfc, wt_sfc

      real (kind=RKIND) :: fac, dep1, dep2, fracAbsorbed, fracAbsorbedRunoff,  &
                           etime

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the
      !   mixing fields have to be added/merged together
      !
      !-----------------------------------------------------------------

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if ( .not. adcOn ) return

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      !
      ! set pointers for fields related to position on sphere
      !
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFlux", tracerGroupSurfaceFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRunoff", tracerGroupSurfaceFluxRunoff)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRemoved", tracerGroupSurfaceFluxRemoved)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFlux", surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFluxRunoff", surfaceThicknessFluxRunoff)

      call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff', inSituThermalExpansionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'inSituSalinityContractionCoeff', inSituSalinityContractionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dThreshMLD', dThreshMLD)
      !
      call mpas_pool_get_array(meshPool, 'fCell', fCell)
      !
      ! set pointers for fields related to vertical mesh
      !
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      !
      ! set pointers for fields related ocean state
      !
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop',BruntVaisalaFreqTop)

      ! need to use RBFs to get stress at cell centers
      call mpas_pool_get_array(forcingPool, 'windStressZonal', stressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', stressMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'stressX', stressX)
      call mpas_pool_get_array(diagnosticsPool, 'stressY', stressY)
      call mpas_pool_get_array(diagnosticsPool, 'stressZ', stressZ)

      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)

      allocate(zonalAverage(nVertLevels),meridionalAverage(nVertLevels))

      nCells = nCellsArray( 2 )
      nEdges = nEdgesArray( 1 )

      allocate(uw_sfc(nCells), vw_sfc(nCells), ws_sfc(nCells), wt_sfc(nCells))
      !$omp do schedule(runtime)
      do iEdge=1,nEdges
         surfaceStress(iEdge) = 0.0_RKIND
      enddo
      !$omp end do

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         surfaceThicknessFlux(iCell) = 0.0_RKIND
         surfaceThicknessFluxRunoff(iCell) = 0.0_RKIND
         surfaceStressMagnitude(iCell) = 0.0_RKIND
         tracerGroupSurfaceFluxRunoff(:, iCell) = 0.0_RKIND
         tracerGroupSurfaceFlux(:, iCell) = 0.0_RKIND
         tracerGroupSurfaceFluxRemoved(:, iCell) = 0.0_RKIND
      end do
      !$omp end do

      call ocn_surface_bulk_forcing_vel(meshPool, forcingPool, surfaceStress, surfaceStressMagnitude, err)

      call ocn_surface_land_ice_fluxes_vel(meshPool, diagnosticsPool, surfaceStress, surfaceStressMagnitude, err)

      if(config_use_activeTracers_surface_bulk_forcing) then
        call ocn_surface_bulk_forcing_tracers(meshPool, "activeTracers", forcingPool, activeTracers,  &
                                            tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
                                            tracerGroupSurfaceFluxRemoved, dt, layerThickness, err)
      endif

      call ocn_surface_land_ice_fluxes_tracers(meshPool, "activeTracers", forcingPool, tracerGroupSurfaceFlux, err)

      call ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)

      !need surface thickness flux too

      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      fac = rfac(config_jerlov_water_type)
      dep1 = depth1(config_jerlov_water_type)
      dep2 = depth2(config_jerlov_water_type)

      !$omp do schedule(runtime)
      do iCell = 1, nCells

         fracAbsorbed = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient, -100.0_RKIND) )
         fracAbsorbedRunoff = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient_runoff, &
                              -100.0_RKIND) )

         wt_sfc(iCell) = tracerGroupSurfaceFlux(indexTemperature, iCell) - fracAbsorbed * surfaceThicknessFlux(iCell) * &
                        activeTracers(indexTemperature,1,iCell) + tracerGroupSurfaceFluxRunoff(indexTemperature,iCell) * &
                        fracAbsorbedRunoff

         ws_sfc(iCell) = tracerGroupSurfaceFlux(indexSalinity, iCell) - fracAbsorbed * surfaceThicknessFlux(iCell) *  &
                        activeTracers(indexSalinity,1,iCell)

         uw_sfc(iCell) = stressZonal(iCell) / rho_sw

         if(abs(stressMeridional(iCell)) < 1e-6_RKIND) then
                 stressMeridional(iCell) = 0.0_RKIND
         endif

         vw_sfc(iCell) = stressMeridional(iCell) / rho_sw

         adcMixing % ze(1,iCell) = 0.0_RKIND
         adcMixing % zm(1,iCell) = -layerThickness(1,iCell)/2.0_RKIND
         do k=2,maxLevelCell(iCell)
            adcMixing % ze(k,iCell) = adcMixing % ze(k-1,iCell) - layerThickness(k-1,iCell)
            adcMixing % zm(k,iCell) = adcMixing % ze(k, iCell) - layerThickness(k,iCell)/2.0_RKIND
         enddo
         k = maxLevelCell(iCell)+1
         adcMixing % ze(k,iCell) = adcMixing % ze(k-1,iCell) - layerThickness(k-1,iCell)
         do k = maxLevelCell(iCell) + 1, nVertLevels
            adcMixing % ze(k+1,iCell) = adcMixing % ze(maxLevelCell(iCell)+1,iCell)
            adcMixing % zm(k,iCell) = adcMixing % ze(maxLevelCell(iCell)+1,iCell)
         enddo

     enddo
     !$omp end do
     !call the individual ADC function and pass the array members maybe?  probably not
     !this will update T in place, so verify it sends time level n+1
     call ADC_main_loop(nCells, nVertLevels, dt, activeTracers, velocityZonal, velocityMeridional, &
          BruntVaisalaFreqTop, uw_sfc, vw_sfc, wt_sfc, ws_sfc, inSituThermalExpansionCoeff,        &
          inSituSalinityContractionCoeff, fCell)

         !--------------------------------------------------------------------

   end subroutine ocn_vmix_adc_mix!}}}

!***********************************************************************
!
!  routine ocn_vmix_adc_init
!
!> \brief   Initializes ocean vertical mixing quantities by using
!> \ get and puts into ADC library
!> \author  Luke Van Roekel
!> \date    10 April 2020
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean. Parameters are set by calling into ADC library
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_adc_init(domain,err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (mpas_pool_type), pointer :: meshPool
      type (block_type), pointer :: block
      type (MPAS_TimeInterval_type) :: timeStep
      integer, intent(out) :: err !< Output: error flag

      integer, pointer :: nVertLevels
      integer :: k, iCell, nCells
      integer, dimension(:), pointer :: nCellsArray

      real(KIND=RKIND) :: dt, numiterations
      real, dimension(:,:), pointer :: KspsU, KspsD, w2New, sigma, tumd, sumd,  wumd, uwNew, &
            vwNew, u2New, v2New, uvNew, utNew, vtNew, vsNew, usNew, length, w2Cur, uwCur,    &
            vwCur, u2Cur, v2Cur, uvCur, utCur, vtCur, vsCur, usCur

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      adcOn = config_adc_enable

      !
      ! only initialize if Assumed Distribution Closure is turned on
      !
      if (.not.config_adc_enable) return

      if(config_use_cvmix) then
        call mpas_log_write( "Cannot use ADC mixing and CVMix at the same time:" // &
               "set either config_use_cvmix or config_adc_enable to .false.", MPAS_LOG_CRIT)
      endif
      !
      ! pull nVertLevels out of the mesh structure
      !

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      if ( mod(dt, config_adc_timestep) .ne. 0 ) then
         numiterations = ceiling(dt, config_adc_timestep)
         config_adc_timestep = dt / numiterations
          call mpas_log_write( "WARNING: timestep for ADC scheme was not an even multiple of" \\
               "config_dt, config_adc_timestep has been reset appropriately", MPAS_LOG_WARN)
      end if

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      call ocn_turbulenceCreate(domain)

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
         nCells = nCellsArray( size(nCellsArray) )
         do iCell = 1, nCells
            do k = 1, nVertLevels
               adcMixing%KspsU(k,iCell) = adcMixing%epsilon
               adcMixing%KspsD(k,iCell) = adcMixing%epsilon
               adcMixing%w2(:,k,iCell) = 0.0_RKIND
               adcMixing%areaFraction(k,iCell) = 0.5_RKIND
               adcMixing%tumd(k,iCell) = 0.0_RKIND
               adcMixing%sumd(k,iCell) = 0.0_RKIND
               adcMixing%wumd(k,iCell) = 0.0_RKIND
               adcMixing%uw(:,k,iCell) = 0.0_RKIND
               adcMixing%vw(:,k,iCell) = 0.0_RKIND
               adcMixing%u2(:,k,iCell) = 0.0_RKIND
               adcMixing%v2(:,k,iCell) = 0.0_RKIND
               adcMixing%uv(:,k,iCell) = 0.0_RKIND
               adcMixing%ut(:,k,iCell) = 0.0_RKIND
               adcMixing%vt(:,k,iCell) = 0.0_RKIND
               adcMixing%wt(:,k,iCell) = 0.0_RKIND
               adcMixing%us(:,k,iCell) = 0.0_RKIND
               adcMixing%vs(:,k,iCell) = 0.0_RKIND
               adcMixing%ws(:,k,iCell) = 0.0_RKIND
            enddo
         enddo
         block => block % next
      enddo
   end subroutine ocn_vmix_adc_init!}}}

!***********************************************************************

   subroutine ocn_vmix_adc_finalize(ierr)

     integer, intent(out) :: ierr

     deallocate(adcMixing%KspsU, adcMixing%KspsD, adcMixing%ze, adcMixing%zm,          &
                adcMixing%eps, adcMixing%length, adcMixing%lenspsD, adcMixing%lenspsU, &
                adcMixing%KhU, adcMixing%KhD, adcMixing%KmU, adcMixing%KmD,            &
                adcMixing%wt_spsU, adcMixing%wt_spsD, adcMixing%ws_spsU, adcMixing%ws_spsD, &
                adcMixing%uw2, adcMixing%vw2, adcMixing%u2w, adcMixing%v2w, adcMixing%w2t,  &
                adcMixing%w2s, adcMixing%wts, adcMixing%uvw, adcMixing%uwt, adcMixing%vwt,  &
                adcMixing%ws2, adcMixing%wt2, adcMixing%areaFraction, adcMixing%Entrainment, &
                adcMixing%Detrainment, adcMixing%tumd, adcMixing%sumd, adcMixing%wumd,       &
                adcMixing%w2tend1, adcMixing%w2tend2, adcMixing%w2tend3, adcMixing%w2tend4,  &
                adcMixing%w2tend5, adcMixing%wttend1, adcMixing%wttend2, adcMixing%wttend3,  &
                adcMixing%wttend4, adcMixing%wttend5, adcMixing%wstend1, adcMixing%wstend2,  &
                adcMixing%wstend3, adcMixing%wstend4, adcMixing%wstend5, adcMixing%uwtend1,  &
                adcMixing%uwtend2, adcMixing%uwtend3, adcMixing%uwtend4, adcMixing%uwtend5,  &
                adcMixing%vwtend1, adcMixing%vwtend2, adcMixing%vwtend3, adcMixing%vwtend4,  &
                adcMixing%vwtend5, adcMixing%u2tend1, adcMixing%u2tend2, adcMixing%u2tend3,  &
                adcMixing%u2tend4, adcMixing%u2tend5, adcMixing%v2tend1, adcMixing%v2tend2,  &
                adcMixing%v2tend3, adcMixing%v2tend4, adcMixing%v2tend5, adcMixing%u2cliptend, &
                adcMixing%v2cliptend, adcMixing%w2cliptend, adcMixing%u2, adcMixing%v2,      &
                adcMixing%w2, adcMixing%wt, adcMixing%ws, adcMixing%uw, adcMixing%vw,        &
                adcMixing%w3, adcMixing%uv, adcMixing%t2, adcMixing%s2, adcMixing%ut,        &
                adcMixing%vt, adcMixing%us, adcMixing%vs, adcMixing%ts)

   end subroutine ocn_vmix_adc_finalize

end module ocn_vmix_adc

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
