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

   subroutine ocn_vmix_adc_mix(meshPool, statePool, forcingPool, diagnosticsPool, adcMixingPool, &
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
        latCell, lonCell, bottomDepth, surfaceBuoyancyForcing, surfaceFrictionVelocity, fCell, &
        boundaryLayerDepth, ssh, indexBoundaryLayerDepth, dcEdge, dvEdge, areaCell, iceFraction, &
        boundaryLayerDepthSmooth, windSpeed10m

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        zMid, zTop, density, displacedDensity, potentialDensity, &
        bulkRichardsonNumber, RiTopOfCell, BruntVaisalaFreqTop, &
        bulkRichardsonNumberBuoy, bulkRichardsonNumberShear, unresolvedShear, normalVelocity

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
      integer, pointer :: index_vertNonLocalFluxTemp, config_cvmix_num_ri_smooth_loops
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnCell, cellMask

      integer :: k, i, iCell, jCell, iNeighbor, iter, timeLevel, kIndexOBL, kav, iEdge, nCells
      integer :: edgeCount, nEdges, topIndex, nsmooth, kpp_stage
      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), allocatable :: surfaceAverageIndex

      real (kind=RKIND) :: r, layerSum, bulkRichardsonNumberStop, sfc_layer_depth, invAreaCell
      real (kind=RKIND) :: normalVelocityAv, factor, delU2, areaSum, blTemp
      real (kind=RKIND) :: sigma, turbulentScalarVelocityScalePoint
      real (kind=RKIND), dimension(:), allocatable :: Nsqr_iface, turbulentScalarVelocityScale, &
                                                      deltaVelocitySquared, normalVelocitySum, &
                                                      potentialDensitySum, RiTemp
      real (kind=RKIND), dimension(:), allocatable, target :: RiSmoothed, BVFSmoothed, OBLDepths, interfaceForcings
      logical :: bulkRichardsonFlag


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
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      !
      ! set pointers for fields related to position on sphere
      !
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'fCell', fCell)

      !
      ! set pointers for fields related to vertical mesh
      !
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      !
      ! set pointers for fields related ocean state
      !
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop',BruntVaisalaFreqTop)

      !
      ! set pointers for fields related forcing at ocean surface
      !
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)

      nCells = nCellsArray( size(nCellsArray) )

      nCells = nCellsArray( 3 )
      do iCell = 1,nCells
               ! fill vertical position of column
         ! CVMix assume top of ocean is at z=0, so building all z-coordinate data based on layerThickness
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
      integer, intent(out) :: err !< Output: error flag

      integer, pointer :: nVertLevels
      integer :: k, iCell, nCells
      integer, dimension(:), pointer :: nCellsArray

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

      !
      ! pull nVertLevels out of the mesh structure
      !
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
