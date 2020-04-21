module ocn_adc_mixing

  use mpas_kind_types
  use ocn_constants
  use ocn_turbulence
  use mpas_log
  use ocn_config

  implicit none

  logical :: defineFirst, stopflag

  integer :: i1,i2

!  integer,parameter :: ntimes = 2

  contains

  subroutine swap_time_levels

    if(i1 == 1) then
      i1 = 2
    else
      i1 = 1
    endif

    if(i2 == 1) then
      i2 = 2
    else
      i2 = 1
    endif

  end subroutine swap_time_levels

  subroutine build_diagnostic_arrays(nCells,nVertLevels,BVF,wtsfc,wssfc, &
        uwsfc, vwsfc, alphaT,betaS, boundaryLayerDepth)
  !construct dTdz, dSdz, dbdz
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, uwsfc, vwsfc
    real,dimension(nCells),intent(inout) :: boundaryLayerDepth
    real,dimension(nVertLevels,nCells),intent(in) :: BVF, alphaT, betaS
    integer :: iCell, k, idx
    real,dimension(nCells) :: wstar
    logical :: first
    real :: maximum, Q

    first = .true.

    do iCell=1,nCells
      maximum = -1.0e-12
      idx = 1

      do k=1,nVertLevels

        if(BVF(k,iCell) > 1.005*maximum .and. first) then
          maximum = BVF(k,iCell)
          idx = k
        elseif(BVF(k,iCell) < maximum) then
          first = .false.
        endif
      enddo

      boundaryLayerDepth(iCell) = abs(adcMixing%ze(idx,iCell))
      print *, boundaryLayerDepth(iCell)
      print *, adcMixing%grav
      print *, 'bs = ',shape(betaS)
      print *, iCell
      print *, alphaT(1,iCell)
      print *, wtsfc(iCell)
      print *, maxval(betaS),minval(betaS)
      print *, wssfc(iCell)
     ! stop
      Q = adcMixing%grav*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))* &
        boundaryLayerDepth(iCell)
      if(Q > 0) then
        wstar(iCell) = abs(Q)**(1.0/3.0)
      else
        wstar(iCell) = 0.0
      endif

      adcMixing%u2(:,1,iCell) = 4.0*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      adcMixing%v2(:,1,iCell) =4.*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      adcMixing%uw(:,1,iCell) = -uwsfc(iCell)
      adcMixing%vw(:,1,iCell) = vwsfc(iCell)
      adcMixing%wt(:,1,iCell) = wtsfc(iCell)
      adcMixing%ws(:,1,iCell) = wssfc(iCell)
    enddo

  end subroutine build_diagnostic_arrays

  subroutine dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(2,nVertLevels,nCells),intent(in) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS
    integer :: i,k, ij

    real,dimension(nVertLevels) :: B, Bup, Bdo
    real,dimension(nVertLevels+1) :: tke, BupEdge, BdoEdge
    real :: sav, tudav, sudav, Tup, Tdo, Sup, Sdo
    real :: s1, z1, zV, sumv, minlen

    !NOTE: will need to convert to some form of displaced density in the mpas framework soon
    !possibly go back to a more traditional length scale formulation
    do i=1,nCells
       tke(:) = 0.5*(adcMixing%u2(i2,:,i) + adcMixing%v2(i2,:,i) + adcMixing%w2(i2,:,i))
       do k=1,nVertLevels
          B(k) = -9.806*(-alphaT(k,i)*(activeTracers(1,k,i) - 15.0_RKIND) + betaS(k,i)*   &
                    (activeTracers(2,k,i) - 35.0_RKIND))

          sav = 0.5_RKIND*(adcMixing%areaFraction(k,i) + adcMixing%areaFraction(k+1,i))
          tudav = 0.5_RKIND*(adcMixing%tumd(k,i) + adcMixing%tumd(k+1,i))
          sudav = 0.5*(adcMixing%sumd(k,i) + adcMixing%sumd(k+1,i))

          Tup = activeTracers(1,k,i) + (1.0_RKIND - sav)*tudav
          Tdo = activeTracers(1,k,i) - sav*tudav
          Sup = activeTracers(2,k,i) + (1.0_RKIND - sav)*sudav
          Sdo = activeTracers(2,k,i) - sav*sudav

          Bup(k) = -9.806*(-alphaT(k,i)*(Tup - 15.0) + betaS(k,i)*(Sup - 35.0))
          Bdo(k) = -9.806*(-alphaT(k,i)*(Tdo - 15.0) + betaS(k,i)*(Sdo - 35.0))

          if(k>1) THEN
             BupEdge(k) = 0.5*(Bup(k-1) + Bup(k))
             BdoEdge(k) = 0.5*(Bdo(k-1) + Bdo(k))
          endif
       enddo

       BdoEdge(nVertLevels+1) = BdoEdge(nVertLevels)
       BupEdge(nVertLevels+1) = BupEdge(nVertLevels)

       BdoEdge(1) = BdoEdge(2)
       BupEdge(1) = BupEdge(2)

       do k=2,nVertLevels

          sumv = 0
          ij=k
          adcMixing%lenup(k,i) = 0.0_RKIND
          do while(sumv <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv + (BupEdge(k) - Bup(ij))*(adcMixing%ze(ij,i)-adcMixing%ze(ij+1,i))**2.0_RKIND/2.0_RKIND
             adcMixing%lenup(k,i) =  adcMixing%lenup(k,i) + abs(adcMixing%ze(ij,i)-adcMixing%ze(ij+1,i))
             ij = ij + 1

             if(sumv > tke(k)) THEN
                ij = ij - 1
                sumv = sumv - (BupEdge(k) - Bup(ij))*(adcMixing%ze(ij,i)-adcMixing%ze(ij+1,i))**2.0/2.0
                adcMixing%lenup(k,i) = adcMixing%lenup(k,i) - abs(adcMixing%ze(ij,i)-adcMixing%ze(ij+1,i))
!                lenup(k,i) = max(0.55,lenup(k,i) + abs((z1-zV)/(s1 - sumv)*(tke(k)-sumv)))
                if(Bup(k-1) - Bup(k) < 0) then
                        minlen = abs(adcMixing%zm(k-1,i) - adcMixing%zm(k,i))
                else
                        minlen = 0.5_RKIND
                endif
                adcMixing%lenup(k,i) = max(minlen, adcMixing%lenup(k,i) + sqrt(2.0_RKIND/(BupEdge(k) -         &
                                Bup(ij-1))*(tke(k) - sumv)))
                exit
             endif

         end do

        !find lendown
        sumv = 0
        ij=k
        adcMixing%lendn(k,i) = 0.0_RKIND
        do while(sumv <= tke(k) .and. ij>1)
           sumv = sumv - (BdoEdge(k) - Bdo(ij-1))*(adcMixing%ze(ij-1,i)-adcMixing%ze(ij,i))**2.0_RKIND/2.0_RKIND
           adcMixing%lendn(k,i) = adcMixing%lendn(k,i) + abs(adcMixing%ze(ij-1,i)-adcMixing%ze(ij,i))
           ij = ij - 1

           if(sumv > tke(k)) THEN
              ij = ij + 1
!              s1 = sumv
!              z1 = zedge(ij,i)
!              zV = zedge(ij-1,i)
              sumv = sumv + (BdoEdge(k) - Bdo(ij-1))*(adcMixing%ze(ij-1,i)-adcMixing%ze(ij,i))**2.0_RKIND/2.0_RKIND
              adcMixing%lendn(k,i) = adcMixing%lendn(k,i) - abs(adcMixing%ze(ij-1,i)-adcMixing%ze(ij,i))
!              lendown(k) = max(0.55,lendown(k) + abs(-(z1-zV)/(sumv)*(tke(k))))
                if(Bdo(k-1) - Bdo(k) < 0) then
                        minlen = abs(adcMixing%zm(k-1,i) - adcMixing%zm(k,i))
                else
                        minlen = 0.55_RKIND
                endif
               adcMixing%lendn(k,i) = max(minlen,adcMixing%lendn(k,i) + sqrt(2.0/(BdoEdge(k) -  &
                                Bdo(ij-1))*(tke(k) - sumv)))
              exit
           endif
        enddo

        adcMixing%length(k,i) = (2.0_RKIND*adcMixing%lenup(k,i)*adcMixing%lendn(k,i))  &
          / (adcMixing%lenup(k,i) + adcMixing%lendn(k,i))
      enddo
   enddo

   adcMixing%length(1,:) = 0.55_RKIND
   adcMixing%length(nVertLevels+1,:) = 0.55_RKIND

  end subroutine dissipation_lengths2

  subroutine build_dissipation_lengths(nCells,nVertLevels,BVF)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nVertLevels+1,nCells),intent(inout) :: BVF

    real :: l, len1, len2, lenmax, KE, integrandTop, integrandBot

    integer :: k,iCell

    do iCell=1,nCells

      integrandTop = 0.0
      integrandBot = 0.0
      do k=1,nvertLevels
        KE = sqrt(0.5_RKIND*(adcMixing%u2(i2,k,iCell) + adcMixing%v2(i2,k,iCell) + adcMixing%w2(i2,k,iCell)))
        integrandTop = integrandTop + 0.5_RKIND*KE*abs(adcMixing%ze(k,iCell)**2.0_RKIND - &
          adcMixing%ze(k+1,iCell)**2.0_RKIND)
        integrandBot = integrandBot + KE*(adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell))
      enddo

      do k=2,nVertLevels
        KE = 0.5_RKIND*(adcMixing%u2(i2,k,iCell) + adcMixing%v2(i2,k,iCell) + adcMixing%w2(i2,k,iCell))

        if(KE > adcMixing%EPSILON) then
          lenmax = 0.53_RKIND*sqrt(2.0_RKIND*KE / (1.0E-15_RKIND + BVF(k,iCell)))
        else
          lenmax = 1.0e6_RKIND
        endif

        len1 = 0.4_RKIND*abs(adcMixing%ze(k,iCell))
        len2 = 0.2_RKIND*integrandTop / (integrandBot + 1.0E-10_RKIND)

        adcMixing%length(k,iCell) = min(1.0_RKIND/(1.0_RKIND / len1 + 1.0_RKIND / len2),lenmax)
        adcMixing%length(k,iCell) = max(adcMixing%length(k,iCell),0.55_RKIND)
      enddo
      adcMixing%length(nVertLevels+1,iCell) = 1e-15_RKIND
!      len(1,iCell) = len(2,iCell) !1e-15
    enddo
  end subroutine build_dissipation_lengths

  subroutine build_areaFraction_updraft_properties(nCells,nVertLevels)
  !builds the updraft area function

  integer,intent(in) :: nCells,nVertLevels
  integer :: iCell, k
  real :: Sw, w3av, lareaFraction, wtav, wsav

  do iCell = 1,nCells
    adcMixing%tumd(1,iCell) = 0.0_RKIND
    adcMixing%wumd(1,iCell) = 0.0_RKIND
    adcMixing%areaFraction(1,iCell) = 0.5_RKIND
    adcMixing%Mc(1,iCell) = 0.0_RKIND
    do k=2,nVertLevels
      w3av = 0.5_RKIND*(adcMixing%w3(i2,k-1,iCell) + adcMixing%w3(i2,k,iCell))

      Sw = w3av / (max(adcMixing%w2(i2,k,iCell)**1.5_RKIND,adcMixing%epsilon))
      lareaFraction = 0.5_RKIND - 0.5_RKIND*Sw / sqrt(4.0_RKIND + Sw**2.0_RKIND)

      if(lareaFraction < 0.01_RKIND) lareaFraction = 0.01_RKIND
      if(lareaFraction > 0.99_RKIND) lareaFraction = 0.99_RKIND

      adcMixing%areaFraction(k,iCell) = lareaFraction
      adcMixing%wumd(k,iCell) = sqrt(adcMixing%w2(i2,k,iCell) / (adcMixing%areaFraction(k,iCell) * &
        (1.0_RKIND - adcMixing%areaFraction(k,iCell))))
      adcMixing%Mc(k,iCell) = adcMixing%areaFraction(k,iCell)*(1.0_RKIND - &
        adcMixing%areaFraction(k,iCell)) * adcMixing%wumd(k,iCell)
    enddo
  enddo

  end subroutine build_areaFraction_updraft_properties

  subroutine calc_scalar_updraft_properties(nCells,nVertLevels,wtsfc, wssfc, alphaT,  &
    betaS, boundaryLayerDepth, tlev)

    integer,intent(in) :: nCells, nVertLevels, tlev
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, boundaryLayerDepth
    real, dimension(nVertLevels,nCells) :: alphaT, betaS

    real :: wtav, McAv, sigav, tumdav, wumdav, sumdav, wb, bld, wstar
    integer :: iCell,k

    do iCell=1,nCells
      do k=2,nVertLevels

        adcMixing%tumd(k,iCell) = adcMixing%wt(tlev,k,iCell) / (1.0E-12_RKIND + adcMixing%Mc(k,iCell))
        adcMixing%sumd(k,iCell) = adcMixing%ws(tlev,k,iCell) / (1.0E-12_RKIND + adcMixing%Mc(k,iCell))
      enddo

      wb = adcMixing%grav*(alphaT(k,iCell)*wtsfc(iCell) - betaS(k,iCell)*wssfc(iCell))
!      wstar = (abs(0.4*boundaryLayerDepth(iCell)*wb))**(1./3.)

      if(wb > 0.0) then
        wb = adcMixing%grav*(alphaT(k,iCell)*wtsfc(iCell) - betaS(k,iCell)*wssfc(iCell))
        wstar = (abs(0.4_RKIND*boundaryLayerDepth(iCell)*wb))**(1.0_RKIND/3.0_RKIND)
        adcMixing%w2t(1,iCell) = -0.3_RKIND*wstar * wtsfc(iCell)
        !Below FIXME!
        adcMixing%w2s(1,iCell) = 0.3_RKIND*wstar * wssfc(iCell)
      else
        adcMixing%w2t(1,iCell) = 0.0_RKIND
        adcMixing%w2s(1,iCell) = 0.0_RKIND
      endif

      !try new boundary condition derived from PDF
      sigav = 0.5_RKIND*(adcMixing%areaFraction(1,iCell) + adcMixing%areaFraction(2,iCell))
      wtav = 0.5_RKIND*(adcMixing%wt(tlev,1,iCell) + adcMixing%wt(tlev,2,iCell))
      McAv = 0.5_RKIND*(adcMixing%w2(tlev,1,iCell) + adcMixing%w2(tlev,2,iCell))
      adcMixing%w2t(1,iCell) = (1.0_RKIND - 2.0_RKIND*sigav)*wtav*sqrt(McAv) / (adcMixing%EPSILON + &
        sigav*(1.0_RKIND-sigav))

      do k=2,nVertLevels
        sigav = 0.5_RKIND*(adcMixing%areaFraction(k,iCell) + adcMixing%areaFraction(k+1,iCell))
        tumdav = 0.5_RKIND*(adcMixing%tumd(k,iCell) + adcMixing%tumd(k+1,iCell))
        sumdav = 0.5_RKIND*(adcMixing%sumd(k,iCell) + adcMixing%sumd(k+1,iCell))
        wumdav = 0.5_RKIND*(adcMixing%wumd(k,iCell) + adcMixing%wumd(k+1,iCell))
        adcMixing%w2t(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*tumdav
        adcMixing%w2s(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*sumdav
      enddo

    enddo
  end subroutine calc_scalar_updraft_properties

  subroutine calc_subplume_fluxes(nCells,nVertLevels,activeTracers,uvel,vvel,BVF,   &
    alphaT,betaS,dt)
  ! builds the subplume tendency terms

  integer,intent(in) :: nCells, nVertLevels
  real,dimension(nVertLevels,nCells),intent(in) :: alphaT,betaS
  real,intent(in) :: dt
  real,dimension(nVertLevels,nCells),intent(in) :: uvel,vvel
  real,dimension(2,nVertLevels,nCells),intent(in) :: activeTracers
  real,dimension(nVertLevels+1,nCells),intent(in) :: BVF

  real,dimension(nVertLevels+1,nCells) :: KspsUtend, KspsDtend
  real :: Uz, Vz, Tz, Sz, B, areaFractionAv,dz,integrandTop,integrandBot, Cval
  !calculate length

  integer :: iCell, k

  do iCell = 1,nCells

    adcMixing%lenspsU(1,iCell) = 0.0_RKIND
    adcMixing%lenspsD(1,iCell) = 0.0_RKIND
    adcMixing%KmU(1,iCell) = 0.0_RKIND
    adcMixing%KhU(1,iCell) = 0.0_RKIND
    adcMixing%KmD(1,iCell) = 0.0_RKIND
    adcMixing%KhD(1,iCell) = 0.0_RKIND
    adcMixing%Entrainment(1,iCell) = 0.0_RKIND
    adcMixing%Detrainment(1,iCell) = 0.0_RKIND

    do k = 2,nVertLevels
      dz = (adcMixing%zm(k-1,iCell) - adcMixing%zm(k,iCell))
      !need to add length scales for Up and Down
      Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dz
      Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dz
      Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dz
      Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dz

      if(BVF(k,iCell) <= 0.0_RKIND) then
        adcMixing%lenspsU(k,iCell) = dz
        adcMixing%lenspsD(k,iCell) = dz
      else
        adcMixing%lenspsU(k,iCell) = min(dz,0.76_RKIND*sqrt(adcMixing%KspsU(k,iCell)/BVF(k,iCell)))
        adcMixing%lenspsD(k,iCell) = min(dz,0.76_RKIND*sqrt(adcMixing%KspsD(k,iCell)/BVF(k,iCell)))
      endif

      adcMixing%KmU(k,iCell) = 0.1_RKIND*adcMixing%lenspsU(k,iCell)*sqrt( adcMixing%KspsU(k,iCell) )
      adcMixing%KhU(k,iCell) = ( 1.0_RKIND+2.0_RKIND*adcMixing%lenspsU(k,iCell)/dz)*adcMixing%KmU(k,iCell)
      adcMixing%wt_spsU(k,iCell) =  -adcMixing%KhU(k,iCell)*Tz
      adcMixing%ws_spsU(k,iCell) =  -adcMixing%KhU(k,iCell)*Sz

      adcMixing%KmD(k,iCell) = 0.1_RKIND*adcMixing%lenspsD(k,iCell)*sqrt( adcMixing%KspsD(k,iCell) )
      adcMixing%KhD(k,iCell) = ( 1.0_RKIND+2.0_RKIND*adcMixing%lenspsD(k,iCell)/dz)*adcMixing%KmD(k,iCell)
      adcMixing%wt_spsD(k,iCell) = -adcMixing%KhD(k,iCell)*Tz
      adcMixing%ws_spsD(k,iCell) = -adcMixing%KhD(k,iCell)*Sz

      adcMixing%Entrainment(k,iCell) = adcMixing%Cww_E*adcMixing%areaFraction(k,iCell)*(1.0_RKIND- &
        adcMixing%areaFraction(k,iCell))*adcMixing%Mc(k,iCell) / ( adcMixing%lendn(k,iCell) + adcMixing%EPSILON )
      adcMixing%Detrainment(k,iCell) = adcMixing%Cww_D*adcMixing%areaFraction(k,iCell)*(1.0_RKIND- &
        adcMixing%areaFraction(k,iCell))*adcMixing%Mc(k,iCell) / ( adcMixing%lenup(k,iCell) + adcMixing%EPSILON )
    enddo

    do k=2,nVertLevels
      adcMixing%eps(k,iCell) = (0.5_RKIND*(adcMixing%u2(i2,k,iCell) + adcMixing%v2(i2,k,iCell) + &
        adcMixing%w2(i2,k,iCell)))**1.5_RKIND/adcMixing%length(k,iCell)

      dz = (adcMixing%zm(k-1,iCell) - adcMixing%zm(k,iCell))
      if(k==2) then
        Cval = 3.96
      else
        Cval = (0.19_RKIND+0.51_RKIND*adcMixing%lenspsU(k,iCell)/dz)
      endif
      KspsUtend(k,iCell) = adcMixing%grav*(alphaT(k,iCell)*adcMixing%wt_spsU(k,iCell) - &
        betaS(k,iCell)*adcMixing%ws_spsU(k,ICell)) + ((adcMixing%KmU(k-1,iCell) +                 &
        adcMixing%KmU(k,iCell))*(adcMixing%KspsU(k-1,iCell) - adcMixing%KspsU(k,iCell)) /         &
        (adcMixing%ze(k-1,iCell) - adcMixing%ze(k,iCell)) - (adcMixing%KmU(k,iCell) +             &
        adcMixing%KmU(k+1,iCell)) * (adcMixing%KspsU(k,iCell) - adcMixing%KspsU(k+1,iCell)) /     &
        (adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell))) / dz -                                 &
        Cval*adcMixing%KspsU(k,iCell)**1.5_RKIND/adcMixing%lenspsU(k,iCell) +                     &
        adcMixing%eps(k,iCell) / (2.0_RKIND*adcMixing%areaFraction(k,iCell))

      if(k==2) then
        Cval = 3.96_RKIND
      else
        Cval = (0.19_RKIND+0.51_RKIND*adcMixing%lenspsD(k,iCell)/dz)
      endif

      KspsDtend(k,iCell) = adcMixing%grav*(alphaT(k,iCell)*adcMixing%wt_spsD(k,iCell) - &
        betaS(k,iCell)*adcMixing%ws_spsD(k,iCell)) + ((adcMixing%KmD(k-1,iCell) +                 &
        adcMixing%KmD(k,iCell))*(adcMixing%KspsD(k-1,iCell) - adcMixing%KspsD(k,iCell)) /         &
        (adcMixing%ze(k-1,iCell) - adcMixing%ze(k,iCell)) - (adcMixing%KmD(k,iCell) +             &
        adcMixing%KmD(k+1,iCell)) * (adcMixing%KspsD(k,iCell) - adcMixing%KspsD(k+1,iCell)) /     &
        (adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell))) / dz -                                 &
        Cval*adcMixing%KspsD(k,iCell)**1.5_RKIND / adcMixing%lenspsD(k,iCell) +                   &
        adcMixing%eps(k,iCell) / (2.0_RKIND*(1.0_RKIND - adcMixing%areaFraction(k,iCell)))
    enddo
  enddo

  do iCell=1,nCells
    do k=2,nVertLevels
      adcMixing%KspsU(k,iCell) = adcMixing%KspsU(k,iCell) + dt*KspsUtend(k,iCell)
      adcMixing%KspsD(k,iCell) = adcMixing%KspsD(k,iCell) + dt*KspsDtend(k,iCell)
    enddo
  enddo

  end subroutine calc_subplume_fluxes

  subroutine diagnose_momentum_fluxes(nCells,nVertLevels,numTracers,activeTracers,uvel,vvel, &
      alphaT,betaS,dt)
! This routine diagnoses all the horizontal related momentum flux components. All assume steady state
! follows a quasi structure function approach

    integer,intent(in) :: numTracers,nCells,nVertLevels
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT,betaS
    real,dimension(nVertLevels,nCells),intent(in) :: uvel,vvel
    real,dimension(numTracers,nVertLevels,nCells),intent(in) :: activeTracers
    real,intent(in) :: dt
    real,dimension(nVertLevels,nCells) :: taupt, taups, taupv
    real,dimension(nVertLevels+1,nCells) :: wttend, wstend, uwtend, vwtend, u2tend, v2tend, &
      uttend, vttend, ustend, vstend, uvtend
    real :: B, Kps, Kpsp1, diff, lenav, Uz, Vz, Tz, Sz, sigav, sumdav
    real :: dz, dzmid, tumdav, Ksps
    integer :: iCell, k

    do iCell=1,nCells
      !compute the TOMs first.
      do k=1,nVertLevels
        Ksps = 0.5_RKIND*((adcMixing%areaFraction(k,iCell)*adcMixing%KspsU(k,iCell) + (1.0_RKIND - &
          adcMixing%areaFraction(k,iCell))*adcMixing%KspsD(k,iCell)) + (adcMixing%areaFraction(k+1,iCell)* &
          adcMixing%KspsU(k+1,iCell) + (1.0_RKIND-adcMixing%areaFraction(k+1,iCell))* &
          adcMixing%KspsD(k+1,iCell)))
        Kps = 0.5_RKIND*(adcMixing%u2(i1,k,iCell) + adcMixing%v2(i1,k,iCell) + adcMixing%w2(i2,k,iCell))
        Kpsp1 = 0.5_RKIND*(adcMixing%u2(i1,k+1,iCell) + adcMixing%v2(i1,k+1,iCell) + adcMixing%w2(i2,k+1,iCell))
        lenav = 0.5_RKIND*(adcMixing%length(k,iCell) + adcMixing%length(k+1,iCell))
        diff = adcMixing%C_mom * sqrt(0.5_RKIND*(Kps + Kpsp1)) / lenav
        dz = adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell)
        adcMixing%uw2(k,iCell) = -diff*(adcMixing%uw(i1,k,iCell) - adcMixing%uw(i1,k+1,iCell)) / dz
        adcMixing%vw2(k,iCell) = -diff*(adcMixing%vw(i1,k,iCell) - adcMixing%vw(i1,k+1,iCell)) / dz
        adcMixing%u2w(k,iCell) = -diff*(adcMixing%u2(i1,k,iCell) - adcMixing%u2(i1,k+1,iCell)) / dz
        adcMixing%v2w(k,iCell) = -diff*(adcMixing%v2(i1,k,iCell) - adcMixing%v2(i1,k+1,iCell)) / dz
        adcMixing%uvw(k,iCell) = -diff*(adcMixing%uv(i1,k,iCell) - adcMixing%uv(i1,k+1,iCell)) / dz

        diff = adcMixing%C_therm*sqrt(0.5*(Kps + Kpsp1)) / lenav
        adcMixing%uwt(k,iCell) = -diff*(adcMixing%ut(i1,k,iCell) - adcMixing%ut(i1,k+1,iCell)) / dz
        adcMixing%vwt(k,iCell) = -diff*(adcMixing%vt(i1,k,iCell) - adcMixing%vt(i1,k+1,iCell)) / dz
        adcMixing%uws(k,iCell) = -diff*(adcMixing%us(i1,k,iCell) - adcMixing%us(i1,k+1,iCell)) / dz
        adcMixing%vws(k,iCell) = -diff*(adcMixing%vs(i1,k,iCell) - adcMixing%vs(i1,k+1,iCell)) / dz
      enddo

      do k=2,nVertLevels
        dz = adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell)
        dzmid = adcMixing%zm(k-1,iCell) - adcMixing%zm(k,iCell)
        Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dz
        Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dz

        Kps = sqrt((adcMixing%u2(i1,k,iCell) + adcMixing%v2(i1,k,iCell) + adcMixing%w2(i2,k,iCell)))
        B = adcMixing%grav*(alphaT(k,iCell)*adcMixing%wt(i1,k,iCell) - betaS(k,iCell)*adcMixing%ws(i1,k,iCell))

        taupt(k,iCell) = Kps / (sqrt(2.0)*adcMixing%c_pt*adcMixing%length(k,iCell))
        taups(k,iCell) = Kps / (sqrt(2.0)*adcMixing%c_ps*adcMixing%length(k,iCell))

        !PDF formulation will allow for fairly easy generalization of vertical tracer fluxes to
        !all tracers, write in tumd, sumd, wumd form, not sure how to do horizontal fluxes though
        adcMixing%wttend1(k,iCell) = -(adcMixing%w2t(k-1,iCell) - adcMixing%w2t(k,iCell)) / dzmid
        adcMixing%wttend2(k,iCell) = -adcMixing%w2(i1,k,iCell)*(activeTracers(1,k-1,iCell) &
          - activeTracers(1,k,iCell)) / dzmid
        adcMixing%wttend3(k,iCell) = (1.0_RKIND - adcMixing%gamma1)*adcMixing%grav*(alphaT(k,iCell) &
          *adcMixing%t2(i2,k,iCell) - betaS(k,iCell)*adcMixing%ts(i2,k,iCell))
        adcMixing%wttend4(k,iCell) = -adcMixing%alpha3/4.0*(adcMixing%ut(i1,k,iCell)*Uz + &
          adcMixing%vt(i1,k,iCell)*Vz)
        adcMixing%wttend5(k,iCell) = adcMixing%kappa_FL*(adcMixing%wt(i1,k-1,iCell) -       &
              adcMixing%wt(i1,k+1,iCell)) / (adcMixing%ze(k-1,iCell) - adcMixing%ze(k+1,iCell))**2.0

        wttend(k,iCell) = adcMixing%wttend1(k,iCell) + adcMixing%wttend2(k,iCell) +        &
          adcMixing%wttend3(k,iCell) + adcMixing%wttend4(k,iCell) + adcMixing%wttend5(k,iCell)

        adcMixing%wstend1(k,iCell) = -(adcMixing%w2s(k-1,iCell) - adcMixing%w2s(k,iCell)) / dz
        adcMixing%wstend2(k,iCell) = -adcMixing%w2(i1,k,iCell)*(activeTracers(2,k-1,iCell) &
          - activeTracers(2,k,iCell)) / dz
        adcMixing%wstend3(k,iCell) = (1.0_RKIND - adcMixing%gamma1)*adcMixing%grav*(alphaT(k,iCell)* &
          adcMixing%ts(i2,k,iCell) - betaS(k,iCell)*adcMixing%s2(i2,k,iCell))
        adcMixing%wstend4(k,iCell) = -adcMixing%alpha3/4.0*(adcMixing%us(i1,k,iCell)*Uz + &
          adcMixing%vs(i1,k,iCell)*Vz)
        adcMixing%wsTend5(k,iCell) = adcMixing%kappa_FL*(adcMixing%ws(i1,k-1,iCell) -     &
          adcMixing%ws(i1,k+1,iCell)) / (adcMixing%ze(k-1,iCell) - adcMixing%ze(k+1,iCell))**2.0_RKIND

        wstend(k,iCell) = adcMixing%wstend1(k,iCell) + adcMixing%wstend2(k,iCell) + &
          adcMixing%wstend3(k,iCell) + adcMixing%wstend4(k,iCell) + adcMixing%wstend5(k,iCell)

        taupv(k,iCell) = Kps / (adcMixing%c_pv*adcMixing%length(k,iCell))

        adcMixing%uwtend1(k,iCell) = -(adcMixing%uw2(k-1,iCell) - adcMixing%uw2(k,iCell)) / dz
        adcMixing%uwtend2(k,iCell) = 0.5*((0.8-4.0/3.0*adcMixing%alpha1)*0.5*Kps**2.0 +  &
          (adcMixing%alpha1 - adcMixing%alpha2)*adcMixing%u2(i1,k,iCell) + (adcMixing%alpha1 +  &
          adcMixing%alpha2 - 2.0_RKIND)*adcMixing%w2(i1,k,iCell))*Uz
        adcMixing%uwtend3(k,iCell) = 0.5_RKIND*(adcMixing%alpha1 - adcMixing%alpha2)*    &
          adcMixing%uv(i1,k,iCell)*Vz
        adcMixing%uwtend4(k,iCell) = adcMixing%beta5*adcMixing%grav*(alphaT(k,iCell)*   &
          adcMixing%ut(i1,k,iCell) - betaS(k,iCell)*adcMixing%us(i1,k,iCell))
        adcMixing%uwtend5(k,iCell) = - 2.0_RKIND*taupv(k,iCell)*adcMixing%uw(i1,k,iCell) + &
          adcMixing%kappa_FL*(adcMixing%uw(i1,k-1,iCell) - adcMixing%uw(i1,k+1,iCell)) / &
          (adcMixing%ze(k-1,iCell) - adcMixing%ze(k+1,iCell))**2.0

        uwtend(k,iCell) = adcMixing%uwtend1(k,iCell) + adcMixing%uwtend2(k,iCell) + &
          adcMixing%uwtend3(k,iCell) + adcMixing%uwtend4(k,iCell) + adcMixing%uwtend5(k,iCell)

        vwtend(k,iCell) = (-(adcMixing%vw2(k-1,iCell) - adcMixing%vw2(k,iCell)) / dz +   &
          0.5_RKIND*((0.8_RKIND-4.0_RKIND/3.0_RKIND*adcMixing%alpha1)*0.5_RKIND*Kps**2.0_RKIND +  &
          (adcMixing%alpha1 - adcMixing%alpha2)*adcMixing%v2(i1,k,iCell) + (adcMixing%alpha1 -   &
          adcMixing%alpha2 - 2.0_RKIND)*adcMixing%w2(i1,k,iCell))*Vz + 0.5_RKIND*(adcMixing%alpha1 &
          - adcMixing%alpha2)*adcMixing%uv(i1,k,iCell)*Uz + adcMixing%beta5*adcMixing%grav*       &
          (alphaT(k,iCell)*adcMixing%vt(i1,k,iCell) - betaS(k,iCell)*adcMixing%vs(i1,k,iCell))) -            &
          taupv(k,iCell)*adcMixing%vw(i1,k,iCell) + adcMixing%kappa_FL*(adcMixing%vw(i1,k-1,iCell) &
          - adcMixing%vw(i1,k+1,iCell)) / (adcMixing%ze(k-1,iCell) - adcMixing%ze(k+1,iCell))**2.0

        uvtend(k,iCell) = (-(adcMixing%uvw(k-1,iCell) - adcMixing%uvw(k,iCell)) / dz + &
          (1.0_RKIND - 0.5_RKIND*(adcMixing%alpha1+adcMixing%alpha2))*(adcMixing%uw(i1,k,iCell)*Vz &
          + adcMixing%vw(i1,k,iCell)*Uz)) - taupv(k,iCell)*adcMixing%uv(i1,k,iCell) +           &
          adcMixing%kappa_VAR*(adcMixing%uv(i1,k-1,iCell) - adcMixing%uv(i1,k+1,iCell)) /       &
          (adcMixing%ze(k-1,iCell) - adcMixing%ze(k+1,iCell))**2.0_RKIND

        adcMixing%u2tend1(k,iCell) = -(adcMixing%u2w(k-1,iCell) - adcMixing%u2w(k,iCell)) / dz
        adcMixing%u2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*adcMixing%alpha1 + adcMixing%alpha2 - &
          2.0_RKIND)*adcMixing%uw(i1,k,iCell)*Uz
        adcMixing%u2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*adcMixing%alpha1*adcMixing%vw(i1,k,iCell)*Vz
        adcMixing%u2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1.0_RKIND-adcMixing%beta5)*B
        adcMixing%u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*adcMixing%eps(k,iCell) + &
          taupv(k,iCell)*(Kps**2.0_RKIND/3.0_RKIND - adcMixing%u2(i1,k,iCell))

        u2tend(k,iCell) = adcMixing%u2tend1(k,iCell) + adcMixing%u2tend2(k,iCell) + &
          adcMixing%u2tend3(k,iCell) + adcMixing%u2tend4(k,iCell) + adcMixing%u2tend5(k,iCell)

        adcMixing%v2tend1(k,iCell) = -(adcMixing%v2w(k-1,iCell) - adcMixing%v2w(k,iCell)) / dz
        adcMixing%v2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*adcMixing%alpha1 +adcMixing%alpha2 - &
          2.0_RKIND)*adcMixing%vw(i1,k,iCell)*Vz
        adcMixing%v2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*adcMixing%alpha1*adcMixing%uw(i1,k,iCell)*Uz
        adcMixing%v2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1-adcMixing%beta5)*B
        adcMixing%v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*adcMixing%eps(k,iCell) + taupv(k,iCell)* &
          (Kps**2/3. - adcMixing%v2(i1,k,iCell))

        v2tend(k,iCell) = adcMixing%v2tend1(k,iCell) + adcMixing%v2tend2(k,iCell) + &
          adcMixing%v2tend3(k,iCell) + adcMixing%v2tend4(k,iCell) + adcMixing%v2tend5(k,iCell)

        !taupt = Kps / (2.0*adcMixing%c_pt*len(k,iCell))
        !taups = Kps / (2.0*adcMixing%c_ps*len(k,iCell))

        Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dz
        Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dz

        uttend(k,iCell) = (-(adcMixing%uwt(k-1,iCell) - adcMixing%uwt(k,iCell))/dz -  &
          adcMixing%uw(i1,k,iCell)*Tz - (1.0_RKIND - adcMixing%alpha3)*adcMixing%wt(i1,k,iCell) &
          *Uz) - adcMixing%ut(i1,k,iCell)*taupt(k,iCell)

        vttend(k,iCell) = (-(adcMixing%vwt(k-1,iCell) - adcMixing%vwt(k,iCell))/dz -  &
          adcMixing%vw(i1,k,iCell)*Tz - (1.0_RKIND - adcMixing%alpha3)*adcMixing%wt(i1,k,iCell) &
          *Vz) - adcMixing%vt(i1,k,iCell)*taupt(k,iCell)

        ustend(k,iCell) = (-(adcMixing%uws(k-1,iCell) - adcMixing%uws(k,iCell))/dz -  &
          adcMixing%uw(i1,k,iCell)*Sz - (1.0_RKIND - adcMixing%alpha3)*adcMixing%ws(i1,k,iCell) &
          *Uz) - adcMixing%us(i1,k,iCell)*taups(k,ICell)

        vstend(k,iCell) = (-(adcMixing%vws(k-1,iCell) - adcMixing%vws(k,iCell))/dz -  &
          adcMixing%vw(i1,k,iCell)*Sz - (1.0_RKIND - adcMixing%alpha3)*adcMixing%ws(i1,k,iCell)* &
          Vz) - adcMixing%vs(i1,k,iCell)*taups(k,iCell)

        adcMixing%t2(i2,k,iCell) = adcMixing%tumd(k,iCell)**2.0_RKIND*adcMixing%areaFraction(k,iCell)*&
          (1.0_RKIND-adcMixing%areaFraction(k,iCell))
        adcMixing%s2(i2,k,iCell) = adcMixing%sumd(k,iCell)**2.0_RKIND*adcMixing%areaFraction(k,iCell)*&
          (1.0_RKIND-adcMixing%areaFraction(k,iCell))
        adcMixing%ts(i2,k,iCell) = adcMixing%tumd(k,iCell)*adcMixing%sumd(k,iCell)* &
          adcMixing%areaFraction(k,iCell)*(1.0_RKIND-adcMixing%areaFraction(k,iCell))

      enddo
    enddo

    do iCell=1,nCells
      adcMixing%u2cliptend(:,iCell) = 0.0_RKIND
      adcMixing%v2cliptend(:,iCell) = 0.0_RKIND
      do k=2,nVertLevels
        adcMixing%u2(i2,k,iCell) = adcMixing%u2(i1,k,iCell) + dt*u2tend(k,iCell)
        if(adcMixing%u2(i2,k,iCell) < 0.0_RKIND) then
          adcMixing%u2cliptend(k,iCell) = -adcMixing%u2(i2,k,iCell)
          adcMixing%u2(i2,k,iCell) = 0.0_RKIND
        endif

        adcMixing%v2(i2,k,iCell) = adcMixing%v2(i1,k,iCell) + dt*v2tend(k,iCell)
        if(adcMixing%v2(i2,k,iCell) < 0.0_RKIND) then
          adcMixing%v2cliptend(k,iCell) = -adcMixing%v2(i2,k,iCell)
          adcMixing%v2(i2,k,iCell) = 0.0_RKIND
        endif

        adcMixing%uw(i2,k,iCell) = adcMixing%uw(i1,k,iCell) + dt*uwtend(k,iCell)
        adcMixing%vw(i2,k,iCell) = adcMixing%vw(i1,k,iCell) + dt*vwtend(k,iCell)
        adcMixing%uv(i2,k,iCell) = adcMixing%uv(i1,k,iCell) + dt*uvtend(k,iCell)
        adcMixing%ut(i2,k,iCell) = adcMixing%ut(i1,k,iCell) + dt*uttend(k,iCell)
        adcMixing%wt(i2,k,iCell) = (adcMixing%wt(i1,k,iCell) + dt*wttend(k,iCell)) / &
          (1.0_RKIND + dt*taupt(k,iCell))
        adcMixing%vt(i2,k,iCell) = adcMixing%vt(i1,k,iCell) + dt*vttend(k,iCell)
        adcMixing%us(i2,k,iCell) = adcMixing%us(i1,k,iCell) + dt*ustend(k,iCell)
        adcMixing%vs(i2,k,iCell) = adcMixing%vs(i1,k,iCell) + dt*vstend(k,iCell)
        adcMixing%ws(i2,k,iCell) = (adcMixing%ws(i1,k,iCell) + dt*wstend(k,iCell)) / &
          (1.0_RKIND + dt*taups(k,iCell))
        if(abs(adcMixing%wt(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: wt out of range, wt = %r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/adcMixing%wt(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

        if(abs(adcMixing%ws(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: ws out of range, ws = %r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/adcMixing%ws(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

        if(abs(adcMixing%u2(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: u2 out of range, u2 = %r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/adcMixing%u2(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

        if(abs(adcMixing%v2(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: v2 out of range, v2 = %r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/adcMixing%v2(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

      enddo

    enddo

  end subroutine diagnose_momentum_fluxes

  subroutine predict_turbulent_quantities(nCells, nVertLevels, nTracers,  &
        dt, activeTracers, uvel,vvel,alphaT,betaS)
    integer,intent(in) :: nCells, nVertLevels, nTracers
    real,intent(in) :: dt
    real,dimension(nTracers,nVertLevels,nCells),intent(in) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS
    real,dimension(nVertLevels,nCells) :: uvel, vvel

    real,dimension(nVertLevels+1,nCells) :: w2tend
    real,dimension(nVertLevels,nCells) :: w3tend

    real :: Sw, St, Ss, Eav, Dav, sigav, sigavp1, wumdAv, tumdAv, sumdAv, wumdAvp1, tumdAvp1, sumdAvp1
    real :: Swup, KspsUav, KspsDav, KspsUavp1, KspsDavp1, KE, Mcav, lenav,u2av,v2av,w2av
    real :: w3temp, w3check, taups, taupt, mval, KEsps, Uz, Vz, dz

    integer :: iCell, k

    do iCell = 1,nCells
      do k=1,nVertLevels
        Eav = 0.5*(adcMixing%Entrainment(k+1,iCell) + adcMixing%Entrainment(k,iCell))
        Dav = 0.5*(adcMixing%Detrainment(k+1,iCell) + adcMixing%Detrainment(k,iCell))
        u2av = 0.5*(adcMixing%u2(i1,k,iCell) + adcMixing%u2(i1,k+1,iCell))
        v2av = 0.5*(adcMixing%v2(i1,k,iCell) + adcMixing%v2(i1,k+1,iCell))
        w2av = 0.5*(adcMixing%w2(i1,k,iCell) + adcMixing%w2(i1,k+1,iCell))

        sigav = 0.5*(adcMixing%areaFraction(k,iCell) + adcMixing%areaFraction(k+1,iCell))
        wumdav = 0.5*(adcMixing%wumd(k,iCell) + adcMixing%wumd(k+1,iCell))
        tumdav = 0.5*(adcMixing%tumd(k,iCell) + adcMixing%tumd(k+1,iCell))
        sumdav = 0.5*(adcMixing%sumd(k,iCell) + adcMixing%sumd(k+1,iCell))
        KspsUav = 0.5*(adcMixing%KspsU(k,iCell) + adcMixing%KspsU(k+1,iCell))
        KspsDav = 0.5*(adcMixing%KspsD(k,iCell) + adcMixing%KspsD(k+1,iCell))
        Mcav = 0.5*(adcMixing%Mc(k,iCell) + adcMixing%Mc(k+1,iCell))
        lenav = 0.5*(adcMixing%length(k,iCell) + adcMixing%length(k+1,iCell))
        if(k==nVertLevels) then
          sigavp1 = 0.5*(adcMixing%areaFraction(k,iCell))
          wumdAvp1 = 0.5*(adcMixing%wumd(k,iCell))
          tumdAvp1 = 0.5*(adcMixing%tumd(k,iCell))
          sumdAvp1 = 0.5*(adcMixing%sumd(k,iCell))
        else
          sigavp1 = 0.5*(adcMixing%areaFraction(k,iCell) + adcMixing%areaFraction(k+1,iCell))
          wumdAvp1 = 0.5*(adcMixing%wumd(k,iCell) + adcMixing%wumd(k+1,iCell))
          tumdAvp1 = 0.5*(adcMixing%tumd(k,iCell) + adcMixing%tumd(k+1,iCell))
          sumdAvp1 = 0.5*(adcMixing%sumd(k,iCell) + adcMixing%sumd(k+1,iCell))
        endif

        KEsps = sigav*KspsUav+ (1.0 - sigav)*KspsDav
        KE = sqrt((u2av+v2av+w2av) + 0.0*KEsps)

        dz = adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell)
        !KE = sqrt(areaFraction(k,iCell)*KspsUav + (1.0 - areaFraction(k,iCell))*KspsDav)
        Swup = - 2.0/3.0*(adcMixing%KspsU(k,iCell) - adcMixing%KspsU(k+1,iCell)) / dz &
           - 2.0/3.0*KspsUav*(log(adcMixing%areaFraction(k,iCell)) -           &
          log(adcMixing%areaFraction(k+1,iCell))) / dz +         &
          2.0/3.0*(adcMixing%KspsD(k,iCell) - adcMixing%KspsD(k+1,iCell)) / dz        &
         + 2.0/3.0*KspsDav*(log(1.0-adcMixing%areaFraction(k,iCell)) -         &
          log(1.0-adcMixing%areaFraction(k+1,iCell))) /  dz

        adcMixing%w3tend1(k,iCell) = wumdav**3.0*(Eav*(3.0*sigav - 2.0) + Dav*(3.0*sigav - 1.0))
        adcMixing%w3tend2(k,iCell) =  wumdav**3.0*(6.0*sigav**2.0 - 6.0*sigav + 1)*             &
              (adcMixing%areaFraction(k,iCell)*(1.0_RKIND - adcMixing%areaFraction(k,iCell))*  &
              adcMixing%wumd(k,iCell) - adcMixing%areaFraction(k+1,iCell)*(1.0_RKIND - &
              adcMixing%areaFraction(k+1,iCell))*adcMixing%wumd(k+1,iCell))/dz
        adcMixing%w3tend3(k,iCell) = -1.5_RKIND*sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)* &
              wumdav**2.0_RKIND*((1.0_RKIND - 2.0_RKIND*adcMixing%areaFraction(k,iCell))*      &
              adcMixing%wumd(k,iCell)**2.0_RKIND - (1.0_RKIND -                                &
              2.0_RKIND*adcMixing%areaFraction(k+1,iCell))*adcMixing%wumd(k+1,iCell)**2.0_RKIND) / dz
        adcMixing%w3tend4(k,iCell) = 3.0_RKIND*(1.0_RKIND - 2.0_RKIND*sigav)*Mcav*wumdav*Swup -          &
              adcMixing%C_mom_w3*KE / (1E-15_RKIND+ sqrt(2.0_RKIND)*lenAv)*adcMixing%w3(i1,k,iCell)
        adcMixing%w3tend5(k,iCell) = 3.0_RKIND*adcMixing%grav*(alphaT(k,iCell)*adcMixing%w2t(k,iCell) - &
              betaS(k,iCell)*adcMixing%w2s(k,iCell))*0.9_RKIND

        w3tend(k,iCell) = adcMixing%w3tend1(k,iCell) + adcMixing%w3tend2(k,iCell) + adcMixing%w3tend3(k,iCell) + &
              adcMixing%w3tend4(k,iCell) + adcMixing%w3tend5(k,iCell)

        if(k>1 .and. k < nVertLevels) then
           w3tend(k,iCell) = w3tend(k,iCell) + adcMixing%kappa_w3*(adcMixing%w3(i1,k-1,iCell) &
              - adcMixing%w3(i1,k+1,iCell)) / (adcMixing%zm(k-1,iCell) - adcMixing%zm(k+1,iCell))**2.0
        endif

      enddo

!      k=1
!      w3check = (w2(i1,k,iCell)+w2(i1,k+1,iCell))**1.5
!      w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)
!      do k=2,nVertLevels
!        w3check = (w2(i1,k,iCell) + w2(i1,k+1,iCell))**1.5
!        w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)
!      enddo

      do k=2,nVertLevels-1
        dz = (adcMixing%zm(k-1,iCell) - adcMixing%zm(k,iCell))

        sigav = 0.5*(adcMixing%areaFraction(k,iCell) + adcMixing%areaFraction(k-1,iCell))
        wumdav = 0.5*(adcMixing%wumd(k,iCell) + adcMixing%wumd(k-1,iCell))
        tumdav = 0.5*(adcMixing%tumd(k,iCell) + adcMixing%tumd(k-1,iCell))
        sumdav = 0.5*(adcMixing%sumd(k,iCell) + adcMixing%sumd(k-1,iCell))
        KspsUav = 0.5*(adcMixing%KspsU(k,iCell) + adcMixing%KspsU(k-1,iCell))
        KspsDav = 0.5*(adcMixing%KspsD(k,iCell) + adcMixing%KspsD(k-1,iCell))
        Mcav = 0.5*(adcMixing%Mc(k,iCell) + adcMixing%Mc(k-1,iCell))

        sigavp1 = 0.5*(adcMixing%areaFraction(k,iCell) + adcMixing%areaFraction(k+1,iCell))
        KspsUavp1 = 0.5*(adcMixing%KspsU(k,iCell) + adcMixing%KspsU(k+1,iCell))
        KspsDavp1 = 0.5*(adcMixing%KspsD(k,iCell) + adcMixing%KspsD(k+1,iCell))

        Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dz
        Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dz

        KEsps = adcMixing%areaFraction(k,iCell)*adcMixing%KspsU(k,iCell) +  &
                  (1.0_RKIND-adcMixing%areaFraction(k,iCell))*adcMixing%KspsD(k,iCell)
        KE = sqrt((adcMixing%u2(i1,k,iCell) + adcMixing%v2(i1,k,iCell) + adcMixing%w2(i1,k,iCell)) + 0.0*KEsps)
        Swup = adcMixing%grav*alphaT(k,iCell)*adcMixing%tumd(k,iCell) - adcMixing%grav*        &
          betaS(k,iCell)*adcMixing%sumd(k,iCell) - 2.0_RKIND/3.0_RKIND*(1.0_RKIND/        &
          adcMixing%areaFraction(k,iCell)*(sigAv*KspsUav - sigavp1*KspsUavp1) / dz -     &
          1.0_RKIND/(1.0_RKIND - adcMixing%areaFraction(k,iCell))*((1.0_RKIND - sigav)* &
          KspsDav - (1.0 - sigavp1)*KspsDavp1) / dz)

        adcMixing%w2tend1(k,iCell) = -adcMixing%wumd(k,iCell)**2.0_RKIND*(         &
          adcMixing%Entrainment(k,iCell) + adcMixing%Detrainment(k,iCell))
        adcMixing%w2tend2(k,iCell) = -(adcMixing%w3(i2,k-1,iCell) - adcMixing%w3(i2,k,iCell)) / dz
        adcMixing%w2tend3(k,iCell) = -adcMixing%C_1*KE / (1.0E-15_RKIND + sqrt(2.0_RKIND)* &
          adcMixing%length(k,iCell))*(adcMixing%w2(i1,k,iCell) - KE**2.0_RKIND/3.0_RKIND)
        adcMixing%w2tend4(k,iCell) = 2.0_RKIND*adcMixing%Mc(k,iCell)*Swup + 4.0_RKIND / &
          3.0_RKIND*adcMixing%C_2*adcMixing%areaFraction(k,iCell)*(1.0_RKIND -          &
          adcMixing%areaFraction(k,iCell))*adcMixing%wumd(k,iCell)*(adcMixing%grav*     &
          alphaT(k,iCell)*adcMixing%tumd(k,iCell) - adcMixing%grav*betaS(k,iCell)*      &
          adcMixing%sumd(k,iCell))
        adcMixing%w2tend5(k,iCell) = (1.0_RKIND/3.0_RKIND*adcMixing%alpha1 -            &
          adcMixing%alpha2)*(adcMixing%uw(i1,k,iCell)*Uz + adcMixing%vw(i1,k,iCell)*Vz)

        w2tend(k,iCell) = adcMixing%w2tend1(k,iCell) + adcMixing%w2tend2(k,iCell) + &
          adcMixing%w2tend3(k,iCell) + adcMixing%w2tend4(k,iCell) + adcMixing%w2tend5(k,iCell)
    enddo
  end do

    do iCell=1,nCells
      do k=2,nVertLevels
        adcMixing%w2(i2,k,iCell) = adcMixing%w2(i1,k,iCell) + dt*w2tend(k,iCell)
        if(adcMixing%w2(i2,k,iCell) < 0.0_RKIND) then
          adcMixing%w2cliptend(k,iCell) = -adcMixing%w2(i2,k,iCell)
          adcMixing%w2(i2,k,iCell) = 0.0_RKIND
        endif

        if(abs(adcMixing%w2(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: w2 out of range, w2 = %r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/adcMixing%w2(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif
      end do

      do k=1,nVertLevels
        w3check = (adcMixing%w2(i2,k,iCell) + adcMixing%w2(i2,k+1,iCell))**1.5
        adcMixing%w3(i2,k,iCell) = min(adcMixing%w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)

        if(abs(adcMixing%w3(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: w3 out of range, w2 = %r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/adcMixing%w3(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif
      enddo
    enddo

  end subroutine predict_turbulent_quantities

  subroutine update_mean_fields(dt,nCells,nVertLevels,numTracers,activeTracers,uvel,vvel,fCell)

    integer,intent(in) :: nCells, nVertLevels, numTracers
    real,dimension(nVertLevels,nCells),intent(out) :: uvel,vvel
    real,dimension(2,nvertLevels,nCells),intent(out) :: activeTracers
    real,intent(in) :: dt
    real,dimension(nCells),intent(in) :: fCell

    real :: utemp, vtemp
    integer :: iCell,k

    do iCell = 1,nCells
      do k = 1,nVertLevels
        utemp = uvel(k,iCell)
        vtemp = vvel(k,iCell)
        uvel(k,iCell) = uvel(k,iCell) - dt*(adcMixing%uw(i2,k,iCell) - adcMixing%uw(i2,k+1,iCell)) /  &
                  (adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell)) !+ dt*fCell(iCell)*vtemp

        vvel(k,iCell) = vvel(k,iCell) - dt*(adcMixing%vw(i2,k,iCell) - adcMixing%vw(i2,k+1,iCell)) /  &
                  (adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell)) !- dt*fCell(iCell)*utemp

        activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt*(adcMixing%wt(i2,k,iCell) - &
          adcMixing%wt(i2,k+1,iCell)) / (adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell))
        activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt*(adcMixing%ws(i2,k,iCell) - &
          adcMixing%ws(i2,k+1,iCell)) / (adcMixing%ze(k,iCell) - adcMixing%ze(k+1,iCell))
      enddo
    enddo

  end subroutine update_mean_fields

  subroutine ADC_main_loop(nCells,nVertLevels, nTracers, dt,activeTracers, uvel, vvel, BVF,  &
      uwsfc, vwsfc, wtsfc, wssfc, alphaT, betaS, fCell, boundaryLayerDepth)

    integer,intent(in) :: nCells, nVertLevels, nTracers
    real,intent(in) :: dt

    real,dimension(nTracers,nVertLevels,nCells),intent(inout) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(inout) :: uvel, vvel, alphaT, betaS
    real,dimension(nCells),intent(in) :: uwsfc, vwsfc, wtsfc, wssfc, fCell
    real,dimension(nCells),intent(inout) :: boundaryLayerDepth
    real,dimension(nVertLevels+1,nCells),intent(inout) :: BVF
    integer :: niter, iIter,iCell,k

    real :: dt_small

    dt_small = config_adc_timestep
    niter = dt / dt_small

    i1 = 1
    i2 = 2
!    call construct_depth_coordinate(ssh,layerThick,nCells,nVertLevels)
    do iIter=1,niter
      call build_diagnostic_arrays(nCells,nVertLevels,BVF,wtsfc,wssfc,  &
        uwsfc,vwsfc,alphaT,betaS,boundaryLayerDepth)
      call predict_turbulent_quantities(nCells, nVertLevels, nTracers, dt_small, &
        activeTracers, uvel, vvel, alphaT, betaS)
      call diagnose_momentum_fluxes(nCells,nVertLevels,nTracers,activeTracers,uvel,vvel,alphaT,betaS,dt_small)
      call build_areaFraction_updraft_properties(nCells, nVertLevels)
      call calc_scalar_updraft_properties(nCells, nVertLevels,wtsfc, wssfc, &
                                alphaT, betaS, boundaryLayerDepth, i2)
      call calc_subplume_fluxes(nCells,nVertLevels,activeTracers,uvel,vvel, BVF,alphaT,betaS,dt_small)
!      call build_dissipation_lengths(nCells,nVertLevels,BVF)
      call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
      call update_mean_fields(dt_small,nCells,nVertLevels, 2, uvel,vvel,activeTracers,fCell)
      call swap_time_levels
    enddo
  end subroutine ADC_main_loop

end module ocn_adc_mixing
