module adc

  use mpas_ocn_constants

  implicit none

  logical :: defineFirst, stopflag

  integer :: record
  !declare all variables up here
  !move all this into registry
  !maybe have it in registry but then store in local memory for length of the run
  !exception would be wt, ws, uw, vw, maybe also consider as kappa + the rest
  real,allocatable,dimension(:,:) :: zmid, zedge, KspsU, KspsD, eps, epsnew, len, lenspsD, &
    lenup, lendn, lenB, lenbuoy, lenblac, lenshea, tau, taunew, &
    KhU, KhD, KmU, KmD, wt_spsU, wt_spsD, ws_spsU, ws_spsD, lenspsU, &
    sigma, kappa, E, D, Ktend, KspsUtend, KspsDtend, Ktend1, &
    Ktend2, Ktend3, w2tend, w2tend1, w2tend2, w2tend3, w2tend4, &
    w2tend5, w3tend, w3tend1, w3tend2, w3tend3, w3tend4, w3tend5, &
    wttend, wttend1, wttend2, wttend3, wttend4, wttend5, &
    wttend6, t2tend, t2tend1, t2tend2, t2tend3,  &
    wstend, s2tend, tstend, tumd, sumd, wumd, Mc, uw2, vw2, u2w, &
    wstend1, wstend2, wstend3, wstend4, wstend5, &
    v2w, w2t, w2s, wts, uvw, uwt, vwt, uws, vws, ws2, wt2,      &
    uwtend,vwtend,u2tend,v2tend,ustend,vstend,uttend,vttend,    &
    uvtend,u2tend1,u2tend2,u2tend3,u2tend4,u2tend5,u2tend6,     &
    uwtend1,uwtend2,uwtend3,uwtend4,uwtend5, u2cliptend,        &
    v2cliptend, w2cliptend,v2tend1,v2tend2,v2tend3,v2tend4,v2tend5

  real,allocatable,dimension(:,:,:) :: u2,v2,uw,vw,uv,w2,w3,w4, &
    wt,t2,ut,vt,ws,s2,us,vs,ts,b2

  real,allocatable,dimension(:) :: boundaryLayerDepth

  real,parameter :: EPSILON = 1.0E-8

  real :: fileTime

  integer :: i1,i2

!  integer,parameter :: ntimes = 2

  contains

  subroutine init_adc(ntimes,nCells, nVertLevels)

  integer, intent(in) :: ntimes, nCells, nVertLevels

  integer :: k,iCell

!  allocate all the variables, set some to zero and such
  allocate(zmid(nVertLevels,nCells))
  allocate(zedge(nVertLevels+1,nCells),lenblac(nVertLevels+1,nCells))
  allocate(KspsU(nVertLevels+1,nCells),KspsD(nVertLevels+1,nCells))
  allocate(eps(nVertLevels+1,nCells),epsnew(nVertLevels+1,nCells))
  allocate(len(nVertLevels+1,nCells),lenspsD(nVertLevels+1,nCells))
  allocate(lenspsU(nVertLevels+1,nCells))
  allocate(lenup(nVertLevels+1,nCells),lendn(nVertLevels+1,nCells))
  allocate(lenbuoy(nVertLevels+1,nCells),lenshea(nVertLevels+1,nCells))
  allocate(tau(nVertLevels+1,nCells),taunew(nVertLevels+1,nCells))
  allocate(KhU(nVertLevels+1,nCells), KhD(nVertLevels+1,nCells))
  allocate(KmU(nVertLevels+1,nCells), KmD(nVertLevels+1,nCells))
  allocate(u2(ntimes,nVertLevels+1,nCells), v2(ntimes,nVertLevels+1,nCells))
  allocate(uw(ntimes,nVertLevels+1,nCells), vw(ntimes,nVertLevels+1,nCells))
  allocate(uv(ntimes,nVertLevels+1,nCells), w2(ntimes,nVertLevels+1,nCells))
  allocate(w3(ntimes,nVertLevels+1,nCells), w4(ntimes,nVertLevels+1,nCells))
  allocate(wt(ntimes,nVertLevels+1,nCells), t2(ntimes,nVertLevels+1,nCells))
  allocate(ut(ntimes,nVertLevels+1,nCells), vt(ntimes,nVertLevels+1,nCells))
  allocate(wt_spsU(nVertLevels+1,nCells), wt_spsD(nVertLevels+1,nCells))
  allocate(ws_spsU(nVertLevels+1,nCells), ws_spsD(nVertLevels+1,nCells))
  allocate(ws(ntimes,nVertLevels+1,nCells), s2(ntimes,nVertLevels+1,nCells))
  allocate(us(ntimes,nVertLevels+1,nCells), vs(ntimes,nVertLevels+1,nCells))
  allocate(ts(ntimes,nVertLevels+1,nCells), b2(ntimes,nVertLevels+1,nCells))
  allocate(uw2(nVertLevels,nCells), vw2(nVertLevels,nCells))
  allocate(uvw(nVertLevels,nCells), u2w(nVertLevels,nCells))
  allocate(v2w(nVertLevels,nCells), uwt(nVertLevels,nCells))
  allocate(vwt(nVertLevels,nCells), w2t(nVertLevels,nCells))
  allocate(uws(nVertLevels,nCells), w2s(nVertLevels,nCells))
  allocate(ws2(nVertLevels,nCells), wts(nVertLevels,nCells))
  allocate(wumd(nVertLevels+1,nCells), tumd(nVertLevels+1,nCells))
  allocate(sumd(nVertLevels+1,nCells), Mc(nVertLevels+1,nCells))
  allocate(vws(nVertLevels,nCells), u2cliptend(nVertLevels+1,nCells))
  allocate(v2cliptend(nVertLevels+1,nCells),w2cliptend(nVertLevels+1,nCells))

  allocate(sigma(nVertLevels+1,nCells),kappa(nVertLevels+1,nCells))
  allocate(E(nVertLevels+1,nCells))
  allocate(D(nVertLevels+1,nCells), Ktend(nVertLevels+1,nCells))
  allocate(KspsUtend(nVertLevels+1,nCells), KspsDtend(nVertLevels+1,nCells))
  allocate(Ktend1(nVertLevels+1,nCells), Ktend2(nVertLevels+1,nCells))
  allocate(ktend3(nVertLevels+1,nCells), w2tend(nVertLevels+1,nCells))
  allocate(w2tend1(nVertLevels+1,nCells), w2tend2(nVertLevels+1,nCells))
  allocate(w2tend3(nVertLevels+1,nCells), w2tend4(nVertLevels+1,nCells))
  allocate(w2tend5(nVertLevels+1,nCells), w3tend1(nVertLevels+1,nCells))
  allocate(w3tend2(nVertLevels+1,nCells), w3tend(nVertLevels+1,nCells))
  allocate(w3tend3(nVertLevels+1,nCells), w3tend4(nVertLevels+1,nCells))
  allocate(w3tend5(nVertLevels+1,nCells), wttend(nVertLevels+1,nCells))
  allocate(wttend1(nVertLevels+1,nCells), wttend2(nVertLevels+1,nCells))
  allocate(wttend3(nVertLevels+1,nCells), wttend4(nVertLevels+1,nCells))
  allocate(wttend5(nVertLevels+1,nCells), wttend6(nVertLevels+1,nCells))
  allocate(t2tend(nVertLevels+1,nCells), t2tend1(nVertLevels+1,nCells))
  allocate(t2tend2(nVertLevels+1,nCells), t2tend3(nVertLevels+1,nCells))
  allocate(wstend(nVertLevels+1,nCells), s2tend(nVertLevels+1,nCells))
  allocate(wstend1(nVertLevels+1,nCells), wstend2(nVertLevels+1,nCells))
  allocate(wstend3(nVertLevels+1,nCells), wstend4(nVertLevels+1,nCells))
  allocate(wstend5(nVertLevels+1,nCells))
  allocate(tstend(nVertLevels+1,nCells), uwtend(nVertLevels+1,nCells))
  allocate(vwtend(nVertLevels+1,nCells), u2tend(nVertLevels+1,nCells))
  allocate(v2tend(nVertLevels+1,nCells), uvtend(nVertLevels+1,nCells))
  allocate(uttend(nVertLevels+1,nCells), vttend(nVertLevels+1,nCells))
  allocate(ustend(nVertLevels+1,nCells), vstend(nVertLevels+1,nCells))
  allocate(u2tend1(nVertLevels+1,nCells), u2tend2(nVertLevels+1,nCells))
  allocate(u2tend3(nVertLevels+1,nCells), u2tend4(nVertLevels+1,nCells))
  allocate(u2tend5(nVertLevels+1,nCells), u2tend6(nVertLevels+1,nCells))
  allocate(boundaryLayerDepth(nCells))
  allocate(v2tend1(nVertLevels+1,nCells), v2tend2(nVertLevels+1,nCells))
  allocate(v2tend3(nVertLevels+1,nCells), v2tend4(nVertLevels+1,nCells))
  allocate(v2tend5(nVertLevels+1,nCells))
  allocate(uwtend1(nVertLevels+1,nCells), uwtend2(nVertLevels+1,nCells))
  allocate(uwtend3(nVertLevels+1,nCells), uwtend4(nVertLevels+1,nCells))
  allocate(uwtend5(nVertLevels+1,nCells))

  do iCell=1,nCells

    do k = 1, nVertLevels
      KspsU(k,iCell) = EPSILON
      KspsD(k,iCell) = EPSILON
      w2(:,k,iCell) = 0.0
      sigma(k,iCell) = 0.5
      tumd(k,iCell) = 0.0
      sumd(k,iCell) = 0.0
      wumd(k,iCell) = 0.0
      uw(:,k,iCell) = 0.0
      vw(:,k,iCell) = 0.0
      u2(:,k,iCell) = 0.0
      v2(:,k,iCell) = 0.0
      uv(:,k,iCell) = 0.0
      ut(:,k,iCell) = 0.0
      vt(:,k,iCell) = 0.0
      us(:,k,iCell) = 0.0
      vs(:,k,iCell) = 0.0
      len(k,iCell)= 2.0*0.4
    enddo

  enddo

  do iCell=1,nCells
    w2(:,1,iCell) = 0.0
    KspsU(nVertLevels+1,iCell) = EPSILON
    KspsD(nVertLevels+1,iCell) = EPSILON
    sigma(nVertLevels+1,iCell) = 0.5
  enddo

  i1 = 1
  i2 = 2
  defineFirst = .true.
  end subroutine init_adc

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

  subroutine construct_depth_coordinate(ssh,layerThick,nCells,nVertLevels)
  ! builds a coordinate where zero is the ssh

    integer :: nCells, nVertLevels
    real,dimension(nVertLevels,nCells), intent(in) :: layerThick
    real,dimension(nCells),intent(in) :: ssh

    integer :: iCell, k

    do iCell=1,nCells
      zedge(1,iCell) = ssh(iCell)
      do k=2,nVertLevels+1
        zedge(k,iCell) = zedge(k-1,iCell) - layerThick(k-1,iCell)
        zmid(k-1,iCell) = 0.5*(zedge(k,iCell) + zedge(k-1,iCell))
      enddo
    enddo

  end subroutine construct_depth_coordinate

  subroutine build_diagnostic_arrays(nCells,nVertLevels,temp,salt,BVF,wtsfc,wssfc, &
        uwsfc, vwsfc, alphaT,betaS,adcConst)
  !construct dTdz, dSdz, dbdz
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, alphaT, betaS, uwsfc, vwsfc
    real,dimension(nVertLevels,nCells),intent(in) :: temp, salt
    type(adc_mixing_constants) :: adcConst
    real,dimension(nVertLevels+1,nCells),intent(out) :: BVF
    integer :: iCell, k, idx
    real,dimension(nCells) :: wstar
    logical :: first
    real :: maximum, Tz, Sz, Bz, Q

    first = .true.

    do iCell=1,nCells
      maximum = -1.0e-12
      idx = 1

      !Change this to pass in....
      BVF(1,iCell) = 0.0
      do k=2,nVertLevels
        Tz = (temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        Sz = (salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        BVF(k,iCell) = max(0.0, adcConst%grav*(alphaT(iCell)*Tz - betaS(iCell)*Sz))

        if(BVF(k,iCell) > 1.005*maximum .and. first) then
          maximum = BVF(k,iCell)
          idx = k
        elseif(BVF(k,iCell) < maximum) then
          first = .false.
        endif
      enddo

      boundaryLayerDepth(iCell) = abs(zedge(idx,iCell))
      Q = adcConst%grav*(alphaT(iCell)*wtsfc(iCell) - betaS(iCell)*wssfc(iCell))* &
        boundaryLayerDepth(iCell)
      if(Q > 0) then
        wstar(iCell) = abs(Q)**(1.0/3.0)
      else
        wstar(iCell) = 0.0
      endif

      u2(:,1,iCell) = 4.0*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      v2(:,1,iCell) =4.*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      uw(:,1,iCell) = -uwsfc(iCell)
      vw(:,1,iCell) = vwsfc(iCell)
      wt(:,1,iCell) = wtsfc(iCell)
      ws(:,1,iCell) = wssfc(iCell)
  !    print *, 'bld = ',boundaryLayerDepth(iCell)
    enddo

  end subroutine build_diagnostic_arrays

  subroutine dissipation_lengths2(nCells,nVertLevels,temp,salt,alphaT,betaS,zedge)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nVertLevels,nCells),intent(in) :: temp,salt
    real,dimension(nVertLevels+1,nCells),intent(in) :: zedge
    real,dimension(nCells),intent(in) :: alphaT, betaS
    integer :: i,k, ij

    real,dimension(nVertLevels) :: B, Bup, Bdo
    real,dimension(nVertLevels+1) :: tke, BupEdge, BdoEdge
    real :: sav, tudav, sudav, Tup, Tdo, Sup, Sdo
    real :: s1, z1, zV, sumv, minlen

    do i=1,nCells
       tke(:) = 0.5*(u2(i2,:,i) + v2(i2,:,i) + w2(i2,:,i))
       do k=1,nVertLevels
          B(k) = -9.806*(-alphaT(i)*(temp(k,i) - 15.0) + betaS(i)*   &
                    (salt(k,i) - 35.0))

          sav = 0.5*(sigma(k,i) + sigma(k+1,i))
          tudav = 0.5*(tumd(k,i) + tumd(k+1,i))
          sudav = 0.5*(sumd(k,i) + sumd(k+1,i))

          Tup = temp(k,i) + (1.0 - sav)*tudav
          Tdo = temp(k,i) - sav*tudav 
          Sup = salt(k,i) + (1.0 - sav)*sudav
          Sdo = salt(k,i) - sav*sudav

          Bup(k) = -9.806*(-alphaT(i)*(Tup - 15.0) + betaS(i)*(Sup - 35.0))
          Bdo(k) = -9.806*(-alphaT(i)*(Tdo - 15.0) + betaS(i)*(Sdo - 35.0))

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
          lenup(k,i) = 0
          do while(sumv <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv + (BupEdge(k) - Bup(ij))*(zedge(ij,i)-zedge(ij+1,i))**2/2.
             lenup(k,i) =  lenup(k,i) + abs(zedge(ij,i)-zedge(ij+1,i))
             ij = ij + 1

             if(sumv > tke(k)) THEN
                ij = ij - 1
!                s1 = sumv
!                z1 = zedge(ij+1,i)
!                zV = zedge(ij,i)
                sumv = sumv - (BupEdge(k) - Bup(ij))*(zedge(ij,i)-zedge(ij+1,i))**2/2.
                lenup(k,i) = lenup(k,i) - abs(zedge(ij,i)-zedge(ij+1,i))
!                lenup(k,i) = max(0.55,lenup(k,i) + abs((z1-zV)/(s1 - sumv)*(tke(k)-sumv)))
                if(Bup(k-1) - Bup(k) < 0) then 
                        minlen = abs(zmid(k-1,i) - zmid(k,i))
                else
                        minlen = 0.5
                endif
                lenup(k,i) = max(minlen, lenup(k,i) + sqrt(2.0/(BupEdge(k) -         &
                                Bup(ij-1))*(tke(k) - sumv)))
                exit   
             endif

         end do

        !find lendown
        sumv = 0
        ij=k
        lendn(k,i) = 0
        do while(sumv <= tke(k) .and. ij>1)
           sumv = sumv - (BdoEdge(k) - Bdo(ij-1))*(zedge(ij-1,i)-zedge(ij,i))**2./2.
           lendn(k,i) = lendn(k,i) + abs(zedge(ij-1,i)-zedge(ij,i))
           ij = ij - 1

           if(sumv > tke(k)) THEN
              ij = ij + 1
!              s1 = sumv
!              z1 = zedge(ij,i)
!              zV = zedge(ij-1,i)
              sumv = sumv + (BdoEdge(k) - Bdo(ij-1))*(zedge(ij-1,i)-zedge(ij,i))**2./2.
              lendn(k,i) = lendn(k,i) - abs(zedge(ij-1,i)-zedge(ij,i))
!              lendown(k) = max(0.55,lendown(k) + abs(-(z1-zV)/(sumv)*(tke(k))))
                if(Bdo(k-1) - Bdo(k) < 0) then
                        minlen = abs(zmid(k-1,i) - zmid(k,i))
                else
                        minlen = 0.55
                endif
               lendn(k,i) = max(minlen,lendn(k,i) + sqrt(2.0/(BdoEdge(k) -  &
                                Bdo(ij-1))*(tke(k) - sumv)))
              exit
           endif
        enddo

        len(k,i) = (2.0*lenup(k,i)*lendn(k,i)) / (lenup(k,i) + lendn(k,i))
      enddo
   enddo

   len(1,:) = 0.55
   len(nVertLevels+1,:) = 0.55
 
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
        KE = sqrt(0.5*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))) 
        integrandTop = integrandTop + 0.5*KE*abs(zedge(k,iCell)**2 - zedge(k+1,iCell)**2)
        integrandBot = integrandBot + KE*(zedge(k,iCell) - zedge(k+1,iCell))
      enddo

      do k=2,nVertLevels
        KE = 0.5*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))

        if(KE > EPSILON) then
          lenmax = 0.53*sqrt(2.0*KE / (1.0E-15 + BVF(k,iCell)))
        else
          lenmax = 1.0e6
        endif

        len1 = 0.4*abs(zedge(k,iCell))
        len2 = 0.2*integrandTop / (integrandBot + 1.0E-10)

        lenbuoy(k,iCell) = len1
        lenblac(k,iCell) = len2
        lenshea(k,iCell) = lenmax
        len(k,iCell) = min(1.0/(1.0 / len1 + 1.0 / len2),lenmax)
        len(k,iCell) = max(len(k,iCell),0.55)
      enddo
      len(nVertLevels+1,iCell) = 1e-15
!      len(1,iCell) = len(2,iCell) !1e-15
    enddo
    print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *, len(:,1)
    print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

  end subroutine build_dissipation_lengths

  subroutine build_sigma_updraft_properties(nCells,nVertLevels)
  !builds the updraft area function

  integer,intent(in) :: nCells,nVertLevels
  integer :: iCell, k
  real :: Sw, w3av, lsigma, wtav, wsav

  do iCell = 1,nCells
    tumd(1,iCell) = 0.0
    wumd(1,iCell) = 0.0
    sigma(1,iCell) = 0.5
    Mc(1,iCell) = 0.0
    do k=2,nVertLevels
      w3av = 0.5*(w3(i2,k-1,iCell) + w3(i2,k,iCell))

      Sw = w3av / (max(w2(i2,k,iCell)**1.5,1e-8))
      lsigma = 0.5 - 0.5*Sw / sqrt(4.0 + Sw**2.0)

      if(lsigma < 0.01) lsigma = 0.01
      if(lsigma > 0.99) lsigma = 0.99

      sigma(k,iCell) = lsigma
      wumd(k,iCell) = sqrt(w2(i2,k,iCell) / (sigma(k,iCell) * (1.0 - sigma(k,iCell))))
      Mc(k,iCell) = sigma(k,iCell)*(1.0 - sigma(k,iCell)) * wumd(k,iCell)
    enddo
  enddo

  end subroutine build_sigma_updraft_properties

  subroutine calc_scalar_updraft_properties(nCells,nVertLevels,wtsfc, wssfc, alphaT, betaS, tlev, adcConst)

    integer,intent(in) :: nCells, nVertLevels, tlev
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, alphaT, betaS
    type(adc_mixing_constants) :: adcConst

    real :: wtav, McAv, sigav, tumdav, wumdav, sumdav, wb, bld, wstar
    integer :: iCell,k

    do iCell=1,nCells
      do k=2,nVertLevels

        tumd(k,iCell) = wt(tlev,k,iCell) / (1.0E-12 + Mc(k,iCell))
        sumd(k,iCell) = ws(tlev,k,iCell) / (1.0E-12 + Mc(k,iCell))
      enddo

      wb = adcConst%grav*(alphaT(iCell)*wtsfc(iCell) - betaS(iCell)*wssfc(iCell))
!      wstar = (abs(0.4*boundaryLayerDepth(iCell)*wb))**(1./3.)

      if(wb > 0.0) then
        wb = adcConst%grav*(alphaT(iCell)*wtsfc(iCell) - betaS(iCell)*wssfc(iCell))
        wstar = (abs(0.4*boundaryLayerDepth(iCell)*wb))**(1./3.)
        w2t(1,iCell) = -0.3*wstar * wtsfc(iCell)
        !Below FIXME!
        w2s(1,iCell) = 0.3*wstar * wssfc(iCell)
      else
        w2t(1,iCell) = 0.0
        w2s(1,iCell) = 0.0
      endif

      !try new boundary condition derived from PDF
      sigav = 0.5*(sigma(1,iCell) + sigma(2,iCell))
      wtav = 0.5*(wt(tlev,1,iCell) + wt(tlev,2,iCell))
      McAv = 0.5*(w2(tlev,1,iCell) + w2(tlev,2,iCell))
      w2t(1,iCell) = (1.0 - 2.0*sigav)*wtav*sqrt(McAv) / (EPSILON + sigav*(1.0-sigav))

      do k=2,nVertLevels
        sigav = 0.5*(sigma(k,iCell) + sigma(k+1,iCell))
        tumdav = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
        sumdav = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
        wumdav = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
        w2t(k,iCell) = sigav*(1.0 - sigav)*(1.0 - 2.0*sigav)*wumdav**2.0*tumdav
        w2s(k,iCell) = sigav*(1.0 - sigav)*(1.0 - 2.0*sigav)*wumdav**2.0*sumdav
      enddo

    enddo
  end subroutine calc_scalar_updraft_properties

  subroutine calc_subplume_fluxes(nCells,nVertLevels,temp,salt,uvel,vvel,BVF,   &
    alphaT,betaS,adcConst,dt)
  ! builds the subplume tendency terms

  integer,intent(in) :: nCells, nVertLevels
  real,dimension(nCells),intent(in) :: alphaT,betaS
  real,intent(in) :: dt
  real,dimension(nVertLevels,nCells),intent(in) :: temp,salt,uvel,vvel
  real,dimension(nVertLevels+1,nCells),intent(in) :: BVF
  type(adc_mixing_constants) :: adcConst

  real :: Uz, Vz, Tz, Sz, B, sigmaAv,integrandTop,integrandBot, Cval
  !calculate length

  integer :: iCell, k

  do iCell = 1,nCells

    lenspsU(1,iCell) = 0.0
    lenspsD(1,iCell) = 0.0
    KmU(1,iCell) = 0.0
    KhU(1,iCell) = 0.0
    KmD(1,iCell) = 0.0
    KhD(1,iCell) = 0.0
    E(1,iCell) = 0.0
    D(1,iCell) = 0.0

    do k = 2,nVertLevels
      !need to add length scales for Up and Down
      Tz = (temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
      Sz = (salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
      Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
      Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))

      if(BVF(k,iCell) <= 0.0) then
        lenspsU(k,iCell) = zmid(k-1,iCell) - zmid(k,iCell)
        lenspsD(k,iCell) = zmid(k-1,iCell) - zmid(k,iCell)
      else
        lenspsU(k,iCell) = min(zmid(k-1,iCell) - zmid(k,iCell),0.76*sqrt(KspsU(k,iCell)/BVF(k,iCell)))
        lenspsD(k,iCell) = min(zmid(k-1,iCell) - zmid(k,iCell),0.76*sqrt(KspsD(k,iCell)/BVF(k,iCell)))
      endif

      KmU(k,iCell) = 0.1*lenspsU(k,iCell)*sqrt( KspsU(k,iCell) )
      KhU(k,iCell) = ( 1.+2.*lenspsU(k,iCell)/( zmid(k-1,iCell) - zmid(k,iCell) ))*KmU(k,iCell)
      wt_spsU(k,iCell) =  -KhU(k,iCell)*Tz
      ws_spsU(k,iCell) =  -KhU(k,iCell)*Sz

      KmD(k,iCell) = 0.1*lenspsD(k,iCell)*sqrt( KspsD(k,iCell) )
      KhD(k,iCell) = ( 1.+2.*lenspsD(k,iCell)/( zmid(k-1,iCell) - zmid(k,iCell) ))*KmD(k,iCell)
      wt_spsD(k,iCell) = -KhD(k,iCell)*Tz
      ws_spsD(k,iCell) = -KhD(k,iCell)*Sz

      E(k,iCell) = adcConst%Cww_E*sigma(k,iCell)*(1.-sigma(k,iCell))*Mc(k,iCell) / ( lendn(k,iCell) + EPSILON )
      D(k,iCell) = adcConst%Cww_D*sigma(k,iCell)*(1.-sigma(k,iCell))*Mc(k,iCell) / ( lenup(k,iCell) + EPSILON )
    enddo

    do k=2,nVertLevels
      eps(k,iCell) = (0.5*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell)))**1.5/len(k,iCell)
     !( sigmaAv*KspsU(k,iCell) + (1.-sigmaAv)*KspsD(k,iCell) )**1.5 / lensps(k,iCell)
      !FIXME we need a ws_spsU part here and shear production!

      if(k==2) then
        Cval = 3.96
      else
        Cval = (0.19+0.51*lenspsU(k,iCell)/(zmid(k-1,iCell) - zmid(k,iCell)))
      endif
      KspsUtend(k,iCell) = adcConst%grav*(alphaT(iCell)*wt_spsU(k,iCell) - betaS(iCell)*ws_spsU(k,ICell)) &
              + ((KmU(k-1,iCell) + KmU(k,iCell))* &
                (KspsU(k-1,iCell) - KspsU(k,iCell)) / (zedge(k-1,iCell) - zedge(k,iCell)) -       &
                (KmU(k,iCell) + KmU(k+1,iCell)) * (KspsU(k,iCell) - KspsU(k+1,iCell)) /           &
                (zedge(k,iCell) - zedge(k+1,iCell))) / (zmid(k-1,iCell) - zmid(k,iCell)) -        &
                Cval*KspsU(k,iCell)**1.5 &
                /lenspsU(k,iCell) + eps(k,iCell) / (2.0*sigma(k,iCell))

      if(k==2) then
        Cval = 3.96
      else
        Cval = (0.19+0.51*lenspsD(k,iCell)/(zmid(k-1,iCell) - zmid(k,iCell)))
      endif

      KspsDtend(k,iCell) = adcConst%grav*(alphaT(iCell)*wt_spsD(k,iCell) - betaS(iCell)*ws_spsD(k,iCell)) &
              + ((KmD(k-1,iCell) + KmD(k,iCell))* &
                (KspsD(k-1,iCell) - KspsD(k,iCell)) / (zedge(k-1,iCell) - zedge(k,iCell)) -       &
                (KmD(k,iCell) + KmD(k+1,iCell)) * (KspsD(k,iCell) - KspsD(k+1,iCell)) /           &
                (zedge(k,iCell) - zedge(k+1,iCell))) / (zmid(k-1,iCell) - zmid(k,iCell)) -        &
                Cval*KspsD(k,iCell)**1.5 / lenspsD(k,iCell) + eps(k,iCell) / (2.0*(1.0 - sigma(k,iCell)))
    enddo
  enddo

  do iCell=1,nCells
    do k=2,nVertLevels
      KspsU(k,iCell) = KspsU(k,iCell) + dt*KspsUtend(k,iCell)
      KspsD(k,iCell) = KspsD(k,iCell) + dt*KspsDtend(k,iCell)
    enddo
  enddo

  end subroutine calc_subplume_fluxes

  subroutine diagnose_momentum_fluxes(nCells,nVertLevels,temp,salt,uvel,vvel,alphaT,betaS,adcConst,dt)
! This routine diagnoses all the horizontal related momentum flux components. All assume steady state
! follows a quasi structure function approach

    integer,intent(in) :: nCells,nVertLevels
    real,dimension(nCells),intent(in) :: alphaT,betaS
    real,dimension(nVertLevels,nCells),intent(in) :: temp,salt,uvel,vvel
    type(adc_mixing_constants) :: adcConst
    real,intent(in) :: dt
    real,dimension(nVertLevels,nCells) :: taupt, taups, taupv
    real :: B, Kps, Kpsp1, diff, lenav, Uz, Vz, Tz, Sz, sigav, sumdav
    real :: tumdav, Ksps
    integer :: iCell, k

    do iCell=1,nCells
      !compute the TOMs first.
      do k=1,nVertLevels
        Ksps = 0.5*((sigma(k,iCell)*KspsU(k,iCell) + (1.0-sigma(k,iCell))*KspsD(k,iCell)) + &
              (sigma(k+1,iCell)*KspsU(k+1,iCell) + (1.0-sigma(k+1,iCell))*KspsD(k+1,iCell)))
        Kps = 0.5*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i2,k,iCell))
        Kpsp1 = 0.5*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i2,k+1,iCell))
        lenav = 0.5*(len(k,iCell) + len(k+1,iCell))
        diff = adcConst%C_mom * sqrt(0.5*(Kps + Kpsp1)) / lenav
        uw2(k,iCell) = -diff*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        vw2(k,iCell) = -diff*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        u2w(k,iCell) = -diff*(u2(i1,k,iCell) - u2(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        v2w(k,iCell) = -diff*(v2(i1,k,iCell) - v2(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        uvw(k,iCell) = -diff*(uv(i1,k,iCell) - uv(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))

        diff = adcConst%C_therm*sqrt(0.5*(Kps + Kpsp1)) / lenav
        uwt(k,iCell) = -diff*(ut(i1,k,iCell) - ut(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        vwt(k,iCell) = -diff*(vt(i1,k,iCell) - vt(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        uws(k,iCell) = -diff*(us(i1,k,iCell) - us(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        vws(k,iCell) = -diff*(vs(i1,k,iCell) - vs(i1,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
      enddo

      do k=2,nVertLevels
        Ksps = sigma(k,iCell)*KspsU(k,iCell) + (1.0 - sigma(k,iCell))*KspsD(k,iCell)
        Kps = sqrt((u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell)))
        B = adcConst%grav*(alphaT(iCell)*wt(i1,k,iCell) - betaS(iCell)*ws(i1,k,iCell))

        taupt(k,iCell) = Kps / (sqrt(2.0)*adcConst%c_pt*len(k,iCell))
        taups(k,iCell) = Kps / (sqrt(2.0)*adcConst%c_ps*len(k,iCell))

        wttend(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) - &
          w2(i1,k,iCell)*(temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + &
          (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*t2(i1,k,iCell) - betaS(iCell)*    &
          ts(i1,k,iCell)) - adcConst%alpha3/4.0*(ut(i1,k,iCell)*Uz + vt(i1,k,iCell)*Vz) +        &
          adcConst%kappa_FL*(wt(i1,k-1,iCell) - wt(i1,k+1,iCell)) / (zedge(k-1,iCell) -          &
          zedge(k+1,iCell))**2.0! - taupt(k,iCell)*wt(i1,k,iCell)

        wttend1(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        wttend2(k,iCell) = -w2(i1,k,iCell)*(temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        wttend3(k,iCell) = (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*t2(i2,k,iCell) - &
            betaS(iCell)* ts(i2,k,iCell))
        wttend4(k,iCell) = - wt(i1,k,iCell) * taupt(k,iCell)
        wttend5(k,iCell) = - adcConst%alpha3/4.0*(ut(i1,k,iCell)*Uz + vt(i1,k,iCell)*Vz)

        wstend(k,iCell) = -(w2s(k-1,iCell) - w2s(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) - &
          w2(i1,k,iCell)*(salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + &
          (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*ts(i2,k,iCell) - betaS(iCell)*    &
          s2(i2,k,iCell)) - adcConst%alpha3/4.0*(us(i1,k,iCell)*Uz + vs(i1,k,iCell)*Vz) +        &
          adcConst%kappa_FL*(ws(i1,k-1,iCell) - ws(i1,k+1,iCell)) / (zedge(k-1,iCell) - &
                zedge(k+1,iCell))**2.0

        wstend1(k,iCell) = -(w2s(k-1,iCell) - w2s(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        wstend2(k,iCell) = -w2(i1,k,iCell)*(salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        wstend3(k,iCell) = (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*ts(i2,k,iCell) - &
                    betaS(iCell)*s2(i2,k,iCell))
        wstend4(k,iCell) = - ws(i1,k,iCell) * taups(k,iCell)
        wstend5(k,iCell) = -adcConst%alpha3/4.0*(us(i1,k,iCell)*Uz + vs(i1,k,iCell)*Vz)

        taupv(k,iCell) = Kps / (adcConst%c_pv*len(k,iCell))
        Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))

        uwtend(k,iCell) = (-(uw2(k-1,iCell) - uw2(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) +   &
                0.5*((0.8-4.0/3.0*adcConst%alpha1)*0.5*Kps**2.0 + (adcConst%alpha1 -                &
                adcConst%alpha2)*u2(i1,k,iCell) + (adcConst%alpha1 + adcConst%alpha2 - 2.0)*        &
                w2(i1,k,iCell))*Uz + 0.5*(adcConst%alpha1 - adcConst%alpha2)*uv(i1,k,iCell)*Vz +    &
                adcConst%beta5*adcConst%grav*(alphaT(iCell)*ut(i1,k,iCell) - betaS(iCell)*          &
                us(i1,k,iCell))) - 2.0*taupv(k,iCell)*uw(i1,k,iCell) + adcConst%kappa_FL*           &
                (uw(i1,k-1,iCell) - uw(i1,k+1,iCell)) / (zedge(k-1,iCell) - zedge(k+1,iCell))**2.0

        uwtend1(k,iCell) = -(uw2(k-1,iCell) - uw2(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        uwtend2(k,iCell) = 0.5*((0.8-4.0/3.0*adcConst%alpha1)*0.5*Kps**2.0 + (adcConst%alpha1 -     &
                adcConst%alpha2)*u2(i1,k,iCell) + (adcConst%alpha1 + adcConst%alpha2 - 2.0)*w2(i1,k,iCell))*Uz
        uwtend3(k,iCell) = 0.5*(adcConst%alpha1 - adcConst%alpha2)*uv(i1,k,iCell)*Vz
        uwtend4(k,iCell) = adcConst%beta5*adcConst%grav*(alphaT(iCell)*ut(i1,k,iCell) - betaS(iCell)*us(i1,k,iCell))
        uwtend5(k,ICell) = - 2.0*taupv(k,iCell)*uw(i1,k,iCell)

        vwtend(k,iCell) = (-(vw2(k-1,iCell) - vw2(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) +  &
          0.5*((0.8-4.0/3.0*adcConst%alpha1)*0.5*Kps**2.0 + (adcConst%alpha1 - adcConst%alpha2)*   &
          v2(i1,k,iCell) + (adcConst%alpha1 - adcConst%alpha2 - 2.0)*w2(i1,k,iCell))*Vz +          &
          0.5*(adcConst%alpha1 - adcConst%alpha2)*uv(i1,k,iCell)*Uz + adcConst%beta5*              &
          adcConst%grav*(alphaT(iCell)*vt(i1,k,iCell) - betaS(iCell)*vs(i1,k,iCell))) -            &
          taupv(k,iCell)*vw(i1,k,iCell) + adcConst%kappa_FL*(vw(i1,k-1,iCell) - vw(i1,k+1,iCell))  &
          / (zedge(k-1,iCell) - zedge(k+1,iCell))**2.0

        uvtend(k,iCell) = (-(uvw(k-1,iCell) - uvw(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) -  &
          (1.0 - 0.5*(adcConst%alpha1+adcConst%alpha2))*(uw(i1,k,iCell)*Vz + vw(i1,k,iCell)*Uz)) - &
           taupv(k,iCell)*uv(i1,k,iCell) + adcConst%kappa_VAR*(uv(i1,k-1,iCell) -                  &
           uv(i1,k+1,iCell)) / (zedge(k-1,iCell) - zedge(k+1,iCell))**2.0

        u2tend(k,iCell) = (-(u2w(k-1,iCell) - u2w(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) +  &
          (1./3.*adcConst%alpha1 + adcConst%alpha2 - 2.0)*uw(i1,k,iCell)*Uz -                      &
          2./3.*adcConst%alpha1*vw(i1,k,iCell)*Vz + 2./3.*(1.-adcConst%beta5)*B - 2./3.*           &
          eps(k,iCell)) + taupv(k,iCell)*(Kps**2/3. - u2(i1,k,iCell)) + adcConst%kappa_VAR*        &
          (u2(i1,k-1,iCell) - u2(i1,k+1,iCell)) / (zedge(k-1,iCell) - zedge(k+1,iCell))**2.0

        u2tend1(k,iCell) = -(u2w(k-1,iCell) - u2w(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        u2tend2(k,iCell) = (1./3.*adcConst%alpha1 + adcConst%alpha2 - 2.0)*uw(i1,k,iCell)*Uz
        u2tend3(k,iCell) = - 2./3.*adcConst%alpha1*vw(i1,k,iCell)*Vz
        u2tend4(k,iCell) = 2./3.*(1.-adcConst%beta5)*B
        u2tend5(k,iCell) = - 2./3.*eps(k,iCell)
        u2tend6(k,ICell) = taupv(k,iCell)*(Kps**2/3. - u2(i1,k,iCell))

        v2tend(k,iCell) = (-(v2w(k-1,iCell) - v2w(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + &
          (1./3.*adcConst%alpha1 +        &
          adcConst%alpha2 - 2.0)*vw(i1,k,iCell)*Vz - 2./3.*adcConst%alpha1*uw(i1,k,iCell)*Uz +  &
          2./3.*(1-adcConst%beta5)*B - 2./3.*eps(k,iCell)) + &
          taupv(k,iCell)*(Kps**2/3. - v2(i1,k,iCell)) + adcConst%kappa_VAR*(v2(i1,k-1,iCell) -  &
          v2(i1,k+1,iCell)) / (zedge(k-1,iCell) - zedge(k+1,iCell))**2.0

        v2tend1(k,iCell) = -(v2w(k-1,iCell) - v2w(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        v2tend2(k,iCell) = (1./3.*adcConst%alpha1 +adcConst%alpha2 - 2.0)*vw(i1,k,iCell)*Vz
        v2tend3(k,iCell) = - 2./3.*adcConst%alpha1*uw(i1,k,iCell)*Uz
        v2tend4(k,iCell) = 2./3.*(1-adcConst%beta5)*B
        v2tend5(k,iCell) = - 2./3.*eps(k,iCell) + taupv(k,iCell)*(Kps**2/3. - v2(i1,k,iCell))

        !taupt = Kps / (2.0*adcConst%c_pt*len(k,iCell))
        !taups = Kps / (2.0*adcConst%c_ps*len(k,iCell))

        Tz = (temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        Sz = (salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))

        uttend(k,iCell) = (-(uwt(k-1,iCell) - uwt(k,iCell))/(zmid(k-1,iCell) - zmid(k,iCell)) - uw(i1,k,iCell)*Tz - &
          (1.0 - adcConst%alpha3)*wt(i1,k,iCell)*Uz) - ut(i1,k,iCell)*taupt(k,iCell)

        vttend(k,iCell) = (-(vwt(k-1,iCell) - vwt(k,iCell))/(zmid(k-1,iCell) - zmid(k,iCell)) - vw(i1,k,iCell)*Tz - &
          (1.0 - adcConst%alpha3)*wt(i1,k,iCell)*Vz) - vt(i1,k,iCell)*taupt(k,iCell)

        ustend(k,iCell) = (-(uws(k-1,iCell) - uws(k,iCell))/(zmid(k-1,iCell) - zmid(k,iCell)) - uw(i1,k,iCell)*Sz - &
          (1.0 - adcConst%alpha3)*ws(i1,k,iCell)*Uz) - us(i1,k,iCell)*taups(k,ICell)

        vstend(k,iCell) = (-(vws(k-1,iCell) - vws(k,iCell))/(zmid(k-1,iCell) - zmid(k,iCell)) - vw(i1,k,iCell)*Sz - &
          (1.0 - adcConst%alpha3)*ws(i1,k,iCell)*Vz) - vs(i1,k,iCell)*taups(k,iCell)

        t2(i2,k,iCell) = tumd(k,iCell)**2.0*sigma(k,iCell)*(1.0-sigma(k,iCell))
        s2(i2,k,iCell) = sumd(k,iCell)**2.0*sigma(k,iCell)*(1.0-sigma(k,iCell))
        ts(i2,k,iCell) = tumd(k,iCell)*sumd(k,iCell)*sigma(k,iCell)*(1.0-sigma(k,iCell))

      enddo
    enddo

    do iCell=1,nCells
      u2cliptend(:,iCell) = 0.0
      v2cliptend(:,iCell) = 0.0
      do k=2,nVertLevels
        u2(i2,k,iCell) = u2(i1,k,iCell) + dt*u2tend(k,iCell)
        if(u2(i2,k,iCell) < 0) then
          u2cliptend(k,iCell) = -u2(i2,k,iCell)
          u2(i2,k,iCell) = 0.0
        endif

        v2(i2,k,iCell) = v2(i1,k,iCell) + dt*v2tend(k,iCell)
        if(v2(i2,k,iCell) < 0) then
          v2cliptend(k,iCell) = -v2(i2,k,iCell)
          v2(i2,k,iCell) = 0.0
        endif

        uw(i2,k,iCell) = uw(i1,k,iCell) + dt*uwtend(k,iCell)
        vw(i2,k,iCell) = vw(i1,k,iCell) + dt*vwtend(k,iCell)
        uv(i2,k,iCell) = uv(i1,k,iCell) + dt*uvtend(k,iCell)
        ut(i2,k,iCell) = ut(i1,k,iCell) + dt*uttend(k,iCell)
        wt(i2,k,iCell) = (wt(i1,k,iCell) + dt*wttend(k,iCell)) / (1.0 + dt*taupt(k,iCell))
        vt(i2,k,iCell) = vt(i1,k,iCell) + dt*vttend(k,iCell)
        us(i2,k,iCell) = us(i1,k,iCell) + dt*ustend(k,iCell)
        vs(i2,k,iCell) = vs(i1,k,iCell) + dt*vstend(k,iCell)
        ws(i2,k,iCell) = (ws(i1,k,iCell) + dt*wstend(k,iCell)) / (1.0 + dt*taups(k,iCell))
        if(abs(wt(i2,k,iCell)) > 1) then
          print *, "ERROR: wt out of range, wt = ",wt(i2,k,iCell)
          print *, "location k,iCell = ", k,iCell
          stop
        endif

        if(abs(ws(i2,k,iCell)) > 1) then
          print *, "ERROR: ws out of range, ws = ",ws(i2,k,iCell)
          print *, "location k,iCell = ", k,iCell
          stop
        endif

      enddo

    enddo

  end subroutine diagnose_momentum_fluxes

  subroutine predict_turbulent_quantities(nCells, nVertLevels, dt, temp, salt, uvel, vvel,alphaT,betaS,adcConst)
    integer,intent(in) :: nCells, nVertLevels
    real,intent(in) :: dt
    real,dimension(nCells),intent(in) :: alphaT, betaS
    type(adc_mixing_constants) :: adcConst
    real,dimension(nVertLevels,nCells) :: temp, salt, uvel, vvel

    real :: Sw, St, Ss, Eav, Dav, sigav, sigavp1, wumdAv, tumdAv, sumdAv, wumdAvp1, tumdAvp1, sumdAvp1
    real :: Swup, KspsUav, KspsDav, KspsUavp1, KspsDavp1, KE, Mcav, lenav,u2av,v2av,w2av
    real :: w3temp, w3check, taups, taupt, mval, KEsps, Uz, Vz

    integer :: iCell, k

    do iCell = 1,nCells
      do k=1,nVertLevels
        Eav = 0.5*(E(k+1,iCell) + E(k,iCell))
        Dav = 0.5*(D(k+1,iCell) + D(k,iCell))
        u2av = 0.5*(u2(i1,k,iCell) + u2(i1,k+1,iCell))
        v2av = 0.5*(v2(i1,k,iCell) + v2(i1,k+1,iCell))
        w2av = 0.5*(w2(i1,k,iCell) + w2(i1,k+1,iCell))

        sigav = 0.5*(sigma(k,iCell) + sigma(k+1,iCell))
        wumdav = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
        tumdav = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
        sumdav = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
        KspsUav = 0.5*(KspsU(k,iCell) + KspsU(k+1,iCell))
        KspsDav = 0.5*(KspsD(k,iCell) + KspsD(k+1,iCell))
        Mcav = 0.5*(Mc(k,iCell) + Mc(k+1,iCell))
        lenav = 0.5*(len(k,iCell) + len(k+1,iCell))
        if(k==nVertLevels) then
          sigavp1 = 0.5*(sigma(k,iCell))
          wumdAvp1 = 0.5*(wumd(k,iCell))
          tumdAvp1 = 0.5*(tumd(k,iCell))
          sumdAvp1 = 0.5*(sumd(k,iCell))
        else
          sigavp1 = 0.5*(sigma(k,iCell) + sigma(k+1,iCell))
          wumdAvp1 = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
          tumdAvp1 = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
          sumdAvp1 = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
        endif

        KEsps = sigav*KspsUav+ (1.0 - sigav)*KspsDav
        KE = sqrt((u2av+v2av+w2av) + 0.0*KEsps)

        !KE = sqrt(sigma(k,iCell)*KspsUav + (1.0 - sigma(k,iCell))*KspsDav)
        Swup = - 2.0/3.0*(KspsU(k,iCell) - KspsU(k+1,iCell)) / (zedge(k,iCell) &
          - zedge(k+1,iCell)) - 2.0/3.0*KspsUav*(log(sigma(k,iCell)) -           &
          log(sigma(k+1,iCell))) / (zedge(k,iCell) - zedge(k+1,iCell)) +         &
          2.0/3.0*(KspsD(k,iCell) - KspsD(k+1,iCell)) / (zedge(k,iCell) -        &
          zedge(k+1,iCell)) + 2.0/3.0*KspsDav*(log(1.0-sigma(k,iCell)) -         &
          log(1.0-sigma(k+1,iCell))) /  (zedge(k,iCell) - zedge(k+1,iCell))

        w3tend(k,iCell) = wumdav**3.0*(Eav*(3.0*sigav - 2.0) + Dav*(3.0*sigav - 1.0)) +    &
                          wumdav**3.0*(6.0*sigav**2.0 - 6.0*sigav + 1)*(sigma(k,iCell)*(1-  &
           sigma(k,iCell))*wumd(k,iCell) -         &
           sigma(k+1,iCell)*(1.0-sigma(k+1,iCell))*wumd(k+1,iCell)) / (zedge(k,iCell) &
            - zedge(k+1,iCell)) - 1.5*sigav*             &
           (1.0 - sigav)*(1.0 - 2.0*sigav)*wumdav**2.0*((1.0 - 2.0*sigma(k,iCell))* &
           wumd(k,iCell)**2.0 -        &
           (1.0 - 2.0*sigma(k+1,iCell))*wumd(k+1,iCell)**2) / (zedge(k,iCell) -   &
           zedge(k+1,iCell)) +       3.0*(1.0 - 2.0*sigav)* &
           Mcav*wumdav*Swup - adcConst%C_mom_w3*KE/(1e-15+sqrt(2.0)*lenAv)*w3(i1,k,iCell) + &
           3.0*adcConst%grav*(alphaT(iCell)*w2t(k,iCell) - betaS(iCell)*w2S(k,iCell))*0.9

        if(k>1 .and. k < nVertLevels) then
           w3tend(k,iCell) = w3tend(k,iCell) + adcConst%kappa_w3*(w3(i1,k-1,iCell) - w3(i1,k+1,iCell)) / (zmid(k-1,iCell) - &
                 zmid(k+1,iCell))**2.0
        endif

        w3tend1(k,ICell) = wumdav**3.0*(Eav*(3.0*sigav - 2.0) + Dav*(3.0*sigav - 1.0))
        w3tend2(k,iCell) = wumdav**3.0*(6.0*sigav**2.0 - 6.0*sigav + 1)*(sigma(k,iCell)*(1-  &
            sigma(k,iCell))*wumd(k,iCell) -         &
            sigma(k+1,iCell)*(1.0-sigma(k+1,iCell))*wumd(k+1,iCell)) / (zedge(k,iCell) &
            - zedge(k+1,iCell))
        w3tend3(k,iCell) = - 1.5*sigav*             &
            (1.0 - sigav)*(1.0 - 2.0*sigav)*wumdav**2.0*((1.0 - 2.0*sigma(k,iCell))* &
            wumd(k,iCell)**2.0 -        &
            (1.0 - 2.0*sigma(k+1,iCell))*wumd(k+1,iCell)**2) / (zedge(k,iCell) -   &
            zedge(k+1,iCell))
       w3tend4(k,ICell) = 3.0*(1.0 - 2.0*sigav)*Mcav*wumdav*Swup- adcConst%C_mom_w3*KE/(1e-15+sqrt(2.0)*lenAv)*w3(i1,k,iCell)
       w3tend5(k,iCell) =  3.0*adcConst%grav*(alphaT(iCell)*w2t(k,iCell) - betaS(iCell)*w2S(k,iCell))*0.9
      enddo

!      k=1
!      w3check = (w2(i1,k,iCell)+w2(i1,k+1,iCell))**1.5
!      w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)
!      do k=2,nVertLevels
!        w3check = (w2(i1,k,iCell) + w2(i1,k+1,iCell))**1.5
!        w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)
!      enddo

      do k=2,nVertLevels-1
        sigav = 0.5*(sigma(k,iCell) + sigma(k-1,iCell))
        wumdav = 0.5*(wumd(k,iCell) + wumd(k-1,iCell))
        tumdav = 0.5*(tumd(k,iCell) + tumd(k-1,iCell))
        sumdav = 0.5*(sumd(k,iCell) + sumd(k-1,iCell))
        KspsUav = 0.5*(KspsU(k,iCell) + KspsU(k-1,iCell))
        KspsDav = 0.5*(KspsD(k,iCell) + KspsD(k-1,iCell))
        Mcav = 0.5*(Mc(k,iCell) + Mc(k-1,iCell))

        sigavp1 = 0.5*(sigma(k,iCell) + sigma(k+1,iCell))
        KspsUavp1 = 0.5*(KspsU(k,iCell) + KspsU(k+1,iCell))
        KspsDavp1 = 0.5*(KspsD(k,iCell) + KspsD(k+1,iCell))

        Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
        Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))

        KEsps = sigma(k,iCell)*KspsU(k,iCell) + (1.0-sigma(k,ICell))*KspsD(k,iCell)
        KE = sqrt((u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell)) + 0.0*KEsps)
        Swup = adcConst%grav*alphaT(iCell)*tumd(k,iCell) - adcConst%grav*        &
          betaS(iCell)*sumd(k,iCell) - 2.0/3.0*(1.0/sigma(k,iCell)*(sigAv*       &
          KspsUav - sigavp1*KspsUavp1) / (zmid(k-1,iCell) - zmid(k,iCell)) -     &
          1.0/(1.0 - sigma(k,iCell))*((1.0 - sigav)*KspsDav - (1.0 - sigavp1)*   &
          KspsDavp1) / (zmid(k-1,iCell) - zmid(k,iCell)))

       w2tend(k,iCell) = -wumd(k,iCell)**2.0*(E(k,iCell) + D(k,iCell)) -         &
!          (Mc(k-1,iCell)*(1.0 - 2.0*sigma(k-1,iCell))*wumd(k-1,iCell)**2 -       &
!          Mc(k+1,iCell)*(1.0 - 2.0*sigma(k+1,iCell))*wumd(k+1,iCell)**2) /       &
          (w3(i1,k-1,iCell) - w3(i1,k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) &
!          (zedge(k-1,iCell) - zedge(k+1,iCell)) + 
          + 2.0*Mc(k,iCell)*Swup - adcConst%C_1* &
          KE / (1.0E-15 + sqrt(2.0)*len(k,iCell))*(w2(i1,k,iCell)-KE**2/3.0) + 4./3.*    &
          adcConst%C_2*sigma(k,iCell)*(1.0-sigma(k,iCell))*wumd(k,iCell)*        &
          (adcConst%grav*alphaT(iCell)*tumd(k,iCell) - adcConst%grav*            &
          betaS(iCell)*sumd(k,iCell)) + (1.0/3.0*adcConst%alpha1 -               &
          adcConst%alpha2)*(uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz) +            &
          adcConst%kappa_VAR*(w2(i1,k-1,iCell) - w2(i1,k+1,iCell)) / (zedge(k-1,iCell) - &
                zedge(k+1,iCell))**2.0

       w2tend1(k,iCell) = -wumd(k,iCell)**2.0*(E(k,iCell) + D(k,iCell))
       w2tend2(k,iCell) = -(w3(i2,k-1,iCell) - w3(i2,k,iCell)) /      &
       (zmid(k-1,iCell) - zmid(k,iCell))
       !-(Mc(k-1,iCell)*(1.0 - 2.0*sigma(k-1,iCell))*wumd(k-1,iCell)**2 -       &
       !   Mc(k+1,iCell)*(1.0 - 2.0*sigma(k+1,iCell))*wumd(k+1,iCell)**2) /       &
       !   (zedge(k-1,iCell) - zedge(k+1,iCell))
        w2tend3(k,ICell) = -adcConst%C_1*KE / (1.0E-15 +   &
            sqrt(2.0)*len(k,iCell))*(w2(i1,k,iCell)-KE**2/3.0)
      w2tend4(k,iCell) = 2.0*Mc(k,iCell)*Swup + 4./3.*adcConst%C_2*sigma(k,iCell)* &
              (1.0-sigma(k,iCell))*wumd(k,iCell)*        &
        (adcConst%grav*alphaT(iCell)*tumd(k,iCell) - adcConst%grav*            &
        betaS(iCell)*sumd(k,iCell))
      w2tend5(k,iCell) =(1.0/3.0*adcConst%alpha1 -               &
      adcConst%alpha2)*(uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz) 

      ! 2.0*Swup*Mc(k,iCell)

!        wttend(k,iCell) = -(E(k,iCell) + D(k,iCell))*wumd(k,iCell)*tumd(k,iCell) - ((1.0-2.0*sigma(k-1,iCell))*       &
!          wumd(k-1,iCell)*Mc(k-1,iCell)*tumd(k-1,iCell) - (1.0-2.0*sigma(k+1,iCell))*wumd(k+1,iCell)*     &
!          tumd(k+1,iCell)*Mc(k+1,iCell)) / (zmid(k-1,iCell) - zmid(k+1,iCell)) - Mcav*wumdav*(temp(k-1,iCell) - &
!          temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + sig(k,iCell)*(1.0-sig(k,iCell))*tumd(k,iCell)*Swup -          &
!          adcConst%C_therm*KE/(len(k,iCell) + 1.0e-15)*wt(i1,k,iCell)
!        wstend(k,iCell) = -(E(k,iCell) + D(k,iCell))*wumd(k,iCell)*sumd(k,iCell) - ((1.0-2.0*sigma(k-1,iCell))*       &
!          wumd(k-1,iCell)*Mc(k-1,iCell)*sumd(k-1,iCell) - (1.0-2.0*sigma(k+1,iCell))*wumd(k+1,iCell)*     &
!          sumd(k+1,iCell)*Mc(k+1,iCell)) / (zmid(k-1,iCell) - zmid(k+1,iCell)) - Mcav*wumdav*(salt(k-1,iCell) - &
!          salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + sigma(k,iCell)*(1.0-sigma(k,iCell))*       &
!          sumd(k,iCell)*Swup - adcConst%C_therm*KE/(len(k,iCell) + 1.0E-15)*ws(i1,k,iCell)

!        taupt = KE / (sqrt(2.0)*adcConst%c_pt*len(k,iCell))
!        taups = KE / (sqrt(2.0)*adcConst%c_ps*len(k,iCell))

!        wttend(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) - &
!          w2(i1,k,iCell)*(temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + &
!          (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*t2(i1,k,iCell) - betaS(iCell)*    &
!          ts(i1,k,iCell)) - wt(i1,k,iCell) * taupt  - adcConst%alpha3/4.0*(ut(i1,k,iCell)*Uz +   &
!          vt(i1,k,iCell)*Vz) + adcConst%kappa_FL*(wt(i1,k-1,iCell) - wt(i1,k+1,iCell)) / (zedge(k-1,iCell) - &
!                zedge(k+1,iCell))**2.0

!        wttend1(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
!        wttend2(k,iCell) = -w2(i1,k,iCell)*(temp(k-1,iCell) - temp(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
!        wttend3(k,iCell) = (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*t2(i2,k,iCell) - &
!            betaS(iCell)* ts(i2,k,iCell))
!        wttend4(k,iCell) = - wt(i1,k,iCell) * taupt
!        wttend5(k,iCell) = - adcConst%alpha3/4.0*(ut(i1,k,iCell)*Uz + vt(i1,k,iCell)*Vz)

!        wstend(k,iCell) = -(w2s(k-1,iCell) - w2s(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) - &
!          w2(i1,k,iCell)*(salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell)) + &
!          (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*ts(i2,k,iCell) - betaS(iCell)*    &
!          s2(i2,k,iCell)) - ws(i1,k,iCell) * taups - adcConst%alpha3/4.0*(us(i1,k,iCell)*Uz +   &
!          vs(i1,k,iCell)*Vz) + adcConst%kappa_FL*(ws(i1,k-1,iCell) - ws(i1,k+1,iCell)) / (zedge(k-1,iCell) - &
!                zedge(k+1,iCell))**2.0

!        wstend1(k,iCell) = -(w2s(k-1,iCell) - w2s(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
!        wstend2(k,iCell) = -w2(i1,k,iCell)*(salt(k-1,iCell) - salt(k,iCell)) / (zmid(k-1,iCell) - zmid(k,iCell))
!        wstend3(k,iCell) = (1.0 - adcConst%gamma1)*adcConst%grav*(alphaT(iCell)*ts(i2,k,iCell) - &
!                    betaS(iCell)*s2(i2,k,iCell))
!        wstend4(k,iCell) = - ws(i1,k,iCell) * taups
!        wstend5(k,iCell) = -adcConst%alpha3/4.0*(us(i1,k,iCell)*Uz + vs(i1,k,iCell)*Vz)  
      !   wstend(k,iCell) = 0.0
      enddo
    enddo

    do iCell=1,nCells
      k=1
      w3check = (w2(i1,k,iCell)+w2(i1,k+1,iCell))**1.5
      w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)
      w2cliptend(iCell,:) = 0.0
      do k=2,nVertLevels
        w2(i2,k,iCell) = w2(i1,k,iCell) + dt*w2tend(k,iCell)
        w3check = (w2(i1,k,iCell) + w2(i1,k+1,iCell))**1.5
        w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)

        if(w2(i2,k,iCell) < 0) then
          w2cliptend(k,iCell) = -w2(i2,k,iCell)
          w2(i2,k,iCell) = 0.0
        endif

        if(abs(w3(i2,k,iCell)) > 1) then
          print *, "ERROR: w3 out of range, w3 = ",w3(i2,k,iCell)
          print *, "location k,iCell = ", k,iCell
          stop
        endif

        if(abs(w2(i2,k,iCell)) > 1) then
          print *, "ERROR: w2 out of range, w2 = ",w2(i2,k,iCell)
          print *, "location k,iCell = ", k,iCell
          stop
        endif

!        wt(i2,k,iCell) = wt(i1,k,iCell) + dt*wttend(k,iCell)
!        ws(i2,k,iCell) = ws(i1,k,iCell) + dt*wstend(k,iCell)
      enddo
    enddo

  end subroutine predict_turbulent_quantities

  subroutine update_mean_fields(dt,nCells,nVertLevels,uvel,vvel,temp,salt,fCor)

    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nVertLevels,nCells),intent(out) :: uvel,vvel,temp,salt
    real,intent(in) :: dt
    real,dimension(nCells),intent(in) :: fCor

    real :: utemp, vtemp
    integer :: iCell,k

    do iCell = 1,nCells
      do k = 1,nVertLevels
        utemp = uvel(k,iCell)
        vtemp = vvel(k,iCell)
        uvel(k,iCell) = uvel(k,iCell) - dt*(uw(i2,k,iCell) - uw(i2,k+1,iCell)) /  &
                  (zedge(k,iCell) - zedge(k+1,iCell)) + dt*fCor(iCell)*vtemp

        vvel(k,iCell) = vvel(k,iCell) - dt*(vw(i2,k,iCell) - vw(i2,k+1,iCell)) /  &
                  (zedge(k,iCell) - zedge(k+1,iCell)) - dt*fCor(iCell)*utemp

        temp(k,iCell) = temp(k,iCell) - dt*(wt(i2,k,iCell) - wt(i2,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
        salt(k,iCell) = salt(k,iCell) - dt*(ws(i2,k,iCell) - ws(i2,k+1,iCell)) / (zedge(k,iCell) - zedge(k+1,iCell))
      enddo
    enddo

  end subroutine update_mean_fields

  subroutine write_turbulent_fields(nCells,nVertLevels,BVF,temp,salt,uvel,vvel)

    integer,intent(in) :: nCells,nVertLevels
    real,dimension(nVertLevels+1,nCells),intent(in),optional :: BVF
    real,dimension(nVertLevels,nCells),intent(in),optional :: temp,salt,uvel,vvel

    integer,parameter :: NDIMS=3
    integer,dimension(NDIMS) :: start,extent
    integer :: ncid,u2_varid,v2_varid,w2_varid,t2_varid,s2_varid, &
      uw_varid,vw_varid,wt_varid,ws_varid,ut_varid,vt_varid,      &
      us_varid,vs_varid,ts_varid,uv_varid,uvw_varid,uw2_varid,    &
      u2w_varid,v2w_varid,vw2_varid,w3_varid,w2t_varid,w2s_varid, &
      sigma_varid,len_varid,ze_varid,zm_varid,kspsu_varid,        &
      kspsd_varid,lspsU_varid, E_varid, D_varid, tumd_varid,       &
      sumd_varid, wumd_varid, MC_varid, bvf_varid, u_varid, v_varid, &
      t_varid, s_varid, l1_varid, l2_varid, l3_varid, lspsD_varid, &
      w2t1_varid,w2t2_varid,w2t3_varid,w2t4_varid, wtt1_varid, &
      wtt2_varid,wtt3_varid,wtt4_varid,w3t1_varid,w3t2_varid, &
      w3t3_varid,w3t4_varid,w3t5_varid, w2t5_varid, u2t1_varid, &
      u2t2_varid,u2t3_varid,u2t4_varid,u2t5_varid,u2t6_varid, &
      uwt1_varid,uwt2_varid,uwt3_varid,uwt4_varid,uwt5_varid, &
      wtt5_varid,u2c_varid,v2c_varid,w2c_varid,v2t1_varid, &
      v2t2_varid,v2t3_varid,v2t4_varid,v2t5_varid,wst1_varid, &
      wst2_varid,wst3_varid,wst4_varid,wst5_varid,lendn_varid, &
      lenup_varid

    integer :: dimids_center(nDIMS),dimids_edges(nDims)
    integer :: ncell_id,nvl_id,ncstat,nvlp1_id,rec_dimid

    integer i,j

    if(defineFirst) then

    ncstat = nf90_create('turb_profiles.nc',NF90_CLOBBER,ncid)

    ncstat = nf90_def_dim(ncid,'nCells',nCells,ncell_id)
    ncstat = nf90_def_dim(ncid,'nVertLevels',nVertLevels,nvl_id)
    ncstat = nf90_def_dim(ncid,'nVertLevelsP1',nvertLevels+1,nvlp1_id)
    ncstat = nf90_def_dim(ncid,'Time',NF90_UNLIMITED,rec_dimid)

    dimids_center = (/nvl_id,ncell_id,rec_dimid/)
    dimids_edges = (/nvlp1_id,ncell_id,rec_dimid/)

    ncstat = nf90_def_var(ncid,'u2',NF90_DOUBLE,dimids_edges,u2_varid)
    ncstat = nf90_def_var(ncid,'v2',NF90_FLOAT,dimids_edges,v2_varid)
    ncstat = nf90_def_var(ncid,'w2',NF90_FLOAT,dimids_edges,w2_varid)
    ncstat = nf90_def_var(ncid,'w2tend1',NF90_FLOAT,dimids_edges,w2t1_varid)
    ncstat = nf90_def_var(ncid,'w2tend2',NF90_FLOAT,dimids_edges,w2t2_varid)
    ncstat = nf90_def_var(ncid,'w2tend3',NF90_FLOAT,dimids_edges,w2t3_varid)
    ncstat = nf90_def_var(ncid,'w2tend4',NF90_FLOAT,dimids_edges,w2t4_varid)
    ncstat = nf90_def_var(ncid,'w2tend5',NF90_FLOAT,dimids_edges,w2t5_varid)
    ncstat = nf90_def_var(ncid,'u2tend1',NF90_FLOAT,dimids_edges,u2t1_varid)
    ncstat = nf90_def_var(ncid,'u2tend2',NF90_FLOAT,dimids_edges,u2t2_varid)
    ncstat = nf90_def_var(ncid,'u2tend3',NF90_FLOAT,dimids_edges,u2t3_varid)
    ncstat = nf90_def_var(ncid,'u2tend4',NF90_FLOAT,dimids_edges,u2t4_varid)
    ncstat = nf90_def_var(ncid,'u2tend5',NF90_FLOAT,dimids_edges,u2t5_varid)
    ncstat = nf90_def_var(ncid,'u2tend6',NF90_FLOAT,dimids_edges,u2t6_varid)
    ncstat = nf90_def_var(ncid,'u2cliptend',NF90_FLOAT,dimids_edges,u2c_varid)
    ncstat = nf90_def_var(ncid,'v2cliptend',NF90_FLOAT,dimids_edges,v2c_varid)
    ncstat = nf90_def_var(ncid,'w2cliptend',NF90_FLOAT,dimids_edges,w2c_varid)
    ncstat = nf90_def_var(ncid,'v2tend1',NF90_FLOAT,dimids_edges,v2t1_varid)
    ncstat = nf90_def_var(ncid,'v2tend2',NF90_FLOAT,dimids_edges,v2t2_varid)
    ncstat = nf90_def_var(ncid,'v2tend3',NF90_FLOAT,dimids_edges,v2t3_varid)
    ncstat = nf90_def_var(ncid,'v2tend4',NF90_FLOAT,dimids_edges,v2t4_varid)
    ncstat = nf90_def_var(ncid,'v2tend5',NF90_FLOAT,dimids_edges,v2t5_varid)

    ncstat = nf90_def_var(ncid,'uwtend1',NF90_FLOAT,dimids_edges,uwt1_varid)
    ncstat = nf90_def_var(ncid,'uwtend2',NF90_FLOAT,dimids_edges,uwt2_varid)
    ncstat = nf90_def_var(ncid,'uwtend3',NF90_FLOAT,dimids_edges,uwt3_varid)
    ncstat = nf90_def_var(ncid,'uwtend4',NF90_FLOAT,dimids_edges,uwt4_varid)
    ncstat = nf90_def_var(ncid,'uwtend5',NF90_FLOAT,dimids_edges,uwt5_varid)
    ncstat = nf90_def_var(ncid,'t2',NF90_FLOAT,dimids_edges,t2_varid)
    ncstat = nf90_def_var(ncid,'s2',NF90_FLOAT,dimids_edges,s2_varid)
    ncstat = nf90_def_var(ncid,'uw',NF90_FLOAT,dimids_edges,uw_varid)
    ncstat = nf90_def_var(ncid,'vw',NF90_FLOAT,dimids_edges,vw_varid)
    ncstat = nf90_def_var(ncid,'wt',NF90_FLOAT,dimids_edges,wt_varid)
    ncstat = nf90_def_var(ncid,'wttend1',NF90_FLOAT,dimids_edges,wtt1_varid)
    ncstat = nf90_def_var(ncid,'wttend2',NF90_FLOAT,dimids_edges,wtt2_varid)
    ncstat = nf90_def_var(ncid,'wttend3',NF90_FLOAT,dimids_edges,wtt3_varid)
    ncstat = nf90_def_var(ncid,'wttend4',NF90_FLOAT,dimids_edges,wtt4_varid)
    ncstat = nf90_def_var(ncid,'wttend5',NF90_FLOAT,dimids_edges,wtt5_varid)
    ncstat = nf90_def_var(ncid,'wstend1',NF90_FLOAT,dimids_edges,wst1_varid)
    ncstat = nf90_def_var(ncid,'wstend2',NF90_FLOAT,dimids_edges,wst2_varid)
    ncstat = nf90_def_var(ncid,'wstend3',NF90_FLOAT,dimids_edges,wst3_varid)
    ncstat = nf90_def_var(ncid,'wstend4',NF90_FLOAT,dimids_edges,wst4_varid)
    ncstat = nf90_def_var(ncid,'wstend5',NF90_FLOAT,dimids_edges,wst5_varid)
    ncstat = nf90_def_var(ncid,'ws',NF90_FLOAT,dimids_edges,ws_varid)
    ncstat = nf90_def_var(ncid,'ut',NF90_FLOAT,dimids_edges,ut_varid)
    ncstat = nf90_def_var(ncid,'vt',NF90_FLOAT,dimids_edges,vt_varid)
    ncstat = nf90_def_var(ncid,'ts',NF90_FLOAT,dimids_edges,ts_varid)
    ncstat = nf90_def_var(ncid,'uv',NF90_FLOAT,dimids_edges,uv_varid)
    ncstat = nf90_def_var(ncid,'us',NF90_FLOAT,dimids_edges,us_varid)
    ncstat = nf90_def_var(ncid,'vs',NF90_FLOAT,dimids_edges,vs_varid)
    ncstat = nf90_def_var(ncid,'E',NF90_FLOAT,dimids_edges,E_varid)
    ncstat = nf90_def_var(ncid,'D',NF90_FLOAT,dimids_edges,D_varid)
    ncstat = nf90_def_var(ncid,'Mc',NF90_FLOAT,dimids_edges,MC_varid)
    ncstat = nf90_def_var(ncid,'tumd',NF90_FLOAT,dimids_edges,tumd_varid)
    ncstat = nf90_def_var(ncid,'sumd',NF90_FLOAT,dimids_edges,sumd_varid)
    ncstat = nf90_def_var(ncid,'wumd',NF90_FLOAT,dimids_edges,wumd_varid)
    ncstat = nf90_def_var(ncid,'len',NF90_FLOAT,dimids_edges,len_varid)
    ncstat = nf90_def_var(ncid,'zedge',NF90_FLOAT,dimids_edges,ze_varid)
    ncstat = nf90_def_var(ncid,'KspsU',NF90_FLOAT,dimids_edges,kspsu_varid)
    ncstat = nf90_def_var(ncid,'KspsD',NF90_FLOAT,dimids_edges,kspsd_varid)
    ncstat = nf90_def_var(ncid,'lenspsD',NF90_FLOAT,dimids_edges,lspsD_varid)
    ncstat = nf90_def_var(ncid,'lenspsU',NF90_FLOAT,dimids_edges,lspsU_varid)
    ncstat = nf90_def_var(ncid,'len1',NF90_FLOAT,dimids_edges,l1_varid)
    ncstat = nf90_def_var(ncid,'len2',NF90_FLOAT,dimids_edges,l2_varid)
    ncstat = nf90_def_var(ncid,'len3',NF90_FLOAT,dimids_edges,l3_varid)
    ncstat = nf90_def_var(ncid,'lendn',NF90_FLOAT,dimids_edges,lendn_varid)
    ncstat = nf90_def_var(ncid,'lenup',NF90_FLOAT,dimids_edges,lenup_varid)

    ncstat = nf90_def_var(ncid,'uvw',NF90_FLOAT,dimids_center,uvw_varid)
    ncstat = nf90_def_var(ncid,'uw2',NF90_FLOAT,dimids_center,uw2_varid)
    ncstat = nf90_def_var(ncid,'u2w',NF90_FLOAT,dimids_center,u2w_varid)
    ncstat = nf90_def_var(ncid,'v2w',NF90_FLOAT,dimids_center,uvw_varid)
    ncstat = nf90_def_var(ncid,'vw2',NF90_FLOAT,dimids_center,vw2_varid)
    ncstat = nf90_def_var(ncid,'w3',NF90_FLOAT,dimids_center,w3_varid)
    ncstat = nf90_def_var(ncid,'w3tend1',NF90_FLOAT,dimids_center,w3t1_varid)
    ncstat = nf90_def_var(ncid,'w3tend2',NF90_FLOAT,dimids_center,w3t2_varid)
    ncstat = nf90_def_var(ncid,'w3tend3',NF90_FLOAT,dimids_center,w3t3_varid)
    ncstat = nf90_def_var(ncid,'w3tend4',NF90_FLOAT,dimids_center,w3t4_varid)
    ncstat = nf90_def_var(ncid,'w3tend5',NF90_FLOAT,dimids_center,w3t5_varid)
    ncstat = nf90_def_var(ncid,'w2t',NF90_FLOAT,dimids_center,w2t_varid)
    ncstat = nf90_def_var(ncid,'w2s',NF90_FLOAT,dimids_center,w2s_varid)
    ncstat = nf90_def_var(ncid,'sigma',NF90_FLOAT,dimids_edges,sigma_varid)
    ncstat = nf90_def_var(ncid,'zmid',NF90_FLOAT,dimids_center,zm_varid)
    ncstat = nf90_def_var(ncid,'TEMP',NF90_FLOAT,dimids_center,t_varid)
    ncstat = nf90_def_var(ncid,'SALT',NF90_FLOAT,dimids_center,s_varid)
    ncstat = nf90_def_var(ncid,'UVEL',NF90_FLOAT,dimids_center,u_varid)
    ncstat = nf90_def_var(ncid,'VVEL',NF90_FLOAT,dimids_center,v_varid)
    ncstat = nf90_def_var(ncid,'BVF',NF90_FLOAT,dimids_edges,bvf_varid)
    ncstat = nf90_enddef(ncid)

    ncstat = nf90_close(ncid)
    record = 1
    defineFirst = .false.
    else

    start(1) = 1
    start(2) = 1
    start(3) = record
    extent(3) = 1
    extent(1) = nVertLevels+1
    extent(2) = nCells

    ncstat = nf90_open('turb_profiles.nc',nf90_write,ncid)
    ncstat = nf90_inq_varid(ncid,'u2',u2_varid)

    ncstat = nf90_put_var(ncid,u2_varid,u2(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2',v2_varid)
    ncstat = nf90_put_var(ncid,v2_varid,v2(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2',w2_varid)
    ncstat = nf90_put_var(ncid,w2_varid,w2(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2tend1',w2t1_varid)
    ncstat = nf90_put_var(ncid,w2t1_varid,w2tend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2tend2',w2t2_varid)
    ncstat = nf90_put_var(ncid,w2t2_varid,w2tend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2tend3',w2t3_varid)
    ncstat = nf90_put_var(ncid,w2t3_varid,w2tend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2tend4',w2t4_varid)
    ncstat = nf90_put_var(ncid,w2t4_varid,w2tend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2tend5',w2t5_varid)
    ncstat = nf90_put_var(ncid,w2t5_varid,w2tend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2tend1',u2t1_varid)
    ncstat = nf90_put_var(ncid,u2t1_varid,u2tend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2tend2',u2t2_varid)
    ncstat = nf90_put_var(ncid,u2t2_varid,u2tend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2tend3',u2t3_varid)
    ncstat = nf90_put_var(ncid,u2t3_varid,u2tend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2tend4',u2t4_varid)
    ncstat = nf90_put_var(ncid,u2t4_varid,u2tend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2tend5',u2t5_varid)
    ncstat = nf90_put_var(ncid,u2t5_varid,u2tend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2tend1',v2t1_varid)
    ncstat = nf90_put_var(ncid,v2t1_varid,v2tend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2tend2',v2t2_varid)
    ncstat = nf90_put_var(ncid,v2t2_varid,v2tend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2tend3',v2t3_varid)
    ncstat = nf90_put_var(ncid,v2t3_varid,v2tend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2tend4',v2t4_varid)
    ncstat = nf90_put_var(ncid,v2t4_varid,v2tend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2tend5',v2t5_varid)
    ncstat = nf90_put_var(ncid,v2t5_varid,v2tend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2tend6',u2t6_varid)
    ncstat = nf90_put_var(ncid,u2t6_varid,u2tend6(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'u2cliptend',u2c_varid)
    ncstat = nf90_put_var(ncid,u2c_varid,u2cliptend(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'v2cliptend',v2c_varid)
    ncstat = nf90_put_var(ncid,v2c_varid,v2cliptend(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2cliptend',w2c_varid)
    ncstat = nf90_put_var(ncid,w2c_varid,w2cliptend(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uwtend1',uwt1_varid)
    ncstat = nf90_put_var(ncid,uwt1_varid,uwtend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uwtend2',uwt2_varid)
    ncstat = nf90_put_var(ncid,uwt2_varid,uwtend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uwtend3',uwt3_varid)
    ncstat = nf90_put_var(ncid,uwt3_varid,uwtend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uwtend4',uwt4_varid)
    ncstat = nf90_put_var(ncid,uwt4_varid,uwtend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uwtend5',uwt5_varid)
    ncstat = nf90_put_var(ncid,uwt5_varid,uwtend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'t2',t2_varid)
    ncstat = nf90_put_var(ncid,t2_varid,t2(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'s2',s2_varid)
    ncstat = nf90_put_var(ncid,s2_varid,s2(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uw',uw_varid)
    ncstat = nf90_put_var(ncid,uw_varid,uw(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'vw',vw_varid)
    ncstat = nf90_put_var(ncid,vw_varid,vw(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wt',wt_varid)
    ncstat = nf90_put_var(ncid,wt_varid,wt(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wttend1',wtt1_varid)
    ncstat = nf90_put_var(ncid,wtt1_varid,wttend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wttend2',wtt2_varid)
    ncstat = nf90_put_var(ncid,wtt2_varid,wttend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wttend3',wtt3_varid)
    ncstat = nf90_put_var(ncid,wtt3_varid,wttend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wttend4',wtt4_varid)
    ncstat = nf90_put_var(ncid,wtt4_varid,wttend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wttend5',wtt5_varid)
    ncstat = nf90_put_var(ncid,wtt5_varid,wttend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wstend1',wst1_varid)
    ncstat = nf90_put_var(ncid,wst1_varid,wstend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wstend2',wst2_varid)
    ncstat = nf90_put_var(ncid,wst2_varid,wstend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wstend3',wst3_varid)
    ncstat = nf90_put_var(ncid,wst3_varid,wstend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wstend4',wst4_varid)
    ncstat = nf90_put_var(ncid,wst4_varid,wstend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wstend5',wst5_varid)
    ncstat = nf90_put_var(ncid,wst5_varid,wstend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'ws',ws_varid)
    ncstat = nf90_put_var(ncid,ws_varid,ws(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'ut',ut_varid)
    ncstat = nf90_put_var(ncid,ut_varid,ut(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'vt',vt_varid)
    ncstat = nf90_put_var(ncid,vt_varid,vt(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'ts',ts_varid)
    ncstat = nf90_put_var(ncid,ts_varid,ts(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'uv',uv_varid)
    ncstat = nf90_put_var(ncid,uv_varid,uv(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'us',us_varid)
    ncstat = nf90_put_var(ncid,us_varid,us(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'vs',vs_varid)
    ncstat = nf90_put_var(ncid,vs_varid,vs(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'lendn',lendn_varid)
    ncstat = nf90_put_var(ncid,lendn_varid,lendn,start,extent)

    ncstat = nf90_inq_varid(ncid,'lenup',lenup_varid)
    ncstat = nf90_put_var(ncid,lenup_varid,lenup,start,extent)

    ncstat = nf90_inq_varid(ncid,'len',len_varid)
    ncstat = nf90_put_var(ncid,len_varid,len,start,extent)

    ncstat = nf90_inq_varid(ncid,'len1',l1_varid)
    ncstat = nf90_put_var(ncid,l1_varid,lenbuoy,start,extent)

    ncstat = nf90_inq_varid(ncid,'len2',l2_varid)
    ncstat = nf90_put_var(ncid,l2_varid,lenblac,start,extent)

    ncstat = nf90_inq_varid(ncid,'len3',l3_varid)
    ncstat = nf90_put_var(ncid,l3_varid,lenshea,start,extent)

    ncstat = nf90_inq_varid(ncid,'zedge',ze_varid)
    ncstat = nf90_put_var(ncid,ze_varid,zedge,start,extent)

    ncstat = nf90_inq_varid(ncid,'KspsU',kspsu_varid)
    ncstat = nf90_put_var(ncid,kspsu_varid,KspsU(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'KspsD',kspsd_varid)
    ncstat = nf90_put_var(ncid,kspsd_varid,KspsD(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'lenspsD',lspsD_varid)
    ncstat = nf90_put_var(ncid,lspsD_varid,lenspsD(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'lenspsU',lspsU_varid)
    ncstat = nf90_put_var(ncid,lspsU_varid,lenspsU(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'E',E_varid)
    ncstat = nf90_put_var(ncid,E_varid,E(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'D',D_varid)
    ncstat = nf90_put_var(ncid,D_varid,D(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'Mc',MC_varid)
    ncstat = nf90_put_var(ncid,MC_varid,Mc(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'tumd',tumd_varid)
    ncstat = nf90_put_var(ncid,tumd_varid,tumd(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'sumd',sumd_varid)
    ncstat = nf90_put_var(ncid,sumd_varid,sumd(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'wumd',wumd_varid)
    ncstat = nf90_put_var(ncid,wumd_varid,wumd(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'sigma',sigma_varid)
    ncstat = nf90_put_var(ncid,sigma_varid,sigma,start,extent)

    if(present(BVF)) then
        ncstat = nf90_inq_varid(ncid,'BVF',bvf_varid)
        ncstat = nf90_put_var(ncid,bvf_varid,BVF,start,extent)
    endif

    extent(1) = nVertLevels

    if(present(temp)) then
      ncstat = nf90_inq_varid(ncid,'TEMP',t_varid)
      ncstat = nf90_put_var(ncid,t_varid,temp,start,extent)
    endif

    if(present(salt)) then
      ncstat = nf90_inq_varid(ncid,'SALT',s_varid)
      ncstat = nf90_put_var(ncid,s_varid,salt,start,extent)
    endif

    if(present(uvel)) then
      ncstat = nf90_inq_varid(ncid,'UVEL',u_varid)
      ncstat = nf90_put_var(ncid,u_varid,uvel,start,extent)
    endif

    if(present(vvel)) then
      ncstat = nf90_inq_varid(ncid,'VVEL',v_varid)
      ncstat = nf90_put_var(ncid,v_varid,vvel,start,extent)
    endif

    ncstat = nf90_inq_varid(ncid,'zmid',zm_varid)
    ncstat = nf90_put_var(ncid,zm_varid,zmid,start,extent)

    ncstat = nf90_inq_varid(ncid,'uvw',uvw_varid)
    ncstat = nf90_put_var(ncid,uvw_varid,uvw,start,extent)

    ncstat = nf90_inq_varid(ncid,'uw2',uw2_varid)
    ncstat = nf90_put_var(ncid,uw2_varid,uw2,start,extent)

    ncstat = nf90_inq_varid(ncid,'u2w',u2w_varid)
    ncstat = nf90_put_var(ncid,u2w_varid,u2w,start,extent)

    ncstat = nf90_inq_varid(ncid,'v2w',v2w_varid)
    ncstat = nf90_put_var(ncid,v2w_varid,v2w,start,extent)

    ncstat = nf90_inq_varid(ncid,'vw2',vw2_varid)
    ncstat = nf90_put_var(ncid,vw2_varid,vw2,start,extent)

    ncstat = nf90_inq_varid(ncid,'w3',w3_varid)
    ncstat = nf90_put_var(ncid,w3_varid,w3(i2,:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w3tend1',w3t1_varid)
    ncstat = nf90_put_var(ncid,w3t1_varid,w3tend1(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w3tend2',w3t2_varid)
    ncstat = nf90_put_var(ncid,w3t2_varid,w3tend2(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w3tend3',w3t3_varid)
    ncstat = nf90_put_var(ncid,w3t3_varid,w3tend3(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w3tend4',w3t4_varid)
    ncstat = nf90_put_var(ncid,w3t4_varid,w3tend4(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w3tend5',w3t5_varid)
    ncstat = nf90_put_var(ncid,w3t5_varid,w3tend5(:,:),start,extent)

    ncstat = nf90_inq_varid(ncid,'w2t',w2t_varid)
    ncstat = nf90_put_var(ncid,w2t_varid,w2t,start,extent)

    ncstat = nf90_inq_varid(ncid,'w2s',w2s_varid)
    ncstat = nf90_put_var(ncid,w2s_varid,w2s,start,extent)

    ncstat = nf90_close(ncid)
    record=record+1
  endif

  if(stopflag) stop
  end subroutine write_turbulent_fields

  subroutine ADC_init(ntimes,nCells,nVertLevels)
    integer,intent(in) :: ntimes,nCells,nVertLevels
    fileTime=0.0
    call init_adc(ntimes,nCells,nVertLevels)
    call write_turbulent_fields(nCells,nVertLevels)
  end subroutine ADC_init

  subroutine ADC_main_loop(nCells,nVertLevels,niter,dt,temp,salt,uvel,vvel,BVF,layerThick,ssh,  &
      uwsfc,vwsfc,wtsfc,wssfc,alphaT,betaS,fCor,fileFrequency,adcConst)

    integer,intent(in) :: nCells,nVertLevels,niter
    real,intent(in) :: dt, fileFrequency
    real,dimension(nVertLevels,nCells),intent(inout) :: temp,salt,uvel,vvel,layerThick
    real,dimension(nCells),intent(in) :: ssh,uwsfc,vwsfc,wtsfc,wssfc,alphaT,betaS,fCor
    real,dimension(nVertLevels+1,nCells),intent(inout) :: BVF
    type(adc_mixing_constants) :: adcConst
    integer :: iIter,iCell,k

    stopflag=.false.
    call construct_depth_coordinate(ssh,layerThick,nCells,nVertLevels)
    do iIter=1,niter
      call build_diagnostic_arrays(nCells,nVertLevels,temp,salt,BVF,wtsfc,wssfc,  &
        uwsfc,vwsfc,alphaT,betaS,adcConst)
      call predict_turbulent_quantities(nCells, nVertLevels, dt, temp, salt, uvel, vvel,  &
        alphaT,betaS,adcConst)
      call diagnose_momentum_fluxes(nCells,nVertLevels,temp,salt,uvel,vvel,alphaT,betaS,adcConst,dt)
      call build_sigma_updraft_properties(nCells, nVertLevels)
      call calc_scalar_updraft_properties(nCells, nVertLevels,wtsfc, wssfc, &
                                alphaT, betaS, i2, adcConst)
      call calc_subplume_fluxes(nCells,nVertLevels,temp,salt,uvel,vvel, BVF,alphaT,betaS,adcConst,dt)
!      call build_dissipation_lengths(nCells,nVertLevels,BVF)
      call dissipation_lengths2(nCells,nVertLevels,temp,salt,alphaT,betaS,zedge)
      call update_mean_fields(dt,nCells,nVertLevels,uvel,vvel,temp,salt,fCor)
      fileTime = fileTime + dt
      if(fileTime >= fileFrequency) then
         call write_turbulent_fields(nCells,nVertLevels,BVF,temp,salt,uvel,vvel)
         fileTime=0.0
      endif
      call swap_time_levels
    enddo
  end subroutine ADC_main_loop

end module adc
