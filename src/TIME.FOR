* Source file TIME.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TmCont(dt,dtMaxW,dtOpt,dMul,dMul2,dtMin,Iter,tPrint,
     !                  tAtm,t,tMax,dtMaxC,ItMin,ItMax,lMinStep,dtInit)
      logical lMinStep
      double precision t,tPrint,tAtm,tMax,tFix

      if(lMinStep) then
        dtMax=amin1(dtMaxW,dtMaxC,dtInit,dtOpt)
        dtOpt=dtMax
        lMinStep=.false.
      else
        dtMax=amin1(dtMaxW,dtMaxC)
      end if
      tFix=dmin1(tPrint,tAtm,tMax)
      if(Iter.le.ItMin.and.(tFix-t).ge.dMul*dtOpt) 
     !  dtOpt=amin1(dtMax,dMul*dtOpt)
      if(Iter.ge.ItMax)
     !  dtOpt=amax1(dtMin,dMul2*dtOpt)
      dt=amin1(dtOpt,sngl(tFix-t))
      iStep=1
      if(dt.gt.0.) iStep=anint(sngl(tFix-t)/dt)
      if(iStep.ge.1.and.iStep.le.10)
     !  dt=amin1(sngl(tFix-t)/iStep,dtMax)
      if(iStep.eq.1) then
        dt=sngl(tFix-t)
        if(dt-dtMax.gt.dtMin) dt=dt/2.
      end if
      if(dt.le.0.0) dt=dtMin/3.

      return
      end

************************************************************************

      double precision function RTime(iMonth,iDay,iHours,iMins,iSecs,
     !                                i100th)

      integer*2 iMonth,iDay,iHours,iMins,iSecs,i100th

      if(iMonth.eq.1.or.iMonth.eq.3.or.iMonth.eq.5.or.iMonth.eq.7.or.
     !   iMonth.eq.8.or.iMonth.eq.10.or.iMonth.eq.12) then
        NoDay=31
      else if(iMonth.eq.4.or.iMonth.eq.6.or.iMonth.eq.9.or.iMonth.eq.11)
     !                                                then
        NoDay=30
      else if(iMonth.eq.2) then
        NoDay=28
      end if
      nMonth=NoDay*24.*60.*60.
      RTime=nMonth+iDay*24.*60.*60.+iHours*60.*60.+iMins*60.+iSecs
      RTime=RTime+i100th/100.

      return
      end

************************************************************************

      subroutine SetBC(tMax,tAtm,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !                 TopInF,BotInF,cT,cBot,NS,tTop,tBot,Ampl,lTemp,
     !                 lChem,KodTop,lVarBC,ierr,lMinStep,lMeteo,Prec,
     !                 rSoil,lLAI,rExtinct,hSink,lCentrif,CosAlf,xConv,
     !                 tConv,iModel,hTopN,iRootIn,xRoot,WLayer,lLinear,
     !                 lActRSU,SPot,iSolAU,rKsTop,nCycles,iCycle,tCycle,
     !                 tInit,jIrrig,crWat,crSal,rLAI,lIsSolPres,lSurf,
     !                 fET,lAPI,lAPI_ET)

      logical TopInF,BotInF,lTemp,lChem,lMinStep,lVarBC,lMeteo,lCentrif,
     !        WLayer,lLAI,lVarGr,lLinear(NS),lActRSU,lIsSolPres,lSurf,
     !        lAPI,lAPI_ET,lDummy,lDummy1,lDummy2,lDummy3
      dimension cBot(NS),cT(NS),cB(11)
      character*3 TheEnd
      double precision tAtm,tMax,tInit,tCycle

      KodTOld=KodTop
      read(31,101) TheEnd
      if(TheEnd.eq.'end') then
        if(iCycle.lt.nCycles) then
          lMinStep=.true.
          if(iCycle.eq.1) tCycle=tAtm-tInit
          iCycle=iCycle+1
          rewind(31)
          iVerA = iGetFileVersion(31,1)
          read(31,*,err=901)
          read(31,*,err=901)
          read(31,*,err=901) !MaxAL
          if(iVerA.eq.4) then
            read(31,*,err=901)
            read(31,*,err=901) lDummy,lDummy,lDummy1,lDummy2,lDummy3
            if(lDummy1) then
              read(31,*,err=901)
              read(31,*,err=901) !rExtinct
            end if
            if(lDummy2) then
              read(31,*,err=901)
              read(31,*,err=901) !nCycles
            end if
            if(lDummy3) then
              read(31,*,err=901)
              read(31,*,err=901) !aInterc
            end if
          end if
          read(31,*,err=901)
          read(31,*,err=901) !hCritS
          read(31,*,err=901)
          read(31,*,err=901)
          crWat=0.
          crSal=0.
        else
          tMax=tAtm
          return
        end if
      end if
      backspace 31
      if(iRootIn.ne.0) then
        if(.not.lChem.and..not.lTemp) then
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT
        else if(lTemp.and..not.lChem) then
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,tBot,
     !                       Ampl
        else
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,tBot,
     !                       Ampl,(cT(jj),cB(jj),jj=1,NS)
        end if
      else if(iRootIn.eq.0) then
        if(.not.lChem.and..not.lTemp) then
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,xRoot
        else if(lTemp.and..not.lChem) then
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,tBot,
     !                       Ampl,xRoot
        else
          read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,tBot,
     !                       Ampl,(cT(jj),cB(jj),jj=1,NS),xRoot
        end if
      end if
      if(iCycle.gt.1) tAtm=(iCycle-1)*tCycle+tAtm
      if(lAPI.and.lAPI_ET) then
        rR=rSoil
        rSoil=0.
      end if
      if(lLAI) then
        rLAI=rR
        rPET=rSoil
        SC=1.
c        SC=hT
        rR=0.
        if(rLAI.gt.0.)
     !    rR=rPET*amax1(0.,1.-exp(-amax1(rExtinct,0.1)*rLAI))*SC
        rSoil=rPET-rR
      end if

      if(lSurf) then
        rSoil=rSoil+(1.-fET)*rR
        rR=fET*rR
      end if

*     Top of the profile
      if(TopInF) then
        lVarGr=.false. ! Variable gravity field - for Scott Jones
        if(lVarGr) then
          if(CosAlf.ne.Prec) lMinStep=.true.
          CosAlf=Prec
          Prec=0.
        end if
        rTopOld=rTop
        hCritA=-abs(hCA)
        if(lVarBC) then
          rTop=Prec
          if(abs(rTopOld-rTop).gt.abs(rTop)*0.2.and.rTop.lt.0.)
     !      lMinStep=.true.
          KodTop=int(rSoil)
          rSoil=0.
          if(KodTop.eq.-1.and.KodTOld.eq.+1.and.Prec.gt.0.
     !      .and.hNewT.gt.0) hNewT=-0.01*xConv
        else
          if(jIrrig.eq.0) then
            if(.not.lMeteo) rTop=abs(rSoil)-abs(Prec)
            if(lMeteo)      rTop=          -abs(Prec)
          end if
          if(abs(rTopOld-rTop).gt.abs(rTop)*0.2.and.rTop.lt.0.)
     !      lMinStep=.true.
          if(rTop.gt.-0.02*rKsTop.and.rTopOld.lt.0..and..not.WLayer)then
            xLimit=0.0
            if(iModel.eq.3) xLimit=-0.03*xConv
            if(KodTop.eq.4.or.hTopN.gt.xLimit) then
              if(iModel.ne.3) xLimit=-0.01*xConv
              hTopN=xLimit
              KodTop=-4
            end if
          end if
        end if
        if(KodTop.eq.3.or.lVarBC) then
          if(abs(hTop-hT).gt.abs(hTop)*0.2) lMinStep=.true.
          hTop=hT
        end if
        rRoot=abs(rR)
      else
        rRoot=abs(rR)
      end if
      hSink=hT ! This should be deleted if Source is not considered

*     Bottom of the profile
      if(BotInF) then
        if(lCentrif) then
          g=9.80665*xConv/tConv/tConv
          CosAlf=hB*hB/g
          hB=0.
          lMinStep=.true.
        end if
        if(abs(rBot-rB).gt.abs(rBot)*0.2) lMinStep=.true.
        rBot=rB
        if(abs(hBot-hB-GWL0L).gt.abs(hBot)*0.2) lMinStep=.true.
        hBot=hB+GWL0L
      end if
      if(lAPI) then
        rBot=rB
        hBot=hB
      end if

      if(lChem) then
        do 11 jj=1,NS
          if(.not.lIsSolPres.and.(cT(jj).gt.0..or.cT(jj).gt.0.))
     !                                          lMinStep=.true.
          if(.not.lLinear(JJ).and.cT(jj).gt.0.) lMinStep=.true.
          cBot(jj)=cB(jj)
          if(lActRSU.and.jj.eq.iSolAU) then
            cBot(jj)=0.
            SPot=cB(jj)
          end if
11      continue
      end if
      return

*     Error when reading from an input file
901   ierr=1
      return

101   format(a3)
      end

************************************************************************

      subroutine SetChemBC(Prec,rSoil,NS,NSD,NumNP,cTop,cT,WLayer,
     !                     hNewT,KodTop,kTopCh,lIsotope,Conc)

      logical WLayer,lIsotope
      dimension cTop(NS),cT(NS),Conc(NSD,NumNP)

      do 11 jj=1,NS
        if(WLayer.and.hNewT.gt.0.) then
          cTop(jj)=cTop(jj)    ! this is handled in the main program
        else
          cTop(jj)=cT(jj)
          if(abs(KodTop).eq.4.and.kTopCh.le.0) then
            if(Prec-rSoil.gt.0.) then
              cTop(jj)=cT(jj)*Prec/(Prec-rSoil)
            else if(rSoil.gt.0.) then
              cTop(jj)=0.
            end if
            if(lIsotope.and.Prec-rSoil.ne.0)
     !        cTop(jj)=(Prec*cT(jj)-rSoil*Conc(jj,NumNP))/(Prec-rSoil)
          end if
        end if
11    continue


      return
      end

************************************************************************

      subroutine DailyVar(tConv,t,rRoot,rRootD)

*     Temperature, max at 1. p.m.
*     Radiation, max at noon, hourly values between 0-6 a.m. and 18-24 p.m.
*     represent 1% of daily value, sinusoid in between

      double precision t

      PI=3.141592654
      tPeriod=1.                 ! one day  !24.*60.*60.*tConv
      tDay=sngl(t)/tConv/86400   ! time in day units

c      if(tPeriod.gt.0.) tTopA=tTop+Ampl*sin(2.*PI*sngl(t)/tPeriod-7.*PI/12.)
      tRemainder=amod(tDay,tPeriod)
      if(tRemainder.le.0.264.or.tRemainder.ge.0.736) then
        rRoot=0.24*rRootD
      else
        rRoot=2.75*rRootD*sin(2.*PI*tDay/tPeriod-6.*PI/12.)
      end if

      return
      end

************************************************************************

      subroutine SinPrec(t,t1,t2,rPrec,rPrecD)

*     Cosinusoidal distribution of precipitation

      double precision t,t1,t2

      PI=3.141592654
      dt=sngl(t2-t1)
      if(rPrecD.gt.0.) then
        rPrec=rPrecD*(1.+1.*cos(2.*PI*sngl(t-t1)/dt-PI))
      else
        rPrec=0.
      end if

      return
      end

************************************************************************

      subroutine Snow(Prec,dt,Temp,SnowMF,SnowSF,SnowLayer,rEvap,xConv,
     !                lMinStep,lChem,cTop,cT,NS,lIsotope)

      logical lMinStep,lChem,lIsotope
      dimension cTop(NS),cT(NS)

      PrecOld=Prec
      rEvapOld=rEvap
      Q=1.
      if(SnowLayer.lt.0.001*xConv) then
        if(Temp.lt.-2.0) then
          Q=1.
        else if(Temp.lt.2.0) then
          Q=1.-((Temp+2.)/4.)
        else
          Q=0.
        end if
      end if

      rTop=Prec*(1.-Q)
      SnowF=Prec*Q

      if(Temp.gt.0..and.SnowLayer.gt.0.) then
        SnowMelt=Temp*SnowMF*dt
      else
        SnowMelt=0.
      end if

      SnowLayer=SnowLayer+SnowF*dt-SnowMelt
      if(SnowLayer.lt.0.) then
        SnowMelt=SnowMelt+SnowLayer
        SnowLayer=0.
      else if(SnowLayer.gt.0..and.rEvap.gt.0.) then
        SnowLayerO=SnowLayer
        if(rEvap*SnowSF*dt.lt.SnowLayer) then
          SnowLayer=SnowLayer-rEvap*SnowSF*dt
          rEvap=0.
        else
          rEvap=(rEvap*SnowSF*dt-SnowLayer)/max(dt,1e-8)
          if(Temp.lt.0.) rEvap=0.
          SnowLayer=0.
        end if
      end if

      Prec=rTop+SnowMelt/max(dt,1e-8)
      if(abs(PrecOld-Prec).gt.abs(Prec)*0.2.and.Prec.gt.0.)
     !      lMinStep=.true.

      if(lChem.and.SnowLayer.gt.0.001*xConv) then
        do 11 jj=1,NS
          if((SnowLayer+dt*(PrecOld-rEvapOld)).gt.0.) then
            if(.not.lIsotope) then
              cTop(jj)=(SnowLayer*cTop(jj)+dt*PrecOld*cT(jj))/
     !                 (SnowLayer+dt*(PrecOld-rEvapOld))
            else
              cTop(jj)=(SnowLayer*cTop(jj)+dt*PrecOld*cT(jj)-
     !         dt*rEvapOld*cTop(jj))/(SnowLayer+dt*(PrecOld-rEvapOld))
            end if
          end if
11      continue
      end if

      return
      end

************************************************************************

*     Irrigation scheduling

      subroutine Irrig(t,dt,tAtmI,rTop,cIrrig,rRoot,hIrrig,jIrrig,
     !                 tIrrig,rIrrig,iIrrigTr,nIrreg,hNode,lMinStep,
     !                 dtMax,dtMin,Evap,Prec)

      double precision t,tAtmI
      logical lMinStep
      dimension tIrrig(4)

      if(t+dtMin.ge.tAtmI) tAtmI=1.d+30
      if(rTop.lt.0.) then
        cIrrig=0.
      else
        cIrrig=cIrrig+rRoot*dt
      end if
      if(hNode.lt.hIrrig.and.jIrrig.eq.0) then    ! Trigger irrigation when pressure drops below some value
c       .and.sngl(t).gt.tIrrig(5)
c       tIrrig(5) - the earliest time irrigation can be triggered (to prevent triggering second
c                   irrigation before moisture front reaches the triggering point)
c       tIrrig(5) is initiated as -1e+30, and updated after each irrigation is triggered as
c       tIrrig(5)=tIrrig(2)+tOffSet (minimum time between irrigations)
        jIrrig=1 ! irrigation is triggered
        nIrreg=nIrreg+1
        if(iIrrigTr.eq.1) tIrrig(3)=(cIrrig+rRoot*tIrrig(4))/rIrrig  ! if calculated from transpiration
        tIrrig(1)=sngl(t)+tIrrig(4)
        tIrrig(2)=sngl(t)+tIrrig(4)+tIrrig(3)
        tAtmI=tIrrig(1)
        write(14,110) nIrreg,t,tIrrig(1),tIrrig(2),tIrrig(3)*rIrrig
      end if                                      ! Irrigate
      if(jIrrig.ge.1) then
        if(dabs(tIrrig(1)-t).lt.0.001*dt) lMinStep=.true.
        if(sngl(t).ge.tIrrig(1).and.sngl(t)+dtMin.lt.tIrrig(2))
     !                                                              then
          jIrrig=2 ! irrigation is active
          rTop=-rIrrig+Evap-Prec   ! to take into acount atmospheric BC
          tAtmI=tIrrig(2)
          if(tIrrig(2)-tIrrig(1).gt.0.) dtMax=(tIrrig(2)-tIrrig(1))/5.
        end if
        if(t+dtMin.ge.tIrrig(2)) then
          jIrrig=0
          rTop=Evap-Prec   ! to take into acount atmospheric BC
          cIrrig=0.
          dtMax=1.e+30
        end if
      end if
      return

110   format(//' Irrigation event ',i3/
     !         ' ================ '/
     !         ' Time when irrigation is triggered: ',f12.5/
     !         ' Beginning of irrigation:           ',f12.5/
     !         ' End of irrigation:                 ',f12.5/
     !         ' Irrigation amount:                 ',f12.5)
      end

************************************************************************

      subroutine ICARDAInit(iCarda,nICarda,iICarda,cDataPath)

      character cFileName*260,cDataPath*260

      iCarda =0   ! =1: ICarda project
      nICarda=1   ! Number of iCarda cycles
      iICarda=1   ! Initialization of iCarda cycles
      if(iCarda.eq.1) then
        write(*,*) "enter the number of computational cycles"
        read(*,*) nICarda
        iLengthPath = Len_Trim(cDataPath)
        cFileName = cDataPath(1:iLengthPath)//'\ICarda.Out'
        open(21,file=cFileName,status='unknown',err=901)
        write(21,110)
      end if

901   return

110   format(//
     !' Cycle    UpperFlux   BottomFlux   PotentialT      ActualT LeachF
     !r   Yield    RootZoneH RootZoneConc   BottomConc  WatRed  SalRed')
      end

************************************************************************

      subroutine ICARDASet(t,tInit,tOld,dt,dtInit,nICarda,iICarda,
     !                     TLevel,ALevel,PLevel,tAtm,tAtm1,rTop,Prec,
     !                     rSoil,rRoot,cTop,cBot,CumQ,ConcB,hRoot,cRoot,
     !                     crWat,crSal)

      double precision t,tInit,tOld,tAtm,tAtm1
      integer PLevel,Alevel,TLevel
      logical lDummy,lDummy1
      dimension CumQ(17)

      rLeach=1.
      if(CumQ(7).gt.0.) rLeach=-CumQ(5)/CumQ(7)
      Yield=0.
      if(CumQ(2).gt.0.) Yield=CumQ(4)/CumQ(2)
      if(CumQ(2).gt.0.) crWat=crWat/CumQ(2)
      if(CumQ(2).gt.0.) crSal=crSal/CumQ(2)
      write(21,110) iICarda,CumQ(7),-CumQ(5),CumQ(2),CumQ(4),rLeach,
     !              Yield,hRoot,cRoot,ConcB,crWat,crSal
      do 11 i=1,12
        CumQ(i)=0.
11    continue
      iICarda=iIcarda+1
      if(iICarda.gt.nICarda) return

      tInit=0.
      tOld=tInit
      dt=dtInit
      t=tInit+dt
      TLevel=2
      ALevel=1
      PLevel=1
      crWat=0.
      crSal=0.

      rewind(31)
      iVerA = iGetFileVersion(31,1)
      read(31,*,err=901)
      read(31,*,err=901)
      read(31,*,err=901) !MaxAL
      if(iVerA.eq.4) then
        read(31,*,err=901)
        read(31,*,err=901) lDummy,lDummy,lDummy1
        if(lDummy1) then
          read(31,*,err=901)
          read(31,*,err=901) !rExtinct
        end if
      end if
      read(31,*,err=901) 
      read(31,*,err=901) !hCritS
      read(31,*,err=901)
901   continue
      read(31,*,err=901) tAtm,Prec,rSoil,rRoot,hCA,rB,hB,hT,tTop,
     !                   tBot,Ampl,cTop,cBot
      rTop=abs(rSoil)-abs(Prec)
      tAtm1=tAtm

      return

110   format(i6,4e13.5,2f8.3,3e13.5,2f8.4)
      end

************************************************************************

      subroutine Root2(rRoot,rRoot2)

      logical lOpen

      inquire(unit=35,opened=lOpen)
      if(lOpen) then
        read(35,*,err=901) xx,Frac
        if(Frac.gt.0.) then
          rRoot2=rRoot*Frac
          rRoot=rRoot*(1.-Frac)
        end if
      end if

901   continue
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
