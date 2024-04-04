*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     HYDRUS   - Numerical model of one-dimensional variably saturated *
*                water flow, heat transport, and transport of solutes  *
*                involved in sequential first-order decay reactions    *
*                                                                      *
*                Version coupled with user-friendly interface for      *
*                MS Windows environment                                *
*                                                                      *
*     Designed by J.Simunek and M. Th. van Genuchten (1996)            *
*                                                                      *
*                                           Last modified: July, 2013  *
*                                                                      *
*           COPYRIGHT (c) 2005, Jiri Simunek, PC-Progress              *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      program HYDRUS

*#######################################################################
*     Procedura pro osetreni nasledujicich signalu: 
*         - signal CTRL+C            
*         - signal CTRK+BREAK   
*         - abnormal termination 
*         - floating point error 
*-----------------------------------------------------------------------

      interface
        function h_sig (signum)
          !MS$ATTRIBUTES C :: h_sig
          integer(4) h_sig
          integer(2) signum
        end function
      end interface

      interface
        function hand_fpe (sigid, except)
           !MS$ATTRIBUTES C :: hand_fpe
          integer(4) hand_fpe
          integer(2)  sigid, except
        end function
      end interface
*#######################################################################

      parameter (NumNPD=5001,
     !           NMatD =100,
     !           NTabD =100,
     !           NObsD =100,
     !           NSD   =11,
     !           NUnitD=7,
     !           NPD   =1000)

      integer PLevel,Alevel,TLevel,PTLevel,err
      logical lSink,WLayer,qGWLF,TopInF,ShortO,lWat,lChem,lTemp,ConvgF,
     !        FreeD,SeepF,BotInF,AtmBC,lRoot,lUpW,lWTDep,lSTDep,lEquil,
     !        lLinear(NSD),lArtD,lMoSink,lSolRed,lSolAdd,lScreen,lTable,
     !        lMinStep,qDrain,lMsSink,lTort,lInitW,lVarBC,lPrint,lFiltr,
     !        lMobIm(NMatD),lMeteo,lBact,lVapor,lDayVar,lEnter,lEnBal,
     !        lVaporOut,lExtrap,lPrintD,lLAI,lDensity,lSnow,lSolHyst,
     !        lCentrif,lSource,lEqInit,lSolub,lSinPrec,lIrrig,lThAds,
     !        lDualNEq,lFlux,lHargr,lMetDaily,lEnd,lMassIni,lSurfact,
     !        lStopConv,lActRSU,lOmegaW,lFluxOut,lIsotope,lFumig,
     !        lAddFum,lHanks,lDblVeg,lYield,lIsSolPres,lPrHeat,lSurf,
     !        lTillage,lAPI,lAPI_ET,lDummy,lOust,lSchwen,lYong,lCond,
     !        lRootIn,lOpen
      double precision P,R,S,Q,RTime1,t,tInit,tOld,RTime,tMax,tAtm,
     !                 tPrint1,tAtm1,tAtm2,tAtmOld,tAtmN,tAtm2O,tAtmI,
     !                 TPrint,tPrintInt,tAddFum,tCycle
      character cFileName*260,cDataPath*260
      integer*2 status
      integer*2 i2,iYear,iMonth,iDay,iHours,iMins,iSecs,i100th
*---- START gettim and getdat gnu alternative function VARABLES ----
      integer dtv(8)
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer hundred
*---- START gettim and getdat gnu alternative function VARABLES ----


      dimension x(NumNPD),hNew(NumNPD),hOld(NumNPD),hTemp(NumNPD),
     !  MatNum(NumNPD),ParD(11,NMatD),TPrint(NPD),Sink(NumNPD),
     !  POptm(NMatD),Beta(NumNPD),LayNum(NumNPD),CumQ(17),SubVol(10),
     !  ConTab(NTabD,NMatD),CapTab(NTabD,NMatD),TheTab(NTabD,NMatD),
     !  Con(NumNPD),Cap(NumNPD),ConSat(NMatD),ths(NMatD),P(NumNPD),
     !  R(NumNPD),S(NumNPD),hSat(NMatD),ParW(11,NMatD),Area(10),
     !  ChPar(NSD*16+4,NMatD),Conc(NSD,NumNPD),vNew(NumNPD),
     !  vOld(NumNPD),Disp(NumNPD),Retard(NumNPD),g0(NumNPD),g1(NumNPD),
     !  Q(NumNPD),wc(NumNPD),ThNew(NumNPD),ThOld(NumNPD),WatIn(NumNPD),
     !  SolIn(NumNPD),Node(NObsD),iUnit(NUnitD),cvCh0(NSD),cvCh1(NSD),
     !  cTop(NSD),cBot(NSD),cRoot(NSD),CumCh(11,NSD),cVOlI(NSD),
     !  cCumA(NSD),cCumT(NSD),cvTop(NSD),cvBot(NSD),cNew(NumNPD),
     !  cTemp(NumNPD),cPrevO(NumNPD),TDep(NSD*16+4),TPar(10,NMatD),
     !  TempO(NumNPD),TempN(NumNPD),SorbN(NumNPD),Sorb(NSD,NumNPD),
     !  q0(NumNPD),q1(NumNPD),aOsm(NSD),thr(NMatD),sSink(NumNPD),
     !  Ah(NumNPD),AK(NumNPD),ATh(NumNPD),AhW(NMatD),AThW(NMatD),
     !  AKW(NMatD),Kappa(NumNPD),AThS(NumNPD),ThRR(NumNPD),cvChR(NSD),
     !  ConR(NumNPD),ConO(NumNPD),AKS(NumNPD),KappaO(NumNPD),
     !  cRootMax(NSD),cvChIm(NSD),hTab(NTabD,NMatD),NTab(NMatD),
     !  rGrowth(1000,5),ConLT(NumNPD),ConVT(NumNPD),ConVh(NumNPD),
     !  ThEq(NumNPD),ThVOld(NumNPD),ThVNew(NumNPD),vVOld(NumNPD),
     !  vVNew(NumNPD),Sorb2(NSD,NumNPD),SorbN2(NumNPD),tIrrig(4),
     !  ThNewIm(NumNPD),ThOldIm(NumNPD),SinkIm(NumNPD),STrans(NumNPD),
     !  DMoist(NMatD,NSD,13,6),Sorb3(NSD,NumNPD),WDep(2+NMatD,NSD*9),
     !  Beta2(NumNPD),Sink2(NumNPD),cT(NSD),iTill(NPD),CumChPrev(2),
     !  SinkDr(NumNPD),cvChDr(NSD)

*     Couvreur model
      logical lCouvreur
      integer StresFun,BCtp
      double precision Krs,Kcomp

      data lStopConv,iNonConv /.true. , 0 /

*-----------------------------------------------------------------------
*      i4ret = SIGNALQQ(SIG$ABORT, h_sig)
*      i4ret = SIGNALQQ(SIG$BREAK, h_sig)
*      i4ret = SIGNALQQ(SIG$INT,   h_sig)
*      i4ret = SIGNALQQ(SIG$FPE,   hand_fpe)
*-----------------------------------------------------------------------

      lOust=.false.     ! DuPont (erosion)
      lTillage=.false.  ! Tillage mixing
      lSchwen=.false.   ! Time-variable scaling factors
      lYong=.false.     ! Time-variabld BC with hCritS and WLayer
      lRootIn=.false.   ! Read parameters for spatial root distribution

*      iCount = NARGS()  ! TODO: changes were here
      iCount = IARGC()
      if(iCount.gt.1) then
        i2=1
        write (*,*) 'unhandled scenario'
*        call GETARG(i2, cDataPath, status)
      else
        cFileName='LEVEL_01.DIR'
        open(10,file=cFileName, status='old',err=901)
        read(10,101,err=904) cDataPath
        close(10)
      end if
      iLengthPath = Len_Trim(cDataPath)
      if(iLengthPath.gt.260-13) goto 930
      cFileName = cDataPath(1:iLengthPath)//'/OPTIONS.IN'
      open(35,file=cFileName, status='old',err=998)
998   continue

*     Initialization
      call Init(CosAlf,NTab,ItCum,TLevel,ALevel,PLevel,hRoot,vRoot,
     !          IterW,IterC,dtMaxC,wCumT,wCumA,err,lVarBC,NSD,cRoot,
     !          cCumT,cCumA,NumNPD,Sink,wc,CumQ,lMeteo,lBact,lVapor,
     !          lEnBal,lDayVar,lEnter,lFiltr,TauW,nPrStep,nTabMod,
     !          iDualPor,dtMaxT,lExtrap,lPrintD,lLAI,rExtinct,lDensity,
     !          ExcesInt,lMinstep,lPrint,lSnow,SnowMF,SnowSF,SnowLayer,
     !          lSolHyst,cTemp,iSunSh,iRelHum,xRoot,lCentrif,Radius,
     !          lSource,iNode,hSink,CumSource,GWL0L,tPrintInt,iTort,
     !          iEnhanc,hSeep,lEqInit,lSolub,lSinPrec,lIrrig,iIrrig,
     !          jIrrig,hIrrig,rIrrig,tIrrig,lThAds,iMoistDep,Sorb3,
     !          OmegaC,WTransf,lDualNEq,lFlux,iCrop,iRootIn,lMassIni,
     !          lMetDaily,iDummy,iDummy,lDummy,Sorb2,lActRSU,OmegaS,
     !          SPot,lOmegaW,OmegaW,lSurfact,lEnd,lVaporOut,lFluxOut,
     !          lIsotope,lFumig,lHanks,nCycles,iCycle,tCycle,iSolAU,
     !          lDblVeg,Sink2,lYield,Yield,PTLevel,crWat,crSal,iInterc,
     !          aInterc,lIsSolPres,lPrHeat,CumHeatFlB,CumHeatFlT,iDrKod,
     !          SinkDr,WDrain,lSurf,fET,lCond,lAPI,lAPI_ET,0)
      data iUnit /50,70,71,75,76,77,78/

      cFileName = cDataPath(1:iLengthPath)//'/SELECTOR.IN'
      open(30,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'/PROFILE.DAT'
      open(32,file=cFileName, status='old',err=901)
      cFileName = cDataPath(1:iLengthPath)//'/I_CHECK.OUT'
      open(50,file=cFileName, status='unknown',err=902)
      if(lPrint) then
        cFileName = cDataPath(1:iLengthPath)//'/RUN_INF.OUT'
        open(70,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/T_LEVEL.OUT'
        open(71,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/NOD_INF.OUT'
        open(75,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/BALANCE.OUT'
        open(76,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/OBS_NODE.OUT'
        open(77,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/PROFILE.OUT'
        open(78,file=cFileName, status='unknown',err=902)
      end if

*     Read input data --------------------------------------------------

      call BasInf(CosAlf,MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,lWat,
     !            lChem,lSink,WLayer,qGWLF,FreeD,SeepF,AtmBC,KodTop,
     !            KodBot,rTop,rRoot,rBot,hCritS,hCritA,GWL0L,Aqh,Bqh,
     !            kTOld,kBOld,NUnitD,iUnit,NMat,NMatD,NLay,lRoot,lTemp,
     !            lWTDep,lEquil,lScreen,qDrain,zBotDr,BaseGW,rSpacing,
     !            iPosDr,rKhTop,rKhBot,rKvTop,rKvBot,Entres,WetPer,
     !            zInTF,GeoFac,lInitW,lVarBC,xConv,tConv,lMeteo,lVapor,
     !            iVer,lPrint,lCentrif,lSnow,hSeep,lFlux,lActRSU,lIrrig,
     !            iDrKod,err)
      if(err.ne.0) goto (905,906,903,916) err
      if(.not.lChem) lFumig=.false.

      call ICARDAInit(iCarda,nICarda,iICarda,cDataPath)
      if(lAPI.and.lChem) then
        cFileName = cDataPath(1:iLengthPath)//'/API_OBS.OUT'
        open(41,file=cFileName, status='unknown',err=902)
        cFileName = cDataPath(1:iLengthPath)//'/API_NOD.OUT'
        open(42,file=cFileName, status='unknown',err=902)
      end if

      if(lPrint) then
        if(lFluxOut) then
          cFileName = cDataPath(1:iLengthPath)//'/T_LEVEL1.OUT'
          open(44,file=cFileName, status='unknown',err=902)
        end if
        if(lVaporOut.and.lVapor) then
          cFileName = cDataPath(1:iLengthPath)//'/NOD_INF_V.OUT'
          open(45,file=cFileName, status='unknown',err=902)
        end if
      end if
      if(TopInF.or.BotInF.or.AtmBC) then
        if(.not.lRoot) lSurf=.false.
        cFileName = cDataPath(1:iLengthPath)//'/ATMOSPH.IN'
        open(31,file=cFileName, status='old',err=901)
        if(lPrint) then
          cFileName = cDataPath(1:iLengthPath)//'/A_LEVEL.OUT'
          open(72,file=cFileName, status='unknown',err=902)
        end if
        if(lMeteo) then
          cFileName = cDataPath(1:iLengthPath)//'/METEO.IN'
          open(33,file=cFileName, status='old',err=901)
          if(lPrint) then
            cFileName = cDataPath(1:iLengthPath)//'/METEO.OUT'
            open(43,file=cFileName, status='unknown',err=902)
          end if
        end if
      end if

      call NodInf(NumNPD,NumNP,hTop,hBot,x,hNew,hOld,MatNum,hTemp,
     !            LayNum,Beta,Ah,AK,ATh,Conc,Sorb,TempN,TempO,NSD,NS,
     !            xSurf,lChem,lTemp,lEquil,lScreen,lBact,Sorb2,lSolub,
     !            err,lAPI,cRootIni,lDualNEq,lDblVeg,Beta2)
      if(err.ne.0) goto (912,903,914,929) err
      call ObsNodIn(NumNP,NObsD,NObs,lPrint,lChem,Node,NS,lFlux,lAPI,
     !              err)
      if(err.ne.0) goto (912,922,915) err
      if(lIrrig) then
        cFileName = cDataPath(1:iLengthPath)//'/IRRIG.OUT'
        open(14,file=cFileName,status='unknown',err=901)
        call IrrigIn(iIrrig,hIrrig,rIrrig,tIrrig,iIrrigTr,nIrreg,err)
        if(err.ne.0) goto 906
      end if
      call MatIn (NMat,ParD,ParW,hTab(1,1),hTab(NTab(1),1),lScreen,err,
     !            NumNP,Ah,iHyst,AhW,AThW,AKW,MatNum,hNew,Kappa,AThS,
     !            ThRR,ConR,AKS,KappaO,iModel,xConv,lTable,IKappa,
     !            nTabMod,iDualPor,0,0,.false.,.false.)
      if(err.ne.0) goto (906,903) err
      if(iHyst.gt.0.and.iHyst.ne.3)
     !call HysterIn(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,Kappa,
     !              AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,iHyst,iModel,
     !              cDataPath)
      if(iModel.eq.nTabMod) then
        cFileName = cDataPath(1:iLengthPath)//'/MATER.IN'
        open(36,file=cFileName, status='old',err=901)
      end if
      if(lInitW.and.iModel.lt.nTabMod) then
        call InitW(NumNP,NMat,Matnum,Kappa,hNew,hOld,hTemp,ParD,ParW,
     !             iModel,hTop,hBot,iDualPor,ThNewIm,err)
        if(err.ne.0) goto 924
      end if

      if(lCond) then   ! Decrease of conductivity with time in M=1
        iCond=0
        if(iCond.eq.0) then
          cA=ParD(6,1)
          if(abs(cA-0.5).lt.0.001) lCond=.false.
          ParD(6,1)=0.5
        else if(iCond.eq.1.or.iCond.eq.2) then
          inquire(unit=35,opened=lOpen)
          if(lOpen) then
            iCond=0
            read(35,*,err=322)
            read(35,*,err=322) iCond
            read(35,*,err=322)
            read(35,*,err=322) t0,cA,cB,cC
            goto 323
322         lCond=.false.
323         continue
          end if
        end if
      end if

      call GenMat(NTab,NTabD,NMat,thr,ths,hSat,ParD,hTab,ConTab,CapTab,
     !            ConSat,TheTab,iModel,lScreen,nTabMod,ConSMax,xConv,
     !            tConv,err)
      if(err.ne.0) goto (903,937) err
      if(iHyst.ne.3) then
       call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,MatNum,
     !            ParD,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,TheTab,ThOld,
     !            thr,ths,lWTDep,TempN,IterW,ConO,Kappa,AThS,ThRR,ConR,
     !            AKS,AhW,AThW,AKW,iModel,lTable,lVapor,ThVOld,ConLT,
     !            ConVT,ConVh,xConv,tConv,nTabMod,hCritA,lDensity,Conc,
     !            NSD,iEnhanc)
      else ! Bob Lenhard hysteresis
        call Hyst(NumNP,NMatD,ParD,ParW,MatNum,Kappa,hNew,hOld,ThOld,
     !            Con,Cap,IKappa,1)
      end if
      call InitDualPor(NumNP,NMat,Matnum,ParD,ThOld,iDualPor,ThNewIm,
     !                 ThOldIm,SinkIm,hNew,STrans,lInitW)
      call TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,TPrint,
     !          t,dtOpt,TopInF,BotInF,lScreen,ItMin,ItMax,MaxAL,hCritS,
     !          NPD,AtmBC,iVer,lPrintD,nPrStep,tPrintInt,lEnter,lDayVar,
     !          lSinPrec,lLAI,rExtinct,nCycles,lChem,lYield,CTp,CRp,
     !          iInterc,aInterc,err)
      if(err.ne.0) goto (907,928) err
      tPrint1=tMax
      if(lPrintD) TPrint1=tInit+tPrintInt
      if(lMeteo) then
        call MeteoIn(rLat,rAlt,ShWRadA,ShWRadB,rLWRadA,rLWRadB,rLWRadA1,
     !               rLWRadB1,WindHeight,TempHeight,iCrop,iLAI,
     !               CropHeight,Albedo,xLAI,xRoot,iInterc,aInterc,
     !               nGrowth,rGrowth,rExtinct,iRadiation,lEnBal,lPrint,
     !               iSunSh,iRelHum,CloudF_Ac,CloudF_Bc,lHargr,
     !               lMetDaily,xConv,err)
        if(err.ne.0) goto (932) err
        if(lEnBal) iCrop=0
      end if
      tAtm2=tMax
      dtInit=dt

      if(lRoot) then
        call RootIn(tRMin,tRHarv,xRMin,xRMax,RGR,lScreen,iver,iRootIn,
     !              nGrowth,rGrowth,tRPeriod,err)
        if(err.ne.0) goto (908,903) err
      end if
      if(lRootIn) call BetaIn(NumNP,Beta,x,cDataPath,lScreen)
      if(lRoot.or.iCrop.gt.1)
     !  call SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR,xRoot,
     !             lRoot,iRootIn,nGrowth,rGrowth,tRPeriod,fET)
      if(lTemp) then
        call TempIn(NMat,TPar,Ampl,tPeriod,kTopT,tTop,kBotT,tBot,TopInF,
     !              BotInF,iCampbell,iVer,lSnow,SnowMF,SnowSF,SnowLayer,
     !              lScreen,err,0)
        if(err.ne.0) goto (910,903) err
        if(lEnBal) kTopT=-1
        if(lPrHeat) then
          cFileName = cDataPath(1:iLengthPath)//'/HEATFL.OUT'
          open(43,file=cFileName, status='unknown',err=902)
        end if
      end if

      if(lChem) then
        call ChemIn(lUpW,lSTDep,NMat,NS,NSD,MaxItC,ChPar,TDep,kTopCh,
     !              cTop,kBotCh,cBot,epsi,tPulse,CumCh,cTolA,cTolR,
     !              lLinear,lEquil,lArtD,PeCr,lScreen,dSurf,cAtm,lTort,
     !              lMobIm,lBact,lFiltr,iMoistDep,WDep,NMatD,iModel,
     !              ParD,iVer,lDualNEq,lMassIni,lEqInit,iTort,err,0)
        if(err.ne.0) goto (911,903,931) err
        if(iDualPor.gt.0) lEquil=.false.
        if(iMoistDep.eq.2)
     !    call MoistDepIn(cDataPath,cFileName,NMat,NMatD,NS,NSD,DMoist,
     !                    iMoistDep)
        if(lFumig) then
          call FumigIn(cDataPath,lFumig,lAddFum,tAddFum,FumMass,zAddFT,
     !                 zAddFB,tInit,tPrev,CumChPrev,lPrint,NS,err)
          if(err.eq.1) goto 902
        end if
        if(lMassIni)
     !    call MassInit(NumNP,NSD,NS,NMat,MatNum,TDep,TempN,ChPar,Conc,
     !                  ThOld,ThOldIm,ths,lLinear,lBact,lFumig)
        if(.not.lEquil.and.(lEqInit.or.lMassIni))
     !    call NonEqInit(NumNP,NSD,NS,NMat,MatNum,TDep,TempO,ChPar,
     !                   Conc,Sorb,lLinear,lMobIm,iDualPor,lBact,Sorb2,
     !                   ThOld,lFumig)
        if(lPrint) then
          call OpenSoluteFiles(NS,cDataPath,iLengthPath,cFileName,err)
          if(err.eq.1) goto 902
        end if
        call IsSolutePresent(lIsSolPres,NSD,NS,NumNP,Conc)
      end if
      if(lSink) then
        call SinkIn(NMat,lChem,lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,
     !              P2L,P3,r2H,r2L,aOsm,c50,P3c,NS,lMsSink,cRootMax,
     !              iVer,OmegaC,lActRSU,OmegaS,SPot1,rKM,cMin,lOmegaW,
     !              iSolAU,lHanks,lDblVeg,cDataPath,rA,rB,rC,rD,lScreen,
     !              err)
        if(err.ne.0) goto (909,903) err
        lCouvreur=.false.
        if(lCouvreur) then
          call CouIn(StresFun,hx_min,Krs,Kcomp,NumNP,x,zSurf,BCtp,BCr,
     !               lCouvreur,lDayVar)
        end if
      end if
      if(TopInF.or.BotInF.or.AtmBC) then
        call SetBC(tMax,tAtm1,rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,
     !             TopInF,BotInF,cT,cBot,NS,tTop,tBot,Ampl,lTemp,lChem,
     !             KodTop,lVarBC,err,lMinStep,lMeteo,Prec,rSoil,lLAI,
     !             rExtinct,hSink,lCentrif,CosAlf,xConv,tConv,iModel,
     !             hNew(NumNP),iRootIn,xRoot,WLayer,lLinear,lActRSU,
     !             SPot,iSolAU,ParD(5,MatNum(NumNP)),nCycles,iCycle,
     !             tCycle,tInit,jIrrig,crWat,crSal,xLAI,lIsSolPres,
     !             lSurf,fET,lAPI,lAPI_ET)
        if(err.eq.1) goto 913
        if(tAtm1.lt.sngl(t)) goto 936
        if(lYield) rRootU=rRoot
        if(lYield.and.lChem) SPotU=SPot
        if(lMeteo) then
          call Meteo(1,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !               tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !               rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !               TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !               nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !               rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !               iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !               EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !               Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !               xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !               Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !               TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !               SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !               RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !               RH_A,EaMean,EaMeanN,rTop,err)
          if(err.ne.0) goto (932,913,933) err
        end if
        if(lDblVeg) call Root2(rRoot,rRoot2)
        if(lChem.and.SnowLayer.le.0)
     !    call SetChemBC(Prec,rSoil,NS,NSD,NumNP,cTop,cT,WLayer,
     !                   hNew(NumNP),KodTop,kTopCh,lIsotope,Conc)
        if(lTemp.and.lSnow) then
          dtSnow=sngl(tAtm1-tInit)
          call Snow(Prec,dtSnow,tTop,SnowMF,SnowSF,SnowLayer,rSoil,
     !              xConv,lMinStep,lChem,cTop,cT,NS,lIsotope)
        end if
        if(.not.lMeteo.and.iInterc.gt.0)
     !    call IntercepS(xConv,tConv,xLAI,rExtinct,rInterc,aInterc,
     !                   Prec,rRoot,ExcesInt)
        tAtm=dmin1(tAtm1,tAtm2)
        if(lDayVar.and..not.lVarBC) then
          rRootD=rRoot
          rSoilD=rSoil
          call DailyVar(tConv,t,rRoot,rRootD)
          call DailyVar(tConv,t,rSoil,rSoilD)
        end if
        if(lSinPrec) then
          PrecD=Prec
          tAtmOld=sngl(tInit)
          call SinPrec(t,tAtmOld,tAtm1,Prec,PrecD)
        end if
        if(KodTop.eq.-4) rTop=abs(rSoil)-abs(Prec)
      end if
      if(lSink) then
        if(.not.TopInF.and..not.BotInF.and..not.AtmBC) SPot=SPot1
        if(lCouvreur) then
          call SetSnkV(NumNP,Beta,zSurf,x,StresFun,Krs,Kcomp,hx_min,
     !                 BCr,BCtp,hNew,hTemp,Sink,rRoot,hRoot,vRoot)
        else
        call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,r2L,
     !              aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot,lMsSink,
     !              ThOld,ParD,dt,OmegaC,iModel,Con,lOmegaW,OmegaW,
     !              lHanks,t,lDblVeg,Sink2,Beta2,rRoot2,crWat,crSal)
        end if
      end if
      close(30)
      close(32)
      close(50)

      if(lPrint) then
        call Profil (NumNP,NMat,x,MatNum,xSurf,Beta,Ah,AK,ATh,thr,ths,
     !               ConSat,hSat,lScreen,err)
        if(err.eq.1) goto 923
        call NodOut (NumNP,NMat,hNew,ThOld,Con,x,xSurf,CosAlf,tInit,
     !               MatNum,Cap,AK,Sink,ConSat,NS,NSD,Conc,TempO,Sorb,
     !               Kappa,lBact,Sorb2,lVapor,lWTDep,ConLT,ConVT,ConVh,
     !               ThOld(NumNP),dt,iDualPor,ThNewIm,SinkIm,STrans,
     !               lDensity,lCentrif,Radius,lVaporOut,lSolub,lAPI,
     !               lChem,ChPar,lDualNEq,TDep,lLinear,lMobIm,lEquil,
     !               ths,lFumig,lDblVeg,Sink2,err)
        if(err.eq.1) goto 920
      end if
      call SubReg (NumNP,NMat,NLay,hNew,ThOld,ThOld,x,MatNum,LayNum,
     !             t-dt,dt,CosAlf,Con,lChem,Conc,ChPar,0,ths,wCumT,
     !             wCumA,cCumT,cCumA,wVolI,cVolI,WatIn,SolIn,lWat,lTemp,
     !             TempN,TPar,TDep,NS,NSD,Sorb,lLinear,lEquil,lMobIm,
     !             err,SubVol,Area,lPrint,lBact,Sorb2,lVapor,ThVOld,
     !             ThVNew,lWTDep,ConLT,ConVh,ConVT,iDualPor,ThNewIm,
     !             ThOldIm,lDensity,lCentrif,Radius,lSolub,lThAds,thr,
     !             lDualNEq,dummy1,dummy2,dummy3,dummy4,cPrevO,lFumig,
     !             lSource,CumSource)
      if(err.eq.1) goto 921
      if(lChem.or.lTemp.or.lFlux)
     !  call Veloc(NumNP,hOld,Con,x,CosAlf,vOld,ThOld,ThOld,Sink,dt,
     !             lVapor,lWTDep,ConLT,ConVT,ConVh,TempO,vVOld,ThVOld,
     !             ThVOld,lDensity,Conc,NSD,lCentrif,Radius,SinkDr)
      do 11 i=1,NumNP
        vNew(i) =vOld(i)
        ThNew(i)=ThOld(i)
        if(lVapor) ThVNew(i)=ThVOld(i)
        if(lVapor) vVNew(i) =vVOld(i)
11    continue
      if(lBact.and..not.lWat)
     !  call Exclusion(NumNP,NMat,NSD,ParD,ChPar,ThNew,vNew,ThOld,vOld)
      if(lOust) then ! soil erosion (DuPont)
        call OustIn(cDataPath,WatEr,SolEr,zEr)
        call OustOut(tInit,NumNP,lWat,lChem,NMat,NLay,ThNew,x,MatNum,
     !               LayNum,Conc,ChPar,TempN,NS,NSD,lLinear,TDep,WatEr,
     !               SolEr,zEr,lPrint)
        rSolM=0.
      end if
      if(lTillage) ! Tillage mixing
     !  call TillageIn(cDataPath,lTillage,rTillDepth,NPD,iTill)
      if(lSchwen) ! Time-variable scaling factors
     !  call ScalFacIn(cDataPath,lSchwen)

      write(*,*)'All input data have been read successfully.'
      write(*,*)'Beginning of numerical solution.'
      if(.not.lScreen) write(*,*)'Further screen output is disabled.'

*---- START gettim and getdat gnu alternative function ----
      call date_and_time(date,time,zone,dtv)
      hundred=100
      iYear=dtv(1)
      iMonth=dtv(2)
      iDay=dtv(3)
      iHours=dtv(5)
      iMins=dtv(6)
      iSecs=dtv(7)
      i100th=mod(dtv(8),hundred)
*      call getdat(iYear,iMonth,iDay)
*      call gettim(iHours,iMins,iSecs,i100th)
*---- END gettim and getdat gnu alternative function ----
      Rtime1=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)

*     Time stepping calculations ---------------------------------------
12    continue

*     Loop between water flow and heat and vapor transport
      iTemp=0
13    continue

      if(lEnBal) then
        EpsiT=0.5
        TempS=EpsiT*TempN(NumNP)+(1.-EpsiT)*TempO(NumNP)
        call Evapor(t,TempS,TMaxA,TMinA,Rad,hNew(NumNP),TempHeight,
     !              WindHeight,Wind,RHMean,HeatFl,rTop,Prec,tPeriod,
     !              rLat,Albedo,SunHours,ThNew(NumNP),xConv,tConv,
     !              iRadiation,Rns,Rnl,Rn,SensFlux,Evap,xLat,Const,
     !              iSunSh,r_H,lMetDaily,Rst,TempA,RH_A,ShWRadA,ShWRadB,
     !              rLWRadA,rLWRadB,iMetHour,rSoil)
        if(lIrrig.and.jIrrig.eq.2) rTop=-rIrrig+rSoil-Prec
      end if

      if(lCond) call Conductivity(NumNP,NMatD,Ak,MatNum,t,iCond,t0,cA,
     !                            cB,cC,ParD)

*     Solve water movement ---------------------------------------------
      if(lWat) then
        if(lSurfact) call Surfact(NumNP,NSD,Conc,ThNew,Ah,AK)
        ItCumO=ItCum
        call WatFlow(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,ParD,ParW,Con,Cap,ConSat,Ah,AK,ATh,hSat,hTemp,
     !               KodTop,KodBot,rTop,rBot,CosAlf,t,dt,x,Sink,P,R,S,
     !               FreeD,SeepF,qGWLF,Aqh,Bqh,GWL0L,hTop,hBot,hCritA,
     !               hCritS,WLayer,IterW,ItCum,TopInf,kTOld,kBOld,TolTh,
     !               TolH,MaxIt,dtMin,tOld,dtOpt,ConvgF,TheTab,ThNew,
     !               ThOld,thr,ths,lWTDep,TempN,Kappa,KappaO,AThS,ThRR,
     !               ConO,ConR,AKS,AhW,AThW,AKW,iHyst,iModel,qDrain,
     !               zBotDr,BaseGW,rSpacing,iPosDr,rKhTop,rKhBot,rKvTop,
     !               rKvBot,Entres,WetPer,zInTF,GeoFac,lTable,lVapor,
     !               xConv,tConv,ConLT,ConVT,ConVh,TauW,ThEq,ThVNew,
     !               ThVOld,nTabMod,iDualPor,ThNewIm,ThOldIm,SinkIm,
     !               vTop,TempO,iTemp,WTransf,lDensity,Conc,NSD,iEnhanc,
     !               lCentrif,Radius,hSeep,lSource,iNode,hSink,
     !               CumSource,lDblVeg,Sink2,iDrKod,SinkDr,WDrain)
        if(.not.ConvgF) then
          iNonConv=iNonConv+1
        else
          iNonConv=0
        end if
        if(lStopConv.and..not.ConvgF.and.iNonConv.ge.10) goto 935

        if(lMeteo.and.lMetDaily.and..not.lEnBal) ! Output for daily variated meteo information
     !    call DayMeteoOut(t,ETcomb,EvapP,TransP,Rns,Rnl,RadTerm,
     !                     AeroTerm,Precc,rInterc,ExcesInt,TempA,RH_A,
     !                     Rst,lPrint)
        if(lEnBal) then
          if(ItCum.gt.ItCumO+MaxIt.and..not.lMetDaily)
     !    call MeteoInt(3,t,tAtm2O,tAtm2,Rad,RadO,RadN,TMaxA,TMaxAO,
     !                  TMaxAN,TMinA,TMinAO,TMinAN,Wind,WindO,WindN,
     !                  RHMean,RHMeanO,RHMeanN,SunHours,SunHoursO,
     !                  SunHoursN,lEnBal)
          TempS=EpsiT*TempN(NumNP)+(1.-EpsiT)*TempO(NumNP)
          M=MatNum(NumNP)
          if(iCampbell.gt.0) then
            AA=(0.57+1.73*TPar(4,M)+0.93*TPar(5,M))/(1.-0.74*TPar(4,M)-
     !          0.49*TPar(5,M))-2.8*TPar(1,M)*(1.-TPar(1,M))
            BB=2.8*TPar(1,M)
            xc=amax1(0.005,TPar(6,M))
            CC=1.+2.6/sqrt(xc)
            DD=0.03+0.7*TPar(1,M)**2
            EE=4.
            rLamb=AA+BB*ThNew(NumNP)-(AA-DD)
     !            *exp(-(CC*ThNew(NumNP))**EE)
          else
            rLamb=amax1(0.,TPar(4,M)+TPar(5,M)*ThNew(NumNP)+
     !                  TPar(6,M)*sqrt(ThNew(NumNP)))
          end if
          dz=x(NumNP)-x(NumNP-1)
          call UpdateEnergy(t,vTop,rTop,HeatFl,TempS,Rns,Rnl,Rn,Evap,
     !                      xLat,SensFlux,xConv,tConv,Const,TLevel,
     !                      nPrStep,r_h,dz,rLamb,iTemp,lPrint,lMetDaily,
     !                      TempA,RH_A,Rst,ShortO,TPrint1,lPrintD,dt)
        end if
      else
        iterW=1
        ItCum=ItCum+1
      end if

*     To calculate the velocities --------------------------------------
      if(lWat.and.(lTemp.or.lChem.or.lFlux))
     !  call Veloc(NumNP,hNew,Con,x,CosAlf,vNew,ThNew,ThOld,Sink,dt,
     !             lVapor,lWTDep,ConLT,ConVT,ConVh,TempN,vVNew,ThVNew,
     !             ThVOld,lDensity,Conc,NSD,lCentrif,Radius,SinkDr)

*     Root zone calculations
      if(lRoot.or.iCrop.gt.1)
     !  call SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR,xRoot,
     !             lRoot,iRootIn,nGrowth,rGrowth,tRPeriod,fET)
      if(lSink) then
        if(lCouvreur) then
          call SetSnkV(NumNP,Beta,zSurf,x,StresFun,Krs,Kcomp,hx_min,
     !                 BCr,BCtp,hNew,hTemp,Sink,rRoot,hRoot,vRoot)
        else
        call SetSnk(NumNP,NMat,MatNum,x,hRoot,vRoot,Sink,rRoot,hNew,
     !              lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,
     !              r2L,aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot,
     !              lMsSink,ThNew,ParD,dt,OmegaC,iModel,Con,lOmegaW,
     !              OmegaW,lHanks,t,lDblVeg,Sink2,Beta2,rRoot2,crWat,
     !              crSal)
        end if
      end if

*     Calculation of heat transport ------------------------------------
      if(lTemp) then
        call Temper(NumNP,NMat,x,dt,t,MatNum,TempO,TempN,TPar,Ampl,P,R,
     !              S,Q,vOld,vNew,ThOld,ThNew,Retard,Disp,Sink,tPeriod,
     !              kTopT,tTop,kBotT,tBot,lVapor,ThVOld,ThVNew,vVOld,
     !              vVNew,g0,lEnBal,HeatFl,xConv,tConv,dtMaxT,iCampbell,
     !              iTemp,HeatFlB)
        if((lVapor.or.lEnBal).and.iTemp.lt.1) then
          iTemp=iTemp+1
          goto 13
        end if
      end if

*     Calculations of the solute transport -----------------------------
      if(lChem) then
        iKod=2
        call Solute(NumNP,NMat,NS,NSD,x,dt,t,tPulse,ChPar,MatNum,ThOld,
     !              ThNew,vOld,vNew,Disp,epsi,kTopCh,cTop,kBotCh,cBot,
     !              Conc,P,R,S,Q,g0,g1,Retard,cvTop,cvBot,cvCh0,cvCh1,
     !              lUpW,wc,Peclet,Courant,dtMaxC,TempO,TempN,cNew,
     !              cPrevO,cTemp,TDep,ths,cTolA,cTolR,IterC,MaxItC,
     !              hTemp,Sorb,SorbN,lLinear,lEquil,lArtD,PeCr,q0,q1,
     !              dSurf,cAtm,lTort,Sink,cRootMax,sSink,cvChR,lMobIm,
     !              cvChIm,TLevel,lBact,Sorb2,SorbN2,dtMin,dtOpt,lWat,
     !              lFiltr,iDualPor,ThOldIm,ThNewIm,SinkIm,STrans,
     !              lSolHyst,iTort,xConv,tConv,lVapor,rBot,err,lSolub,
     !              lThAds,thr,iMoistDep,NMatD,DMoist,WDep,iKod,Sorb3,
     !              Beta,lDualNEq,AtmBC,lSink,lActRSU,OmegaS,OmegaW,
     !              SPot,iSolAU,rKM,cMin,lDensity,lIsotope,lFumig,
     !              hNew(NumNP),lIsSolPres,rA,rB,rC,rD,qDrain,iDrKod,
     !              SinkDr,cvChDr)
        if(err.ne.0) goto 927
      end if

*     Output ------------------------------------------------------------
*     T-level information
      if(abs(t-tMax).le.0.5*dtMin.or.t.gt.tMax) lEnd=.true.
      jPrint=0
      if((abs(TPrint(PLevel)-t).lt.0.001*dt.or.
     !    (lPrintD.and.(abs(TPrint1-t).lt.0.001*dt)).or.
     !    (.not.ShortO.and.abs(float((TLevel+nPrStep-1)/nPrStep)-
     !    (TLevel+nPrStep-1)/float(nPrStep)).lt.0.0001))) jPrint=1

      call TLInf(NumNP,Con,x,CosAlf,t,dt,IterW,IterC,TLevel,rTop,rRoot,
     !           vRoot,hNew,hRoot,CumQ,ItCum,KodTop,KodBot,ConvgF,lWat,
     !           lChem,cRoot,NS,NSD,Conc,cvTop,cvBot,cvCh0,cvCh1,Peclet,
     !           Courant,wCumT,wCumA,cCumT,cCumA,CumCh,ThNew,ThOld,Sink,
     !           lScreen,err,cvChR,cvChIm,lPrint,lVapor,lWTDep,ConLT,
     !           ConVh,ConVT,TempN,rSoil,Prec,nPrStep,ThVOld,ThVNew,
     !           xConv,iDualPor,SinkIm,WTransf,lDensity,SnowLayer,
     !           lCentrif,Radius,ThNewIm,WLayer,hCritS,lEnd,lFluxOut,
     !           jPrint,dummy1,dummy2,NObs,Node,vNew,cNew,cTop,lDblVeg,
     !           rRoot2,PTLevel,lTemp,lPrHeat,HeatFl,HeatFlB,CumHeatFlT,
     !           CumHeatFlB,qDrain,iDrKod,WDrain,jIrrig,rIrrig,
     !           LayNum,NLay)  ! Adam S: calculation of volumes in subregions
      if(err.ne.0) goto (919,918,926) err
      if(NObs.gt.0.and.lPrint.and.jPrint.eq.1) then
        call ObsNod(t,NumNP,NObs,NS,NSD,Node,Conc,hNew,ThNew,TempN,
     !              lChem,ThNewIm,vNew,vVNew,lFlux,Sorb,lAPI,err)
        if(err.eq.1) goto 922
      end if
      if(lChem.and.
     !          lFumig.and.lAddFum.and.abs(t-tAddFum).le.0.5*dtMin) then
        call AddFumig(NumNP,NMat,ThNew,x,MatNum,Conc,ChPar,ths,TempN,
     !                NS,NSD,lLinear,lEquil,FumMass,zAddFT,zAddFB,
     !                TDep,cVolI,lMinStep,tAddFum,tMax)
      end if
      if(lPrintD.and.dabs(TPrint1-t).lt.0.001*dt)
     !                                         TPrint1=TPrint1+tPrintInt
      if(lYield)
     !  call Yield1(lChem,Yield,rRoot,vRoot,rSoil,cvChR,NS,SPot,CTp,CRp,
     !              dt,t,rRootU,SPotU,TLevel,0)

*     P-level information ----------------------------------------------
      if(abs(TPrint(PLevel)-t).lt.0.001*dt) then
*       Tillage mixing at print times
        if(lTillage.and.iTill(PLevel).eq.1) ! Tillage mixing
     !    call Tillage(NumNP,NMat,hNew,hOld,hTemp,ThNew,ThOld,x,MatNum,
     !                 ParD,iModel,lChem,NS,NSD,lLinear,lEquil,ChPar,
     !                 Conc,ths,TempN,TDep,Sorb,rTillDepth)

        if(lPrint) then
          call NodOut(NumNP,NMat,hNew,ThNew,Con,x,xSurf,CosAlf,
     !                TPrint(PLevel),MatNum,Cap,AK,Sink,ConSat,NS,NSD,
     !                Conc,TempN,Sorb,Kappa,lBact,Sorb2,lVapor,lWTDep,
     !                ConLT,ConVT,ConVh,ThOld(NumNP),dt,iDualPor,
     !                ThNewIm,SinkIm,STrans,lDensity,lCentrif,Radius,
     !                lVaporOut,lSolub,lAPI,lChem,ChPar,lDualNEq,TDep,
     !                lLinear,lMobIm,lEquil,ths,lFumig,lDblVeg,Sink2,
     !                err)
          if(err.eq.1) goto 920
          if(lFumig) then
            call FumigOut(cvTop,CumCh,tPrev,CumChPrev,TPrint(PLevel),NS,
     !                    NSD,NumNP,NMat,NObs,ChPar,TDep,Conc,TempN,
     !                    MatNum,Node,err)
            if(err.eq.1) goto 934
          end if
        end if
        call SubReg(NumNP,NMat,NLay,hNew,ThNew,ThOld,x,MatNum,LayNum,t,
     !              dt,CosAlf,Con,lChem,Conc,ChPar,PLevel,ths,wCumT,
     !              wCumA,cCumT,cCumA,wVolI,cVolI,WatIn,SolIn,lWat,
     !              lTemp,TempN,TPar,TDep,NS,NSD,Sorb,lLinear,lEquil,
     !              lMobIm,err,SubVol,Area,lPrint,lBact,Sorb2,lVapor,
     !              ThVOld,ThVNew,lWTDep,ConLT,ConVh,ConVT,iDualPor,
     !              ThNewIm,ThOldIm,lDensity,lCentrif,Radius,lSolub,
     !              lThAds,thr,lDualNEq,dummy1,dummy2,dummy3,dummy4,
     !              cPrevO,lFumig,lSource,CumSource)
        if(err.eq.1) goto 921
        if(lAPI.and.lChem) then
          cMin=0.0005
          cMaxR=0.
          cMax=0.
          do 14 i=1, NumNP
            if(Conc(1,i).gt.cMax) cMax=Conc(1,i)
            if(Conc(1,i).gt.cMaxR.and.Beta(i).gt.0.) cMaxR=Conc(1,i)
14        continue
          if(cMax.lt.cMin) tMax=TPrint(PLevel)
          if(cMaxR.lt.cRootIni/10.) lAPI_ET=.true.
        end if
        if(lSource) write(76,110) CumSource
        PLevel=PLevel+1
      end if

*     A-level information ----------------------------------------------
      if(dabs(t-tAtm).le.0.001*dt.and.(TopInF.or.BotInF.or.AtmBC)) then
        if(lPrint) then
          call ALInf(t,CumQ,hNew(NumNP),hRoot,hNew(1),ALevel,iCycle,err)
          if(err.ne.0) goto (925) err
        end if
        if(dabs(t-tAtm1).le.0.001*dt) then
          tAtmOld=tAtm1
          call SetBC(tMax,tAtm1,rTop,rR,rBot,hCritA,hBot,hTop,GWL0L,
     !               TopInF,BotInF,cT,cBot,NS,tTop,tBot,Ampl,lTemp,
     !               lChem,KodTop,lVarBC,err,lMinStep,lMeteo,Prec,rS,
     !               lLAI,rExtinct,hSink,lCentrif,CosAlf,xConv,tConv,
     !               iModel,hNew(NumNP),iRootIn,xRoot,WLayer,lLinear,
     !               lActRSU,SPot,iSolAU,ParD(5,MatNum(NumNP)),nCycles,
     !               iCycle,tCycle,tInit,jIrrig,crWat,crSal,xLAI,
     !               lIsSolPres,lSurf,fET,lAPI,lAPI_ET)
          if(err.eq.1) goto 913
          ii=0 ! switch whether to do adjustment (=1) or not (=0)
          if(ii.eq.1.and.KodBot.ge.0) then
            if(TLevel.eq.1) wAdded=0. ! amount of water addded due to this modification
c            call ModFlowBot(NumNP,hNew,hOld,hTemp,ThNew,ThOld,x,hBot,
c     !                      iModel,MatNum,ParD,NMat,wVolI,Conc,ChPar,
c     !                      NSD,wAdded)
            call ModFlowBot1(NumNP,hNew,hOld,hTemp,ThNew,ThOld,x,hBot,
     !                       NTab(1),NTabD,NMat,iModel,hTab,Contab,
     !                       hSat,ParD,MatNum,Con,wVolI,Conc,ChPar,NSD,
     !                       wAdded)
c            write(78,*) wAdded
          end if

*         Special option for Yong Li (dynamic ponding or runoff)
          if(lYong) then
            inquire(unit=35,opened=lOpen)
            if(lOpen) then
              read(35,*,err=321) rDummy,WLayer,hCritS
321           write(*,*) tAtm1,rDummy,WLayer,hCritS
            end if
          end if
          if(lOust) then            ! Dry Deposition
            if(cT(1).gt.0.) then !accumulate BC when no rain
              rSolM=rSolM+cT(1)
              cT(1)=0.
            end if
            if(rTop.lt.0..and.SnowLayer.le.0..and.tTop.gt.+2.) then ! apply BC when rain
              if(abs(tAtm1-tAtmOld).gt.0..and.abs(Prec).gt.0.)
     !          cT(1)=rSolM/Prec/(tAtm1-tAtmOld)
              rSolM=0.
            end if
            if(rBot.gt.0.) then     ! Wind Erosion
              call Erosion(tAtm1,NumNP,NMat,NSD,x,hNew,hOld,hTemp,ThNew,
     !                     Conc,TempN,MatNum,ChPar,hBot,wVolI,cVolI(1),
     !                     WatEr,SolEr,zEr)
            else if(dabs(t-tMax).gt.dt/2.) then
              call ErosionOut(tAtm1,cVolI(1),WatEr,SolEr,zEr)
            end if
          end if
          if(.not.lMeteo) then
            rRoot=rR
            rSoil=rS
          end if
          if(lSchwen) ! Time-variable scaling factors
     !      call ScalFac(NumNP,AK,Ah,ATh,MatNum)
        end if
        if(lMeteo.and.dabs(t-tAtm2).le.0.001*dt) then
          call Meteo(2,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !               tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !               rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !               TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !               nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !               rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !               iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !               EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !               Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !               xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !               Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !               TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !               SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !               RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !               RH_A,EaMean,EaMeanN,rTop,err)
          if(err.ne.0) goto (932,913,933) err
        end if
        if(lYield)
     !    call Yield1(lChem,Yield,rRoot,vRoot,rSoil,cvChR,NS,SPot,CTp,
     !                CRp,dt,t,rRootU,SPotU,TLevel,1)
        if(lDblVeg) call Root2(rRoot,rRoot2)
        tAtm=dmin1(tAtm1,tAtm2)
        if(lChem.and.SnowLayer.le.0)
     !    call SetChemBC(Prec,rSoil,NS,NSD,NumNP,cTop,cT,WLayer,
     !                   hNew(NumNP),KodTop,kTopCh,lIsotope,Conc)
        if(lTemp.and.lSnow) then
          dtSnow=sngl(tAtm1-tAtmOld)
          call Snow(Prec,dtSnow,tTop,SnowMF,SnowSF,SnowLayer,rSoil,
     !              xConv,lMinStep,lChem,cTop,cT,NS,lIsotope)
        end if
        if(.not.lMeteo.and.iInterc.gt.0)
     !    call IntercepS(xConv,tConv,xLAI,rExtinct,rInterc,aInterc,
     !                   Prec,rRoot,ExcesInt)
        if(lDayVar.and..not.lVarBC) then
          rRootD=rRoot
          rSoilD=rSoil
        end if
        if(lSinPrec) PrecD=Prec
        if(.not.lVarBC.and.KodTop.eq.-4) rTop=abs(rSoil)-abs(Prec)
        if(.not.lVarBC) rTop=abs(rSoil)-abs(Prec)
        ALevel=ALevel+1
      end if

      if(lOust.and.jPrint.eq.1)
     !  call OustOut(t,NumNP,lWat,lChem,NMat,NLay,ThNew,x,MatNum,LayNum,
     !               Conc,ChPar,TempN,NS,NSD,lLinear,TDep,WatEr,SolEr,
     !               zEr,lPrint)

      if(lIrrig) then
        call Irrig(t,dt,tAtmI,rTop,cIrrig,rRoot,hIrrig,jIrrig,tIrrig,
     !             rIrrig,iIrrigTr,nIrreg,hNew(Node(iIrrig)),lMinStep,
     !             dtMaxT,dtMin,rSoil,Prec)
        tAtm=dmin1(tAtm1,tAtmI)
      end if
      if(WLayer.and.hNew(NumNP).gt.0.) then ! mass balance in the surface layer
        hT=hNew(NumNP)
        do 15 jj=1,NS
          if((hT+dt*(Prec-rSoil)).gt.0.)
     !      cTop(jj)=(hT*cTop(jj)+dt*Prec*cT(jj))/(hT+dt*(Prec-rSoil))
15      continue
      end if

*     Time governing ---------------------------------------------------
      if(abs(t-tMax).le.0.5*dtMin.or.t.gt.tMax) then
        if(iCarda.eq.1) then
          call ICARDASet(t,tInit,tOld,dt,dtInit,nICarda,iICarda,TLevel,
     !                   ALevel,PLevel,tAtm,tAtm1,rTop,Prec,rSoil,rRoot,
     !                   cTop(1),cBot(1),CumQ,Conc(1,1),hRoot,cRoot(1),
     !                   crWat,crSal)
          if(iICarda.le.nICarda) goto 12
        end if
        call CloseOutput(RTime1,NS,TopInF,BotInF,lChem,lScreen,lMeteo,
     !                   lPrint,0,lAPI)
        write(*,*)'Calculations have finished successfully.'
        if(lEnter) then
          write(*,*) 'Calculation complete'
        end if
        stop
      else
        tOld=t
        dtOld=dt
        kTOld=KodTop
        kBOld=KodBot
        if(.not.lWat) IterW=1
        Iter=max0(IterW,IterC)
        if(lChem.and.abs(tPulse-t).lt.dtMin) lMinStep=.true.
        dtMaxA=min(dtMaxC,dtMaxT)
        if(lSinPrec.and.PrecD.gt.0)
     !    dtMaxA=min(dtMaxA,sngl(tAtm1-tAtmOld)/20.)
        if(lChem.and.lFumig.and.lAddFum) tAtm=dmin1(tAtm,tAddFum)
        call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,
     !              min(TPrint(PLevel),TPrint1),tAtm,t,tMax,dtMaxA,
     !              ItMin,ItMax,lMinStep,dtInit)
        t=t+dt
        if(lDayVar.and..not.lVarBC) then
          call DailyVar(tConv,t,rRoot,rRootD)
          call DailyVar(tConv,t,rSoil,rSoilD)
          rTop=abs(rSoil)-abs(Prec)
          if(jIrrig.eq.2) rTop=-rIrrig+abs(rSoil)-abs(Prec)
        end if
        if(lSinPrec) then
          call SinPrec(t,tAtmOld,tAtm1,Prec,PrecD)
          rTop=abs(rSoil)-abs(Prec)
          if(jIrrig.eq.2) rTop=-rIrrig+abs(rSoil)-abs(Prec)
        end if
        if(lMeteo) then
          call Meteo(3,lMetDaily,lDayVar,t,dt,tInit,tMax,tAtm2,
     !               tAtmN,tAtm2O,dtMax,rLat,rAlt,ShWRadA,ShWRadB,
     !               rLWRadA,rLWRadB,rLWRadA1,rLWRadB1,WindHeight,
     !               TempHeight,iCrop,iLAI,rRoot,xConv,tConv,rGrowth,
     !               nGrowth,iInterc,rInterc,aInterc,ExcesInt,lEnBal,
     !               rExtinct,lPrint,lHargr,iRadiation,iSunSh,iRelHum,
     !               iMetHour,CloudF_Ac,CloudF_Bc,Prec,Precc,rSoil,
     !               EvapP,TransP,Rns,Rnl,RadTerm,AeroTerm,Rst,ETcomb,
     !               Rad,RadN,RadO,Wind,WindN,WindO,Albedo,AlbedoN,
     !               xLAI,xLAIN,xRoot,xRootN,CropHeight,CropHeightN,
     !               Ampl,tTop,TMaxAN,TMinAN,TMax1,TMaxN,TMaxO,TMin1,
     !               TMinN,TMinO,TempA,TMaxA,TMaxAO,TMinA,TMinAO,
     !               SunHours,SunHoursN,SunHoursO,RHMean,RHMeanN,
     !               RHMeanO,RHMax,RHMaxN,RHMaxO,RHMin,RHMinN,RHMinO,
     !               RH_A,EaMean,EaMeanN,rTop,err)
          if(err.ne.0) goto (932,913,933) err
        end if
        TLevel=TLevel+1
        if(TLevel.gt.999999) TLevel=2
      end if

*     New updated values
      call Update(NumNP,lWat,lChem,lTemp,lVapor,iDualPor,lExtrap,dt,
     !            dtOld,hTemp,hNew,hOld,ThOld,ThNew,vOld,vNew,ThVOld,
     !            ThVNew,vVOld,vVNew,ThOldIm,ThNewIm,TempO,TempN,rTop,
     !            xConv,ConSMax,KodTop,KodBot,WLayer)
      goto 12

* --- End of time loop -------------------------------------------------

*     Error messages
901   ierr=1
      goto 1000
902   ierr=2
      goto 1000
903   ierr=3
      goto 1000
904   ierr=4
      goto 1000
905   ierr=5
      goto 1000
906   ierr=6
      goto 1000
907   ierr=7
      goto 1000
908   ierr=8
      goto 1000
909   ierr=9
      goto 1000
910   ierr=10
      goto 1000
911   ierr=11
      goto 1000
912   ierr=12
      goto 1000
913   ierr=13
      goto 1000
914   ierr=14
      goto 1000
915   ierr=15
      goto 1000
916   ierr=16
      goto 1000
c917   ierr=17  ! TODO: changes were here: no c
c      goto 1000
918   ierr=18
      goto 1000
919   ierr=19
      goto 1000
920   ierr=20
      goto 1000
921   ierr=21
      goto 1000
922   ierr=22
      goto 1000
923   ierr=23
      goto 1000
924   ierr=24
      goto 1000
925   ierr=25
      goto 1000
926   ierr=26
      goto 1000
927   ierr=27
      goto 1000
928   ierr=28
      goto 1000
929   ierr=29
      goto 1000
930   ierr=30
      goto 1000
931   ierr=31
      goto 1000
932   ierr=32
      goto 1000
933   ierr=33
      goto 1000
934   ierr=34
      goto 1000
935   ierr=35
      goto 1000
936   ierr=36
      goto 1000
937   ierr=37
      goto 1000

1000  call ErrorOut(ierr,cFileName,cDataPath,iLengthPath)
      call CloseOutput(RTime1,NS,TopInF,BotInF,lChem,lScreen,lMeteo,
     !                 lPrint,1,lAPI)
      if(lEnter) then
        write(*,*) 'Press Enter to continue'
        read(*,*)
      end if
      stop

101   format(a)
110   format('Cumulative flux at the nodal source/sink: ',f12.3)
      end

************************************************************************

      subroutine ErrorOut(ierr,cFileName,cDataPath,iLengthPath)

      character*260 cErr(60),cFileName,cDataPath,cFileNameErr

      cErr( 1)='Open file error in file :'
      cErr( 2)='File already exists or hard disk is full ! Open file err
     !or in output file : '
      cErr( 3)='Error when writing to an output file !'
      cErr( 4)='Error when reading from an input file Level_01.dir data
     !pathway !'
      cErr( 5)='Error when reading from an input file Selector.in Basic
     !Informations !'
      cErr( 6)='Error when reading from an input file Selector.in Water
     !Flow Informations !'
      cErr( 7)='Error when reading from an input file Selector.in Time I
     !nformations !'
      cErr( 8)='Error when reading from an input file Selector.in Root G
     !rowth Informations !'
      cErr( 9)='Error when reading from an input file Selector.in Sink I
     !nformations !'
      cErr(10)='Error when reading from an input file Selector.in Heat T
     !ransport Informations !'
      cErr(11)='Error when reading from an input file Selector.in Solute
     ! Transport Informations !'
      cErr(12)='Error when reading from an input file Profile.dat !'
      cErr(13)='Error when reading from an input file Atmosph.in !'
      cErr(14)='Dimension in NumNPD is exceeded !'
      cErr(15)='Dimension in NObsD is exceeded !'
      cErr(16)='Dimension in NMatD or NLay is exceeded !'
      cErr(17)='Error when writing into an output file I_Check.out !'
      cErr(18)='Error when writing into an output file Run_Inf.out !'
      cErr(19)='Error when writing into an output file T_Level.out !'
      cErr(20)='Error when writing into an output file Nod_Inf.out !'
      cErr(21)='Error when writing into an output file Balance.out !'
      cErr(22)='Error when writing into an output file Obs_Node.out !'
      cErr(23)='Error when writing into an output file Profile.out !'
      cErr(24)='Initial water content condition is lower than Qr !'
      cErr(25)='Error when writing into an output file A_Level.out !'
      cErr(26)='Error when writing into an output file Solute.out !'
      cErr(27)='Does not converge in the solute transport module ! Check
     ! iteration criteria in the Solute Transport - General Information
     !window!'
      cErr(28)='Number of Print-Times is exceeded !'
      cErr(29)='Dimension in NSD is exceeded !'
      cErr(30)='The path to the project is too long !!!'
      cErr(31)='Bulk density can not be equal to zero !'
      cErr(32)='Error when reading from an input file Meteo.in !'
      cErr(33)='Crop Height must be smaller than the height of wind and
     !temperature measurements!'
      cErr(34)='Error when writing into an output file Fumig.out !'
      cErr(35)='Numerical solution has been stopped since it has not con
     !verged in 10 consecutive time steps !'
      cErr(36)='The first time-variable BC record is at time smaller tha
     !n tInit+dtInit ! You may want to lower the initial time step.'
      cErr(37)='Dimension in NTabD in Mater.in exceeded !'

      cFileNameErr = cDataPath(1:iLengthPath)//'/ERROR.MSG'
      open(99,file=cFileNameErr,status='unknown',err=901)
      if(ierr.le.2) then
c        if(lScreen) write(*,*) cErr(ierr),cFileName
        write( *,*) cErr(ierr),cFileName
        write(99,*) cErr(ierr),cFileName
      else
c        if(lScreen) write(*,*) cErr(ierr)
        write( *,*) cErr(ierr)
        write(99,*) cErr(ierr)
      end if
      close(99)
      return

901   write(*,*) 'Folder with input data of the specified project does n
     !ot exist or pathway is too long or corrupted'
      write(*,*) cFileName
      return
      end

************************************************************************

*     iKod=0: Closing after finishing calculations
*     iKod=1: Closing due to error

      subroutine CloseOutput(RTime1,NS,TopInF,BotInF,lChem,lScreen,
     !                       lMeteo,lPrint,iKod,lAPI)

*      use MSFLIB
      logical TopInF,BotInF,lChem,lScreen,lMeteo,lPrint,lOpen,lAPI
      integer*2 iYear,iMonth,iDay,iHours,iMins,iSecs,i100th
      double precision RTime1,RTime2,RTime

*---- START gettim and getdat gnu alternative function ----
      integer dtv(8)
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer hundred

      call date_and_time(date,time,zone,dtv)
      hundred=100
      iYear=dtv(1)
      iMonth=dtv(2)
      iDay=dtv(3)
      iHours=dtv(5)
      iMins=dtv(6)
      iSecs=dtv(7)
      i100th=mod(dtv(8),hundred)
*      call getdat(iYear,iMonth,iDay)
*      call gettim(iHours,iMins,iSecs,i100th)
*---- END gettim and getdat gnu alternative function ----

      Rtime2=RTime(iMonth,iDay,iHours,iMins,iSecs,i100th)
      if(iKod.eq.0.and.lScreen) write(*,*)'Run time [sec]',Rtime2-RTime1
      if(lPrint) then
        write(70,'(''end'')')
        if(iKod.eq.0) then
          write(76,*)
          write(76,*) 'Calculation time [sec]',Rtime2-RTime1
        end if
        write(71,'(''end'')')
        if(TopInF.or.BotInF) then
          inquire(unit=72,opened=lOpen)
          if(lOpen) write(72,'(''end'')')
          inquire(unit=43,opened=lOpen)
          if(lOpen.and.lMeteo) write(43,'(''end'')')
        end if
        write(77,'(''end'')')
        if(lChem) then
          do 13 jj=1,NS
            inquire(unit=80+jj,opened=lOpen)
            if(lOpen) write(80+jj,'(''end'')')
13        continue
        end if
        if(lAPI.and.lChem) write(41,'(''end'')')
      end if
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

*-----------------------------------------------------------------------
*     Signal handler routines 
*-----------------------------------------------------------------------
*      integer(4) function h_sig (signum)
*        use MSFLIB
*        !MS$ATTRIBUTES C :: h_sig
*        integer(2) signum
*        select case (signum)
*          case (SIG$ABORT)
*            write(*,*) 'Abnormal termination !'
*          case (SIG$INT)
*            write(*,*) 'Program terminated !'
*          case (SIG$FPE)
*            write(*,*) 'Floating point error !'
*          case (SIG$BREAK)
*            write(*,*) 'Program terminated !'
*          case default
*            write(*,*) 'Unknown signal handler !'
*        end select
*        h_sig = 1
*        call CloseFiles
*        stop
*      end
*
*      function hand_fpe (signum, excnum)
*        !MS$ATTRIBUTES C :: hand_fpe
*        use MSFLIB
*        integer(2)  signum, excnum
*        select case(excnum)
*          case(FPE$INVALID )
*            write(*,*) 'Floating Point Error - Invalid number !'
*          case( FPE$DENORMAL )
*            write(*,*) 'Floating Point Error - Denormalized number !'
*          case( FPE$ZERODIVIDE )
*            write(*,*) 'Floating Point Error - Zero divide !'
*          case( FPE$OVERFLOW )
*            write(*,*) 'Floating Point Error - Overflow !'
*          case( FPE$UNDERFLOW )
*            write(*,*) 'Floating Point Error - Underflow !'
*          case( FPE$INEXACT )
*            write(*,*) 'Floating Point Error - Inexact precision !'
*          case default
*            write(*,*) 'Floating Point Error - Non-IEEE type !'
*        end select
*        hand_fpe = 1
*        call CloseFiles
*      end
*
************************************************************************

      subroutine CloseFiles
      logical lOpen

      inquire(unit=43,opened=lOpen)
      if(lOpen) then
        write(43,'(''end'')')
        close(43)
      end if
      inquire(unit=44,opened=lOpen)
      if(lOpen) then
        write(44,'(''end'')')
        close(44)
      end if
      inquire(unit=70,opened=lOpen)
      if(lOpen) then
        write(70,'(''end'')')
        close(70)
      end if
      inquire(unit=71,opened=lOpen)
       if(lOpen) then
        write(71,'(''end'')')
        close(71)
      end if
      inquire(unit=72,opened=lOpen)
      if(lOpen) then
        write(72,'(''end'')')
        close(72)
      end if
      inquire(unit=77,opened=lOpen)
      if(lOpen) then
        write(77,'(''end'')')
        close(77)
      end if
      inquire(unit=75,opened=lOpen)
      if(lOpen) close(75)
      inquire(unit=76,opened=lOpen)
      if(lOpen) close(76)
      inquire(unit=78,opened=lOpen)
      if(lOpen) close(78)
      inquire(unit=31,opened=lOpen)
      if(lOpen) close(31)
      inquire(unit=81,opened=lOpen)
      if(lOpen) then
        write(81,'(''end'')')
        close(81)
      end if
      inquire(unit=82,opened=lOpen)
      if(lOpen) then
        write(82,'(''end'')')
        close(82)
      end if
      inquire(unit=83,opened=lOpen)
      if(lOpen) then
        write(83,'(''end'')')
        close(83)
      end if
      inquire(unit=41,opened=lOpen)
      if(lOpen) close(41)
      inquire(unit=42,opened=lOpen)
      if(lOpen) close(42)

      return
      end

************************************************************************
