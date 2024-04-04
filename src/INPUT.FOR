* Source file INPUT.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine BasInf(CosAlf,MaxIt,TolTh,TolH,TopInF,BotInF,ShortO,
     !                  lWat,lChem,SinkF,WLayer,qGWLF,FreeD,SeepF,AtmBC,
     !                  KodTop,KodBot,rTop,rRoot,rBot,hCritS,hCritA,
     !                  GWL0L,Aqh,Bqh,kTOld,kBOld,NUnitD,iUnit,NMat,
     !                  NMatD,NLay,lRoot,lTemp,lWDep,lEquil,lScreen,
     !                  qDrain,zBotDr,BaseGW,rSpacing,iPosDr,rKhTop,
     !                  rKhBot,rKvTop,rKvBot,Entres,WetPer,zInTF,GeoFac,
     !                  lInitW,lVarBC,xConv,tConv,lMeteo,lVapor,iVer,
     !                  lPrint,lCentrif,lSnow,hSeep,lFlux,lActRSU,
     !                  lIrrig,iDrKod,ierr)

      character*72 Hed
      character*5 LUnit,TUnit,MUnit
      logical TopInF,BotInF,ShortO,lWat,lChem,lTemp,SinkF,WLayer,qGWLF,
     !        FreeD,SeepF,AtmBC,lRoot,lWDep,lEquil,lScreen,qDrain,
     !        lInitW,lVarBC,lPrint,lMeteo,lVapor,lCentrif,lSnow,lFlux,
     !        lActRSU,lIrrig,lDummy
      dimension iUnit(NUnitD)

*---- START gettim and getdat gnu alternative function PARAMETERS ----
      integer dtv(8)
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer hundred
*---- START gettim and getdat gnu alternative function PARAMETERS----

      iVer = iGetFileVersion(30,1)
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,'(a)',err=901) Hed
      read(30,*,err=901)
      read(30,'(a)',err=901) LUnit
      read(30,'(a)',err=901) TUnit
      read(30,'(a)',err=901) MUnit
      call Conversion(LUnit,TUnit,xConv,tConv)
      read(30,*,err=901)
      read(30,*,err=901) lWat,lChem,lTemp,SinkF,lRoot,ShortO,lWDep,
     !                   lScreen,AtmBC,lEquil
      if(iVer.eq.3) then
        read(30,*,err=901)
        read(30,*,err=901) lSnow,lDummy,lDummy,lDummy
      else if(iVer.ge.4) then
        read(30,*,err=901)
        read(30,*,err=901) lSnow,lDummy,lMeteo,lVapor,lActRSU,lFlux,
     !                     lIrrig
      end if
      if(lSnow.and..not.lTemp) lSnow=.false.
      read(30,*,err=901)
      read(30,*,err=901) NMat,NLay,CosAlf
      if(NMat.gt.NMatD.or.NLay.gt.10) then
        ierr=4
        return
      end if
      read(30,*,err=902)
      read(30,*,err=902)
      read(30,*,err=902) MaxIt,TolTh,TolH
      read(30,*,err=902)
      read(30,*,err=902) TopInF,WLayer,KodTop,lInitW
      if(KodTop.eq.0) lVarBC=.true.
      read(30,*,err=902)
      if(iVer.le.3) then
        ii=1
        read(30,*,err=904) BotInF,qGWLF,FreeD,SeepF,KodBot,qDrain
        ii=0
904     if(ii.eq.1) qDrain=.false.
        hSeep=0.
      else
        read(30,*,err=902) BotInF,qGWLF,FreeD,SeepF,KodBot,qDrain,hSeep
      end if
      if((.not.TopInF.and.KodTop.eq.-1).or.
     !   (.not.BotInF.and.KodBot.eq.-1.and.
     !   .not.qGWLF.and..not.FreeD.and..not.SeepF.and..not.qDrain)) then
        read(30,*,err=902)
        read(30,*,err=902) rTop,rBot,rRoot  !,hCritS,hCritA
        if(qGWLF.or.FreeD.or.SeepF.or.qDrain) rBot=0.
      else
        rTop=0.
        rBot=0.
        rRoot=0.
      end if
      if(qGWLF) then
        read(30,*,err=902)
        read(30,*,err=902) GWL0L,Aqh,Bqh
      end if
      if(qDrain) then
        read(30,*,err=902)
        read(30,*,err=902) iPosDr
        iDrKod=0
        if(iPosDr.gt.10) then
          iDrKod=1
          iPosDr=iPosDr-10
        end if
        read(30,*,err=902)
        read(30,*,err=902) zBotDr,rSpacing,Entres
        zBotDr=-abs(zBotDr)
        read(30,*,err=902)
        if ((iPosDr.eq.1).or.(iPosDr.eq.6)) then
          read(30,*,err=902) rKhTop
        else if(iPosDr.eq.2) then
          read(30,*,err=902) BaseGW,rKhTop,WetPer
        else if(iPosDr.eq.3) then
          read(30,*,err=902) BaseGW,rKhTop,rKhBot,WetPer
        else if(iPosDr.eq.4) then
          read(30,*,err=902) BaseGW,rKvTop,rKvBot,rKhBot,WetPer,zInTF
        else if(iPosDr.eq.5) then
          read(30,*,err=902) BaseGW,rKhTop,rKvTop,rKhBot,WetPer,zInTF,
     !                       GeoFac
        end if
        BaseGW=-abs(BaseGW)
        zInTF=-abs(zInTF)
      end if

*     Input modifications
      rRoot=abs(rRoot)
      hCritA=1.e+10
      hCritA=-abs(hCritA)
      if(TopInF) KodTop=isign(3,KodTop)
      if(BotInF) KodBot=isign(3,KodBot)
      if(AtmBC.and.KodTop.lt.0) then
        hCritS=0
        KodTop=-4
      end if
      if(WLayer) KodTop=-iabs(KodTop)
      if(qGWLF)  KodBot=-7
      if(FreeD)  KodBot=-5
      if(SeepF)  KodBot=-2
      kTOld=KodTop
      kBOld=KodBot

      if(lScreen) then
        write(*,*)'----------------------------------------------------'
        write(*,*)'|                                                  |'
        write(*,*)'|                    HYDRUS                        |'
        write(*,*)'|                                                  |'
        write(*,*)'|   Code for simulating one-dimensional variably   |'
        write(*,*)'|    saturated water flow, heat transport, and     |'
        write(*,*)'|   transport of solutes involved in sequential    |'
        write(*,*)'|         first-order decay reactions              |'
        write(*,*)'|                                                  |'
        write(*,*)'|         Based on updated version 4.14            |'
        write(*,*)'|     enhanced to run on linux-based systems       |'
        write(*,*)'|                                                  |'
        write(*,*)'|           Last modified: Dec, 2023               |'
        write(*,*)'|                                                  |'
        write(*,*)'----------------------------------------------------'
        write(*,*)
        write(*,*) Hed
        write(*,*)
c        write(*,*) 'Press Enter to continue'
c        read(*,*)
      end if
      ii=1
      if(lPrint) ii=NUnitD
      do 11 i=1,ii
        write(iUnit(i),*,err=903)'******* Program HYDRUS'
        write(iUnit(i),*,err=903)'******* ',Hed
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
*        call getdat(ii,imonth,iday)
*        call gettim(ihours,mins,isecs,ii)
*---- END gettim and getdat gnu alternative function ----
        write(iUnit(i),100,err=903) iday,imonth,iYear,ihours,mins,isecs
        write(iUnit(i),*,err=903)'Units: L = ',LUnit,', T = ',TUnit,
     !                                ', M = ',MUnit
11    continue

      if(CosAlf.le.1.) lCentrif=.false.
      if(lCentrif) then
        g=9.80665*xConv/tConv/tConv 
        CosAlf=CosAlf*CosAlf/g
      end if

      write(50,*,err=903)
      write(50,*,err=903) 'CosAlf,MaxIt,TolTh,  TolH'
      write(50,110,err=903) CosAlf,MaxIt,TolTh,TolH
      write(50,*,err=903)
      write(50,*,err=903) 'TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,FreeD,
     !SeepF,lWat,lChem,lTemp,lRoot,lWDep'
      write(50,120,err=903) TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,
     !                      FreeD,SeepF,lWat,lChem,lTemp,lRoot,lWDep
      return

*     Error when reading from an input file
901   ierr=1
      return
902   ierr=2
      return
*     Error when writing into an output file
903   ierr=3
      return

100   format(' Date: ',i3,'.',i2,'.',i4,'    Time: ',i3,':',i2,':'
     !,i2)
110   format(f6.3,i5,f8.3,f8.5)
120   format(13l6)
      end

************************************************************************

      subroutine Conversion(LUnit,TUnit,xConv,tConv)

*     conversions from m and s to Hydrus units

      character LUnit*5,TUnit*5
      xConv=1.
      tConv=1.
      if     (LUnit.eq."cm  ") then
        xConv=100.
      else if(LUnit.eq."mm  ") then
        xConv=1000.
      end if
      if     (TUnit.eq."min ") then
        tConv=1./60.
      else if(TUnit.eq."hours") then
        tConv=1./(60.*60.)
      else if(TUnit.eq."days") then
        tConv=1./(60.*60.*24.)
      else if(TUnit.eq."years") then
        tConv=1./(60.*60.*24.*365.)
      end if
      return
      end

************************************************************************

      subroutine NodInf(NumNPD,NumNP,hTop,hBot,x,hNew,hOld,MatNum,hTemp,
     !                  LayNum,Beta,Ah,AK,ATh,Conc,Sorb,TempN,TempO,NSD,
     !                  NS,xSurf,lChem,lTemp,lEquil,lScreen,lBact,Sorb2,
     !                  lSolub,ierr,lAPI,cRootIni,lDualNEq,lDblVeg,
     !                  Beta2)

      dimension x(NumNPD),hNew(NumNPD),hOld(NumNPD),MatNum(NumNPD),
     !          hTemp(NumNPD),LayNum(NumNPD),Beta(NumNPD),Ah(NumNPD),
     !          AK(NumNPD),ATh(NumNPD),Conc(NSD,NumNPD),TempN(NumNPD),
     !          TempO(NumNPD),Sorb(NSD,NumNPD),SConc(5),
     !          SSorb(5),C(5),S(5),Sorb2(NSD,NumNPD),Beta2(NumNPD)
      logical lChem,lTemp,lEquil,lScreen,lBact,lSolub,lAPI,lDualNEq,
     !        lDblVeg

      if(lScreen) write(*,*) 'reading nodal information'
      iVer = iGetFileVersion(32,1)
      read(32,*,err=901) n
      do 11 i=1,n
        read(32,*,err=901)
11    continue
      read(32,*,err=901) NumNP,ii,NS
      if(NumNP.gt.NumNPD) then
        ierr=3
        return
      end if
      if(NS.gt.NSD) then
        ierr=4
        return
      end if

*     Read nodal point information
      j=NumNP+1
12    continue
      j=j-1
      if(.not.lChem.and..not.lTemp) then
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx
        Te=20.
      else if(.not.lChem) then
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx,Te
      else if(lEquil) then
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx,Te,(C(ii),ii=1,NS)
      else
        read(32,*,err=901) n,x1,h,M,L,B,Ax,Bx,Dx,Te,(C(ii),ii=1,NS),
     !                    (S(ii),ii=1,NS)
      end if
      n=NumNP-n+1
      x(n)=x1
      hOld(n)=h
      MatNum(n)=M
      LayNum(n)=L
      Beta(n)=B
      Ah(n)=Ax
      AK(n)=Bx
      if(.not.lDblVeg.and.abs(Dx).lt.1.e-4) Dx=1.
      ATh(n)=Dx
      TempO(n)=Te
      do 1 ii=1,NS
        if(lChem) then
          Conc(ii,n)=C(ii)
          if(.not.lEquil)                 Sorb (ii,n)=S(ii)
          if(lBact.or.lSolub.or.lDualNEq) Sorb2(ii,n)=0.
          if(lSolub) then
            Sorb2(ii,n)=Sorb(ii,n)
            Sorb(ii,n)=0.
          end if
        end if
1     continue

      if(j-n) 13,18,14
13    write(*,*)'ERROR in NodInf at node =', n
      stop
14    continue
      dx=x(nOld)-x(n)
      ShOld=(hOld(nOld)-hOld(n))/dx
      SBeta=(Beta(nOld)-Beta(n))/dx
      SAh=(Ah(nOld)-Ah(n))/dx
      SAK=(AK(nOld)-AK(n))/dx
      SATh=(ATh(nOld)-ATh(n))/dx
      STemp=(TempO(nOld)-TempO(n))/dx
      if(lChem) then
        do 15 ii=1,NS
          SConc(ii)=(Conc(ii,nOld)-Conc(ii,n))/dx
          SSorb(ii)=(Sorb(ii,nOld)-Sorb(ii,n))/dx
15      continue
      end if
      do 17 i=nOld-1,n+1,-1
        dx=x(nOld)-x(i)
        hOld(i)=hOld(nOld)-ShOld*dx
        Beta(i)=Beta(nOld)-SBeta*dx
        Ah(i)=Ah(nOld)-SAh*dx
        AK(i)=AK(nOld)-SAK*dx
        ATh(i)=ATh(nOld)-SATh*dx
        TempO(i)=TempO(nOld)-STemp*dx
        if(lChem) then
          do 16 ii=1,NS
            Conc(ii,i)=Conc(ii,nOld)-SConc(ii)*dx
            Sorb(ii,i)=Sorb(ii,nOld)-SSorb(ii)*dx
            if(lBact.or.lSolub.or.lDualNEq) 
     !        Sorb2(ii,i)=Sorb2(ii,nOld)-SSorb(ii)*dx
16        continue
        end if
        MatNum(i)=MatNum(i+1)
        LayNum(i)=LayNum(i+1)
17    continue
      j=n
18    continue
      nOld=n
      if(j.gt.1) goto 12

      SBeta1=0.
      if(Beta(NumNP).gt.0.) SBeta1=Beta(NumNP)*(x(NumNP)-x(NumNP-1))/2.
      SBeta2=0.
      if(lDblVeg.and.Ath(NumNP).gt.0.)
     !                      SBeta2=Ath(NumNP)*(x(NumNP)-x(NumNP-1))/2. 
      do 19 i=2,NumNP-1
        if(Beta(i).gt.0.) SBeta1=SBeta1+Beta(i)*(x(i+1)-x(i-1))/2.
        if(lDblVeg.and.Ath(i).gt.0.) 
     !                    SBeta2=SBeta2+Ath(i)*(x(i+1)-x(i-1))/2.
19    continue
      do 20 i=2,NumNP
        if(SBeta1.gt.0.) then
          if(Beta(i).gt.0.) Beta(i)=Beta(i)/SBeta1
        else
          Beta(i)=0.
        end if
        if(lDblVeg.and.SBeta2.gt.0.) then
          Beta2(i)=Ath(i)/SBeta2
        else
          Beta2(i)=0.
        end if
20    continue
      xSurf=x(NumNP)
      if(lDblVeg) then
        do 23 i=1,NumNP
          Ath(i)=1.
23      continue
      end if

*     Print nodal information
      write(50,110,err=902)
      cRootIni=0.
      do 21 n=NumNP,1,-1
        if(.not.lChem.and..not.lTemp) then
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n)
        else if(.not.lChem) then
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n),
     !                          TempO(n)
        else if(lEquil) then
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n),
     !                          TempO(n),(Conc(ii,n),ii=1,NS)
        else
          write(50,120,err=902) NumNP-n+1,x(n),hOld(n),MatNum(n),
     !                          LayNum(n),Beta(n),Ah(n),AK(n),ATh(n),
     !                          TempO(n),(Conc(ii,n),ii=1,NS),
     !                          (Sorb(ii,n),ii=1,NS)
        end if
        hNew(n) =hOld(n)
        hTemp(n)=hOld(n)
        TempN(n)=TempO(n)
        if(lAPI.and.lChem.and.Beta(n).gt.0..and.Conc(1,n).gt.cRootIni) 
     !     cRootIni=Conc(1,n)  
21    continue
      write(50,'(''end'')',err=902)
      hBot=hNew(1)
      hTop=hNew(NumNP)
      write(50,130,err=902) NS
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format (/'Nodal point information'//
     !'node      x         hOld    MatN LayN  Beta      Ah       AK ',
     !'     ATh     Temp    Conc(1...NS)         Sorb(1...NS)'/)
120   format (i4,2f11.3,2i5,f8.3,3f9.3,f8.2,10e12.4,10e12.4)
130   format (/' Number of species in the chain : ',i3)
      end

************************************************************************

      subroutine ObsNodIn(NumNP,NObsD,NObs,lPrint,lChem,Node,NS,lFlux,
     !                    lAPI,ierr)

      dimension Node(NObsD)
      logical lPrint,lChem,lFlux,lAPI
      character*30 Text1,Text2,Text3

      read(32,*,err=901) NObs
      if(NObs.gt.NObsD) then
        ierr=3
        return
      end if
      if(NObs.gt.0) then
        read(32,*,err=901) (Node(i),i=1,NObs)
        do 11 i=1,NObs
          Node(i)=NumNP-Node(i)+1
11      continue
        if(lPrint) then
          NObsA=min(10,NObs)
          Text3='Node('
          if(lAPI.and.lChem) then
            Text1='    Flux         Conc    '
            write(41,301,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            write(41,302,err=902) (Text1,i=1,NObsA)
          end if
          Text1='      h        theta    Temp   '
          if(lFlux) Text1='      h        theta    Flux   '
          Text2='   Conc     '
          if(.not.lChem) then
            write(77,140,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            write(77,200,err=902) (Text1,i=1,NObsA)
          else 
            if(NS.eq.1)
     !        write(77,150,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.2)
     !        write(77,160,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.3)
     !        write(77,170,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.4)
     !        write(77,180,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.5)
     !        write(77,190,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.6)
     !        write(77,260,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.7)
     !        write(77,261,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.8)
     !        write(77,262,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.9)
     !        write(77,263,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.10)
     !        write(77,264,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.11)
     !        write(77,265,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.12)
     !        write(77,266,err=902) (Text3,NumNP-Node(j)+1,j=1,NObsA)
            if(NS.eq.1)
     !        write(77,210,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)
            if(NS.eq.2)
     !        write(77,220,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)
            if(NS.eq.3)
     !        write(77,230,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)
            if(NS.eq.4)
     !        write(77,240,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)
            if(NS.eq.5)
     !        write(77,250,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.6)
     !        write(77,270,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.7)
     !        write(77,271,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.8)
     !        write(77,272,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.9)
     !        write(77,273,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.10)
     !        write(77,274,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.11)
     !        write(77,275,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
            if(NS.eq.12)
     !        write(77,276,err=902) (Text1,(Text2,j=1,NS),i=1,NObsA)  
          end if
        end if
      end if
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

140   format (///16x,10(15x,a5,i3,')', 7x))
150   format (///16x,10(15x,a5,i3,')',18x))
160   format (///16x,10(15x,a5,i3,')',29x))
170   format (///16x,10(15x,a5,i3,')',40x))
180   format (///16x,10(15x,a5,i3,')',51x))
190   format (///16x,10(15x,a5,i3,')',62x))
260   format (///16x,10(15x,a5,i3,')',73x))
261   format (///14x,10(15x,a5,i3,')',84x))
262   format (///14x,10(15x,a5,i3,')',95x))
263   format (///14x,10(15x,a5,i3,')',106x))
264   format (///14x,10(15x,a5,i3,')',117x))
265   format (///14x,10(15x,a5,i3,')',128x))
266   format (///14x,10(15x,a5,i3,')',139x))
200   format (/'         time     ',10(a31,     2x))
210   format (/'         time     ',10(a31, a12,2x))
220   format (/'         time     ',10(a31,2a12,2x))
230   format (/'         time     ',10(a31,3a12,2x))
240   format (/'         time     ',10(a31,4a12,2x))
250   format (/'         time     ',10(a31,5a12,2x))
270   format (/'         time     ',10(a31,6a12,2x))
271   format (/'       time     ',10(a31,7a12,2x))
272   format (/'       time     ',10(a31,8a12,2x))
273   format (/'       time     ',10(a31,9a12,2x))
274   format (/'       time     ',10(a31,10a12,2x))
275   format (/'       time     ',10(a31,11a12,2x))
276   format (/'       time     ',10(a31,12a12,2x))
301   format (///16x,10(9x,a5,i3,')', 8x))
302   format (/'         time     ',10(a23,3x))
      end

************************************************************************

      subroutine InitW(NumNP,NMat,Matnum,Kappa,hNew,hOld,hTemp,ParD,
     !                 ParW,iModel,hTop,hBot,iDualPor,ThNewIm,ierr)

      dimension MatNum(NumNP),Kappa(NumNP),ParD(11,NMat),ParW(11,NMat),
     !          hNew(NumNP),hOld(NumNP),hTemp(NumNP),ThNewIm(NumNP)

      do 11 i=1,NumNP
        M=MatNum(i)
        ThTotal=hNew(i)
        if(iDualPor.gt.0) then
          hNew(i)=hNew(i)*ParD(2,M)/(ParD(2,M)+ParD(8,M))
          ThMobile=hNew(i)
          ThNewIm(i)=ThTotal-ThMobile
        end if
        if(Kappa(i).eq.-1) then
          Qe=min((hNew(i)-ParD(1,M))/(ParD(2,M)-ParD(1,M)),1.)
	    if(Qe.lt.0.) goto 901
          hNew(i)=FH(iModel,Qe,ParD(1,M))
        else
          Qe=min((hNew(i)-ParW(1,M))/(ParW(2,M)-ParW(1,M)),1.)
	    if(Qe.lt.0.) goto 901
          hNew(i)=FH(iModel,Qe,ParW(1,M))
        end if
        hOld(i)=hNew(i)
        hTemp(i)=hNew(i)
11    continue
      hBot=hNew(1)
      hTop=hNew(NumNP)
      return

901   ierr=1
      return
      end

************************************************************************

      subroutine InitDualPor(NumNP,NMat,MatNum,Par,theta,iDualPor,
     !                       ThNewIm,ThOldIm,SinkIm,hNew,STrans,lInitW)

      logical lInitW
      dimension MatNum(NumNP),Par(11,NMat),theta(NumNP),ThNewIm(NumNP),
     !          ThOldIm(NumNP),SinkIm(NumNP),hNew(NumNP),STrans(NumNP)

      do 11 i=1,NumNP
        M=MatNum(i)
        if(iDualPor.eq.0) ThNewIm(i)=0.
        if(.not.lInitW) then
          if(iDualPor.eq.1) then
            Se=(Theta(i)-Par(1,M))/(Par(2,M)-Par(1,M))
	      ThNewIm(i)=Par(7,M)+Se*(Par(8,M)-Par(7,M))
          else if(iDualPor.eq.2) then
            ThNewIm(i)=FQ(0,hNew(i),Par(7,M))
          end if
        end if
        ThOldIm(i)=ThNewIm(i)
        SinkIm(i)=0.
        STrans(i)=0.
11    continue

      return
      end

************************************************************************

      subroutine MatIn(NMat,ParD,ParW,hTab1,hTabN,lScreen,ierr,NumNP,Ah,
     !                 iHyst,AhW,AThW,AKW,MatNum,hNew,Kappa,AThS,ThRR,
     !                 ConR,AKS,KappaO,iModel,xConv,lTable,IKappa,
     !                 nTabMod,iDualPor,NVarW,iQSame,lAw2Ad,lKSame)

      dimension ParD(11,NMat),ParW(11,NMat),Ah(NumNP),AhW(NMat),
     !          AKW(NMat),AThW(NMat),MatNum(NumNP),hNew(NumNP),
     !          Kappa(NumNP),AThS(NumNP),ThRR(NumNP),ConR(NumNP),
     !          AKS(NumNP),KappaO(NumNP)
      logical lScreen,lInverse,lAw2Ad,lKSame,lTable

      lInverse=.false.
      if(NVarW.gt.0) lInverse=.true.
      if(lScreen) write(*,*) 'reading material information'
      read(30,*,err=901)
      read(30,*,err=901) hTab1,hTabN
      read(30,*,err=901)
      read(30,*,err=901) iModel,iHyst
*     iModel = 0: van Genuchten
*              1: modified van Genuchten (Vogel and Cislerova)
*              2: Brooks and Corey
*              3: van Genuchten with air entry value of 2 cm
*              4: log-normal (Kosugi)
*              5: dual-porosity function (Durner)
*              6: dual-porosity system with transfer proportional to water content
*              7: dual-porosity system with transfer proportional to pressure head
*              8: dual-permeability system, not handled by this program
*              9: nTabMod: general tables (nTabMod)
      if(iModel.eq.8) then
        write(*,*) 'Dual-permeability models are implemented in differen
     !t code !!'
        write(*,*) 'Press Enter to continue'
        read(*,*)
        stop
      end if
      if(iModel.lt.nTabMod) then
        hTab1=-amin1(abs(hTab1),abs(hTabN))
        hTabN=-amax1(abs(hTab1),abs(hTabN))
        lTable=.true.
        if((hTab1.gt.-0.00001.and.hTabN.gt.-0.00001).or.hTab1.eq.hTabN)
     !                                                         then
          lTable=.false.
          hTab1=-0.0001*xConv
          hTabN=-100.  *xConv
        end if
      else
        lTable=.true.
      end if

      if(iHyst.gt.0) then
        read(30,*,err=901)
        read(30,*,err=901) IKappa
      else
        IKappa=-1
      endif
      do 11 i=1,NumNP
        Kappa(i)=IKappa
        KappaO(i)=IKappa
11    continue
      if(iModel.eq.2.or.iModel.eq.4.or.
     !   (iModel.eq.0.and.iHyst.eq.0)) then
        write(50,110,err=902)
      else if(iModel.eq.1.or.iModel.eq.3) then
        write(50,111,err=902)
      else if(iModel.eq.0.and.iHyst.gt.0) then
        write(50,112,err=902)
      else if (iModel.eq.5) then
        write(50,113,err=902)
      else if (iModel.eq.6) then
        write(50,115,err=902)
      else if (iModel.eq.7) then
        write(50,116,err=902)
      else if (iModel.eq.nTabMod) then
        write(50,114,err=902)
      end if
      read(30,*,err=901)
      rHEntry=0.02*xConv
      if(iModel.eq.0.or.iModel.eq.2.or.iModel.eq.3.or.iModel.eq.4) then
        NPar=6
      else if(iModel.eq.1) then
        NPar=10
      else if(iModel.eq.5) then
        NPar=9
      else if(iModel.eq.6) then
        NPar=9
        iModel=0
        iDualPor=1
      else if(iModel.eq.7) then
        NPar=11
        iModel=0
        iDualPor=2
      else if(iModel.eq.nTabMod) then
        NPar=3
      else
        NPar=6
      end if
      do 12 M=1,NMat
        if(iHyst.eq.0) then
          if(lInverse) then
            read(30,*,err=901)
          else
            read(30,*,err=901) (ParD(i,M),i=1,NPar)
          end if
          if(iModel.eq.1) then
            ParD(7,M)=amax1(ParD(7,M),ParD(2,M))
            ParD(8,M)=amin1(ParD(8,M),ParD(1,M))
            ParD(9,M)=amin1(ParD(9,M),ParD(2,M))
            ParD(10,M)=amin1(ParD(10,M),ParD(5,M))
          else if(iModel.eq.3) then
            ParD(7,M)=ParD(1,M)+(ParD(2,M)-ParD(1,M))*
     !            (1.+(ParD(3,M)*rHEntry)**ParD(4,M))**(1.-1./ParD(4,M))
          end if
          if(iModel.eq.nTabMod) then
            write(50,121,err=902) M,(ParD(i,M),i=1,NPar)
          else
            write(50,120,err=902) M,(ParD(i,M),i=1,NPar)
          end if
        else  ! Hysteresis
          if(lInverse) then
            read(30,*,err=901)
            ParW(2,M)=ParD(8,M)
            ParW(3,M)=ParD(9,M)
            ParW(5,M)=ParD(10,M)
            if(iQSame.gt.0) then
              if(iQSame.eq.2) then
                ParW(2,M)=ParD(2,M)
              else
                ParW(2,M)=amin1(ParD(2,M),ParW(2,M))
              end if
              ParD(7,M)=ParD(2,M)
              ParD(9,M)=ParD(2,M)
              ParW(7,M)=ParW(2,M)
              ParW(9,M)=ParW(2,M)
            else
              ParW(2,M)=amin1(ParD(2,M),ParW(2,M))
              ParW(7,M)=amin1(ParD(7,M),ParW(7,M))
              ParW(9,M)=amin1(ParD(9,M),ParW(9,M))
            end if
            if(lAw2Ad) then
              ParW(3,M)=2.*ParD(3,M)
            else
              ParW(3,M)=amax1(ParD(3,M),ParW(3,M))
            end if
            if(lKSame) then
              ParD(10,M)=ParD(5,M)
              ParW(5,M) =ParD(5,M)
              ParW(10,M)=ParD(5,M)
            else
              ParW(5,M) =amin1(ParW(5,M), ParD(5,M))
              ParW(10,M)=amin1(ParW(10,M),ParD(10,M))
            end if
          else
            read(30,*,err=901) (ParD(i,M),i=1,7),ParW(2,M),ParW(3,M),
     !                          ParW(5,M)
          end if
          ParD(7,M)=amax1(ParD(7,M),ParD(2,M))
          write(50,120,err=902) M,(ParD(i,M),i=1,7),ParW(2,M),
     !                          ParW(3,M),ParW(5,M)
          ParW(1,M)=ParD(1,M)
          ParW(4,M)=ParD(4,M)
          AhW(M)=ParD(3,M)/ParW(3,M)
          AThW(M)=(ParW(2,M)-ParW(1,M))/(ParD(2,M)-ParD(1,M))
          AKW(M)=1.0
          if(iHyst.eq.2) AKW(M)=ParW(5,M)/ParD(5,M)
          ParW(7,M)=ParW(1,M)+AThW(M)*(ParD(7,M)-ParD(1,M))
          ParD(8,M)=ParD(1,M)
          ParD(9,M)=ParD(2,M)
          ParD(10,M)=ParD(5,M)
          ParW(8,M)=ParW(1,M)
          ParW(9,M)=ParW(2,M)
          ParW(10,M)=ParW(5,M)
          ParW(6,M)=ParD(6,M)
        end if
12    continue

*     Hysteresis Update for Initial Pressure Head Distributions

      do 13 i=1,NumNP
        M=MatNum(i)
        if(iModel.lt.nTabMod)
     !    hNew(i)=amax1(hNew(i),Ah(i)*FH(iModel,0.00000001,ParD(1,M)))
        AThS(i)=1.
        AKS(i)=1.
        ThRR(i)=ParD(1,M)
        ConR(i)=0.
13    continue
      return

*     Calculate the air-water interfacial area
      iTab=101
      g=9.81*xConv            ! Gravitational acceleration from [m/s2] to [L/s2]
      Temp=20.
      Sigma=75.6-0.1425*Temp-2.38e-4*Temp**2           ! Surface tension [g/s2]
      row=1.-7.37e-6*(Temp-4.)**2+3.79e-8*(Temp-4.)**3 ! Density of soil water [g/cm3]
      row=row*1.e+06/xConv/xConv/xConv                 ! to [g/L3]

      do 16 m=1,NMat
        Ae=0.
        Se=1.
        hh=0.
        write(50,117)
        write(50,122) Se,ParD(2,m),hh,Ae
        do 14 i=2,iTab-1   ! Se from 0.99 to 0.01
          Se=1.-(i-1)*0.01
          if(Ae.lt.1.e+6) call qromb(1.,Se,sInt,iModel,ParD(1,m))
          Ae=sInt*ParD(2,M)*row*g/sigma
          Th=ParD(1,m)+Se*(ParD(2,m)-ParD(1,m))
          hh=FH(iModel,Th,ParD(1,M))
          write(50,122) Se,Th,hh,Ae
14      continue
        do 15 i=1,10   ! Se from 0.009 to 0.000
          Se=1.-(iTab-2)*0.01-i*0.001
          if(Ae.lt.1.e+6) call qromb(1.,Se,sInt,iModel,ParD(1,m))
          Ae=sInt*ParD(2,M)*row*g/sigma
          Th=ParD(1,m)+Se*(ParD(2,m)-ParD(1,m))
          hh=FH(iModel,Th,ParD(1,M))
          write(50,122) Se,Th,hh,Ae
15      continue
16    continue

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l'/)
111   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l       Qm     Qa     Qk       Kk
     !'/)
112   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l       Qm     QsW  AlfaW      Ks
     !W'/)
113   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l       W2       Alfa2        n2'
     !/)
114   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        Ks'
     !/)
115   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l      QrIm   QsIm  Omega'/)
116   format(//'MatNum, Param. array:'//'   Mat     Qr     Qs        ',
     !'Alfa         n          Ks      l      QrIm   QsIm   Alfa2       
     ! n2  Omega'/)
117   format(//'   Saturation       Theta   Pressure Head  Air-Water Are
     !a'/)
120   format(i5,2x,2f7.3,4e12.3,3f7.3,2e12.3)
121   format(i5,2x,2f7.3,3e12.3,2f7.3,3e12.3)
122   format(2f13.5,2e15.6,f13.5)
      end

************************************************************************

      subroutine HysterIn(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,
     !                    Kappa,AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,
     !                    iHyst,iModel,cDataPath)

      character cFileName*260,cDataPath*260
      dimension MatNum(NumNP),ThOld(NumNP),hOld(NumNP),ParD(11,NMat),
     !          ThNew(NumNP),Kappa(NumNP),AThS(NumNP),ThRR(NUmNP),
     !          ConO(NumNP),ConR(NumNP),AKS(NumNP),KappaO(NumNP),
     !          Ah(NumNP),AK(NumNP),ParW(11,NMat)

      iLengthPath = Len_Trim(cDataPath)
      cFileName = cDataPath(1:iLengthPath)//'/HYSTERESIS.IN'
      open(35,file=cFileName,status='old',err=901)
      read(35,*,err=901) ! Input file "Hyster.in"
      read(35,*,err=901) ! n,Theta,Kappa
      do 11 i=1,NumNP
        read(35,*,err=901) n,ThOld(i),KappaO(i)
        ThNew(i)=ThOld(i)
        Kappa(i)=KappaO(i)
11    continue
      call Hyster(NumNP,NMat,hOld,MatNum,ParD,ParW,ThNew,ThOld,
     !            Kappa,AThS,ThRR,ConO,ConR,AKS,KappaO,Ah,AK,
     !            iHyst,iModel,-1.)
901   continue

      return
      end
      
************************************************************************

      subroutine GenMat(NTab,NTabD,NMat,thr,ths,hSat,Par,hTab,ConTab,
     !                  CapTab,ConSat,TheTab,iModel,lScreen,nTabMod,
     !                  ConSMax,xConv,tConv,ierr)

      dimension hTab(NTabD,NMat),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          Par(11,NMat),ConSat(NMat),thr(NMat),hSat(NMat),
     !          ths(NMat),TheTab(NTabD,NMat),NTab(NMat)
      logical lScreen

      if(lScreen) write(*,*) 'generating materials'
      if(iModel.lt.nTabMod) then
        write(50,110,err=901)
        write(50,120,err=901)
        hTab1=hTab(1,1)
        hTabN=hTab(NTab(1),1)
        dlh=(alog10(-hTabN)-alog10(-hTab1))/(NTab(1)-1)
        do 11 i=1,NTab(1)
          alh=alog10(-hTab1)+(i-1)*dlh
          hTab(i,1)=-10**alh
11     continue
      else
        read(36,*,err=901)
        read(36,*,err=901) iCap !(=1; input capacity; =0: do not input capacity)
        write(50,111,err=901)
        write(50,120,err=901)
      end if
      ConSMax=0.
c      ConSMax=1e+30
      do 13 M=1,NMat
        if(iModel.lt.nTabMod) then
          hSat(M)  =FH(iModel,1.0,Par(1,M))
          ConSat(M)=Par(5,M)
          if(ConSat(M).gt.ConSMax) ConSMax=ConSat(M)
c          if(ConSat(M).lt.ConSMax) ConSMax=ConSat(M)
          thr(M)   =Par(1,M)
          ths(M)   =Par(2,M)
          write(50,*,err=901)
          do 12 i=1,NTab(1)
            ConTab(i,M)=FK(iModel,hTab(i,1),Par(1,M))
            CapTab(i,M)=FC(iModel,hTab(i,1),Par(1,M))
            TheTab(i,M)=FQ(iModel,hTab(i,1),Par(1,M))
            Qe         =FS(iModel,hTab(i,1),Par(1,M))
            ConV=ConVh(hTab(i,1),TheTab(i,M),ths(M),xConv,tConv)
            a10h=alog10(max(-hTab(i,1),1e-30))
            a10K=alog10(ConTab(i,M))
            write(50,130,err=901) TheTab(i,M),hTab(i,1),a10h,
     !                            CapTab(i,M),ConTab(i,M),a10K,Qe,ConV
12        continue
        else if(iModel.eq.nTabMod) then ! Table
          read(36,*)
          read(36,*) NTab(M)
          if(NTab(M).gt.NTabD) goto 902
          read(36,*)
          hSat(M)  =0.0
          ConSat(M)=Par(3,M)
          thr(M)   =Par(1,M)
          ths(M)   =Par(2,M)
          do 15 i=1,NTab(M)
            if(iCap.eq.1) then
              read(36,*) TheTab(i,M),hTab(i,M),ConTab(i,M),CapTab(i,M)
            else
              read(36,*) TheTab(i,M),hTab(i,M),ConTab(i,M)
            end if
15        continue
          write(50,*,err=901)
          do 16 i=1,NTab(M)
            if(iCap.eq.0) then
              if(i.eq.1) then
                CapTab(i,M)=(TheTab(2,M)-TheTab(1,M))/
     !                      (  hTab(2,M)-  hTab(1,M))
              else if(i.eq.NTab(M)) then
                CapTab(i,M)=(TheTab(NTab(M),M)-TheTab(NTab(M)-1,M))/
     !                      (  hTab(NTab(M),M)-  hTab(NTab(M)-1,M))
              else
                CapTab(i,M)=(TheTab(i+1,M)-TheTab(i-1,M))/
     !                      (  hTab(i+1,M)-  hTab(i-1,M))
              end if
            end if
            Qe  =(TheTab(i,M)-Par(1,M))/(Par(2,M)-Par(1,M))
            a10h=alog10(max(-hTab(i,M),1e-30))
            a10K=alog10(ConTab(i,M))
            write(50,130,err=901) TheTab(i,M),hTab(i,M),a10h,
     !                            CapTab(i,M),ConTab(i,M),a10K,Qe
16        continue
        end if
        write(50,140,err=901)
13    continue
      return

*     Error when writing into an output file
901   ierr=1
      return

*     The dimensions in Mater.in are larger than NTabD
902   ierr=2
      return

110   format(/7x,'Table of Hydraulic Properties which are interpolated i
     !n simulation'/7x,65('=')/)
111   format(/7x,'Hydraulic Properties which are interpolated from input 
     ! tables in simulation'/7x,75('=')/)
120   format('  theta         h        log h        C             K',
     !'        log K          S          Kv')
130   format(f8.4,e12.3,e12.4,e12.4,e12.4,e12.4,f10.4,e12.4)
140   format('end')
      end

************************************************************************

      subroutine TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !                TPrint,t,dtOpt,TopInF,BotInF,lScreen,ItMin,
     !                ItMax,MaxAL,hCritS,NPD,AtmBC,iVer,lPrintD,nPrStep,
     !                tPrintInt,lEnter,lDayVar,lSinPrec,lLAI,rExtinct,
     !                nCycles,lChem,lYield,CTp,CRp,iInterc,aInterc,ierr)

      logical TopInF,BotInF,AtmBC,lPrintD,lEnter,lScreen,lDayVar,lLAI,
     !        lSinPrec,lBCCycles,lChem,lYield,lInterc
      double precision t,tInit,tOld,tMax,tAtm,TPrint,tPrintInt
      dimension TPrint(NPD)

      if(lScreen) write(*,*) 'reading time information'
      read(30,*,err=901)
      read(30,*,err=901)
      read(30,*,err=901) dt,dtMin,dtMax,dMul,dMul2,ItMin,ItMax,MPL
      if(MPL.gt.NPD) then
        ierr=2
        return
      end if
      read(30,*,err=901)
      read(30,*,err=901) tInit,tMax
      if(iVer.gt.2) then
        read(30,*,err=901)
        read(30,*,err=901) lPrintD,nPrStep,tPrintInt,lEnter
      end if
      read(30,*,err=901)
      read(30,*,err=901) (TPrint(i),i=1,MPL)
      if ((MPL>0).and.(TPrint(1)<=dt)) then ! modified by Adam Szymkiewicz
        dt = TPrint(1)/2.
      end if
      dtOpt=dt
      if(TopInF.or.BotInF.or.AtmBC) then
        nRead=1
        if(lYield) nRead=2
        do 12 iRead=1,nRead
          iVerA = iGetFileVersion(31,1)
          read(31,*,err=901)
          read(31,*,err=901)
          read(31,*,err=901) MaxAL
          if(iVerA.ge.4) then
            read(31,*,err=901)
            read(31,*,err=901) lDayVar,lSinPrec,lLAI,lBCCycles,lInterc
            if(lInterc) iInterc=1
            if(lLAI) then
              read(31,*,err=901)
              read(31,*,err=901) rExtinct
            end if
            if(lBCCycles) then
              read(31,*,err=901)
              read(31,*,err=901) nCycles
              nCycles=max0(1,nCycles)
            end if
            if(lInterc) then
              read(31,*,err=901)
              read(31,*,err=901) aInterc
            end if
          end if
          read(31,*,err=901) 
          read(31,*,err=901) hCritS
          read(31,*,err=901)
          if(lYield.and.iRead.eq.1) then
            CTp=0.   ! Cumulative potential transpiration
            CRp=0.   ! Cumulative potential nutrient demand
            tPrev=tInit
            do 11 iAtm=1,MaxAL
*             reading setup for two solutes, the second being nutrient
              if(lChem) then
                read(31,*,err=901) tAtm,Prec,rSoil,rR,hCA,rB,hB,hT,tTop,
     !                             tBot,Ampl,cT1,cB1,cT2,cB2
              else
                read(31,*,err=901) tAtm,Prec,rSoil,rR
              end if
              dtAtm=tAtm-tPrev
              CTp=CTp+dtAtm*rR
              if(lChem) CRp=CRp+dtAtm*cB2
              tPrev=tAtm
11          continue
            rewind(31)
          end if
12      continue
      else
        tAtm=tMax
      end if
      TPrint(MPL+1)=tMax
      tOld=tInit
      t=tInit+dt
      return

*     Error when reading from an input file
901   ierr=1
      return
      end

************************************************************************

      subroutine MeteoIn(Latitude,Altitude,ShortWaveRadA,ShortWaveRadB,
     !                   LongWaveRadA,LongWaveRadB,LongWaveRadA1,
     !                   LongWaveRadB1,WindHeight,TempHeight,iCrop,iLAI,
     !                   CropHeight,Albedo,LAI,xRoot,iInterc,aInterc,
     !                   iGrowth,rGrowth,rExtinct,iRadiation,lEnBal,
     !                   lPrint,iSunSh,iRelHum,CloudF_Ac,CloudF_Bc,
     !                   lHargr,lMetDaily,xConv,ierr)

      implicit real(A-H,L-Z)
      logical lEnBal,lPrint,lHargr,lMetDaily
      dimension rGrowth(1000,5)

*     Latitude             ! Latitude of the location, [degree, N=+, S=-]
*     Altitude             ! Altitude of the location above mean see level [m]
*     ShortWaveRadA=0.25   ! fraction of extraterrestrial readiation on 
*                          ! overcast days, first Angstrom coefficient, a_s, eq.55
*     ShortWaveRadB=0.50   ! Input, fraction of extraterrestrial readiation on 
*                          !   overcast days, second Angstrom coefficient, b_s, eq.55
*     LAI                  ! Leaf area index [-]
*     LongWaveRadA=0.90    ! cloudiness factor, a_c, eq. 59
*     LongWaveRadB=0.10    ! cloudiness factor, b_c, eq. 59
*     LongWaveRadA1=0.34   ! emissivity correlation coefficient, a_l, eq. 60
*     LongWaveRadB1=-0.139 ! emissivity correlation coefficient, b_l, eq. 60
*     WindHeight=200.      ! Measurement hight of wind, [cm]
*     TempHeight=190.      ! Measurement hight of temperature, [cm]

*     lEnBal               ! Calculate energy balance
*     lHargr               ! Use Hargreaves equation to calculate potential ET
*     lMetDaily            ! Generate daily variations of meteo variables from daily values
*     iRadiation           ! =0: Potential radiation; =1: solar radiation; =2: net radiation
*     iCrop                ! =1: constant; =2: read from a table; =3: read daily growth data
*     iSunSh               ! =0: Sunshine hours; =1: cloudiness; =2: transmission coefficient
*     iRelHum	             ! =0: Relative humidity; =1: Vapor pressure

      iVerM = iGetFileVersion(33,1)
      read(33,*,err=901)
      read(33,*,err=901)
      read(33,*,err=901) MaxALMet,iRadiation,lHargr
      if(iVerM.ge.4) then
        read(33,*,err=901)
        read(33,*,err=901) lEnBal,lMetDaily
      end if
      Altitude=0.
      if(iRadiation.ne.2) then
        read(33,*,err=901)
        read(33,*,err=901) Latitude,Altitude
        read(33,*,err=901)
        read(33,*,err=901) ShortWaveRadA,ShortWaveRadB
        read(33,*,err=901)
        read(33,*,err=901) LongWaveRadA,LongWaveRadB
        read(33,*,err=901)
        read(33,*,err=901) LongWaveRadA1,LongWaveRadB1
      end if
      read(33,*,err=901)
      read(33,*,err=901) WindHeight,TempHeight
      read(33,*,err=901)
      read(33,*,err=901) iCrop,iSunSh,iRelHum
      if(iRadiation.eq.1.and.iSunSh.eq.3) then
        read(33,*,err=901)
        read(33,*,err=901) CloudF_Ac,CloudF_Bc
      end if
      if(iCrop.ge.1) then
        read(33,*,err=901)
        read(33,*,err=901) iLAI,rExtinct
        read(33,*,err=901)
        read(33,*,err=901) iInterc
        if(iCrop.eq.1) then
          read(33,*,err=901)
          read(33,*,err=901) CropHeight,Albedo,LAI,xRoot
          CropHeight=CropHeight*100./xConv ! conversion to cm
        else if (iCrop.eq.2) then
          read(33,*,err=901)
          read(33,*,err=901) iGrowth
          if(iGrowth.gt.1000) then
            write(*,*) 'Number of crop growth data is larger than 1000'
            write(*,*) 'Press Enter to continue'
            read(*,*)
            stop
          end if
          read(33,*,err=901)
          do 11 i=1,iGrowth
            read(33,*,err=901) (rGrowth(i,j),j=1,5)
            rGrowth(i,2)=rGrowth(i,2)*100./xConv ! conversion to cm
11        continue
        end if
        if(iInterc.eq.1) then
          read(33,*,err=901)
          read(33,*,err=901) aInterc
        end if
      else
        read(33,*,err=901)
        read(33,*,err=901) Albedo
        iLai=0.
        iInterc=0.
      end if
      read(33,*,err=901)
      read(33,*,err=901)
      read(33,*,err=901)
      if(lPrint) then
        if(.not.lMetDaily) then
          if(lEnBal) then
            write(43,110)
          else
            write(43,120)
          end if
        else
          if(lEnBal) then
            write(43,130)
          else
            write(43,140)
          end if
        end if
      end if

      return

*     Error when reading from an input file
901   ierr=1
      return

110   format(/'     Time      Short.Rad.    Long.Rad.    Radiation    Se
     !nsible      Latent     Heat Flux     Balance'/
     !        '      [d]       [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]    [M
     !J/m2/d]   [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]'/)
120   format(/'   Time       ET       EvapP   TranspP     Rns        Rnl
     !    RadTerm  AeroTerm     Prec     Interc   ExInterc'/
     !        '    [d]     [mm/d]    [mm/d]    [mm/d] [MJ/m2/d]  [MJ/m2/
     !d]   [mm/d]    [mm/d]    [mm/d]    [mm/d]    [mm/d]'/)
130   format(/'     Time      Short.Rad.    Long.Rad.    Radiation    Se
     !nsible      Latent     Heat Flux     Balance        AirTemp       
     !AirRH      SolarRad'/
     !        '      [d]       [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]    [M
     !J/m2/d]   [MJ/m2/d]    [MJ/m2/d]    [MJ/m2/d]         [C]         
     ! [%]       [MJ/m2/d]'/)
140   format(/'   Time       ET        Evap    Transp     Rns        Rnl
     !    RadTerm  AeroTerm     Prec     Interc   ExInterc   AirTemp    
     !AirRH   SolarRad'/
     !        '    [d]     [mm/d]    [mm/d]    [mm/d] [MJ/m2/d]  [MJ/m2/
     !d]   [mm/d]    [mm/d]    [mm/d]    [mm/d]    [mm/d]      [C]      
     ![%]    [MJ/m2/d]'/)
      end

************************************************************************

      subroutine SinkIn(NMat,lChem,lMoSink,lSolRed,lSolAdd,P0,POptm,
     !                  P2H,P2L,P3,r2H,r2L,aOsm,c50,P3c,NS,lMsSink,
     !                  cRootMax,iVer,OmegaC,lActRSU,OmegaS,SPot,rKM,
     !                  cMin,lOmegaW,iSolAU,lHanks,lDblVeg,cDataPath,
     !                  rA,rB,rC,rD,lScreen,ierr)

      character cFileName*260,cDataPath*260
      dimension POptm(NMat),aOsm(NS),cRootMax(NS)
      logical lChem,lMoSink,lSolRed,lSolAdd,lScreen,lMsSink,lActRSU,
     !        lOmegaW,lHanks,lDblVeg
      common /DoubleVeg/ xP0,xP2H,xP2L,xP3,xr2H,xr2L,xPOptm

      if(lScreen) write(*,*) 'reading sink information'
      read(30,*,err=901)
      read(30,*,err=901)
      if(iVer.le.2) then
        read(30,*,err=901) iMoSink,(cRootMax(i),i=1,NS)
      else
        read(30,*,err=901) iMoSink,(cRootMax(i),i=1,NS),OmegaC
      end if
      if(iMoSink.eq.0) then
        lMoSink=.true.
      else
        lMoSink=.false.
      end if
      read(30,*,err=901)
      if(lMoSink) then
        read(30,*,err=901) P0,P2H,P2L,P3,r2H,r2L
        read(30,*,err=901)
        read(30,*,err=901) (POptm(i),i=1,NMat)
        P0 =-abs(P0)
        P2L=-abs(P2L)
        P2H=-abs(P2H)
        P3 =-abs(P3)
        if(lDblVeg) then
          read(35,*,err=903)
          read(35,*,err=903)
          read(35,*,err=903) xP0,xP2H,xP2L,xP3,xr2H,xr2L
          read(35,*,err=903)
          read(35,*,err=903) xPOptm
          read(35,*,err=903)
          read(35,*,err=903)
          xP0 =-abs(xP0)
          xP2L=-abs(xP2L)
          xP2H=-abs(xP2H)
          xP3 =-abs(xP3)
          xPOptm=-abs(xPOptm)
          goto 904
903       lDblVeg=.false.
904       continue
        end if
      else
        read(30,*,err=901) P0,P3
      end if
      if(lChem) then
        read(30,*,err=901)
        read(30,*,err=901) lSolRed
        if(lSolRed) then
          read(30,*,err=901)
          read(30,*,err=901) lSolAdd
          read(30,*,err=901)
          if(lSolAdd) then
            read(30,*,err=901) (aOsm(i),i=1,NS)
          else
            read(30,*,err=901) c50,P3c,(aOsm(i),i=1,NS),iMsSink
            if(iMsSink.eq.0) then
              lMsSink=.false.
            else
              lMsSink=.true.
            end if
          end if
        end if
        if(lActRSU) then              ! Active uptake only for a selected solute
          read(30,*,err=901)
          read(30,*,err=901) OmegaS,SPot,rKM,cMin,lOmegaW,iSolAU
          rA=1.
          rB=0.
          rC=0.
          rD=0.
c          read(30,*,err=901)               !Naftali
c          read(30,*,err=901) rA,rB,rC,rD   !Naftali
        end if
      end if

      if(lHanks) then
        iLengthPath = Len_Trim(cDataPath)
        cFileName = cDataPath(1:iLengthPath)//'/ROOT_H.OUT'
        open(62,file=cFileName,status='unknown',err=901)
        write(62,110)
110     format(/'         Time       hRoot')
      end if      
      return

*     Error when reading from an input file
901   ierr=1
      return
      end

************************************************************************

      subroutine RootIn(tRMin,tRHarv,xRMin,xRMax,RGR,lScreen,iver,
     !                  iRootIn,nGrowth,rGrowth,tRPeriod,ierr)

      logical lScreen
      dimension rGrowth(1000,5)

      tRPeriod=1.e+30
      if(lScreen) write(*,*) 'reading of root growth information'
      read(30,*,err=901)
      if(iVer.ge.4) then
        read(30,*,err=901)
        read(30,*,err=901) iRootIn 
*           =1: Table; =2: Logistic growth; =0: With atmospheric data
        if(iRootIn.eq.1) then ! Table
          read(30,*,err=901)
          read(30,*,err=901) nGrowth
          if(nGrowth.gt.1000) then
            write(*,*) 'Number of crop growth data is larger than 1000'
            write(*,*) 'Press Enter to continue'
            read(*,*)
            stop
          end if
          read(30,*,err=901)
          do 11 i=1,nGrowth
            read(30,*,err=901) rGrowth(i,1),rGrowth(i,5)
11        continue
        end if
      end if
      if(iVer.lt.4.or.iRootIn.eq.2) then  ! Logistic growth
        iRootIn=2
        read(30,*,err=901)
        if(iVer.lt.4) then
          read(30,*,err=901) iRFak,tRMin,tRMed,tRHarv,xRMin,xRMed,xRMax
        else
          read(30,*,err=901) iRFak,tRMin,tRMed,tRHarv,xRMin,xRMed,xRMax,
     !                       tRPeriod
        end if
        if(iRFak.eq.1) then
          tRMed=(tRHarv+tRMin)/2.
          xRMed=(xRMax+xRMin)/2.
        end if
        rtm=tRMed-tRMin
        if(rtm.lt.1.e-20.or.xRMed.lt.1.e-10) then
          write(*,*) 'Time(depth) Root Data must be larger then Initial 
     !Root Growth Time(depth) !!'
          goto 901
        end if
        RGR=0.
        if(xRMed*(xRMax-xRMin).gt.0.)
     !    RGR=-(1./rtm)*alog(amax1(.0001,(xRMin*(xRMax-xRMed)))/
     !                                   (xRMed*(xRMax-xRMin)))
        write(50,110,err=902) tRMin,tRHarv,xRMin,xRMax,RGR
      end if
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(/' Root growth information'/1x,23('=')/' tRMin = ',f10.3,
     !' tRHarv = ',f10.3/' xRMin = ',f10.3,' xRMax = ',f10.3/
     !' Root growth rate = ',e11.3)
      end

************************************************************************

*     Reading root distribution information (when lRootIn=true)

      subroutine BetaIn(NumNP,Beta,z,cDataPath,lScreen)

      logical lScreen
      character cFileName*260,cDataPath*260
      dimension Beta(NumNP),z(NumNP)

      if(lScreen) write(*,*) 'reading of root distribution information'

      iLengthPath = Len_Trim(cDataPath)
      cFileName = cDataPath(1:iLengthPath)//'/ROOT.IN'

      open(35,file=cFileName,status='old',err=901)
      read(35,*,err=901) ! Input file "Root.in"
      read(35,*,err=901) ! iRShape
      read(35,*,err=901) iRShape  !=0: constant, =1: exponential, 
                                  !=2: linear,   =3: VG&H,        =4: Vrugt
      read(35,*,err=901) ! Max rooting depth
      read(35,*,err=901) rZmax
      read(35,*,err=901) ! Parameters
      if     (iRShape.eq.1) then
        read(35,*,err=901) alpha
      else if(iRShape.eq.4) then
        read(35,*,err=901) rZ0,rA
      end if

      do 11 i=1,NumNP
        Beta(i)=0.
11    continue
      if(rZmax.le.0.) return  ! no root zone

      zSurf=z(NumNP)
      do 12 i=1,NumNP
        if(z(i).ge.-rZmax) then
          if     (iRShape.eq.0) then
            Beta(i)=1.
          else if(iRShape.eq.1) then
            Beta(i)=exp(alpha*z(i))
          else if(iRShape.eq.2) then
            Beta(i)=1.+z(i)/rZmax
          else if(iRShape.eq.3) then
            if(z(i).lt.zSurf-rZmax) then
              Beta(i)=0.
            else if(z(i).lt.zSurf-0.2*rZmax) then
              Beta(i)=2.08333/rZmax*(1.-(zSurf-z(i))/rZmax)
            else
              Beta(i)=1.66667/rZmax
            end if
          else if(iRShape.eq.4) then
            r1=(rZmax-(zSurf-z(i)))/rZmax
            rRootA=rA
            if((zSurf-z(i)).gt.rZ0) rRootA=1.
            r2=rRootA/(rZmax)*abs(rZ0-(zSurf-z(i)))
            r2=exp(-r2)
            Betai=amax1(r1*r2,0.)
            Beta(i)=Beta(i)+Betai
          end if
        end if
12    continue


*     Normalize Beta Function
      SBeta=0.
      if(Beta(NumNP).gt.0.) SBeta=Beta(NumNP)*(z(NumNP)-z(NumNP-1))/2.
      do 19 i=2,NumNP-1
        if(Beta(i).gt.0.) SBeta=SBeta+Beta(i)*(z(i+1)-z(i-1))/2.
19    continue
      do 20 i=2,NumNP
        if(SBeta.gt.0.) then
          if(Beta(i).gt.0.) Beta(i)=Beta(i)/SBeta
        else
          Beta(i)=0.
        end if
20    continue
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return
      end

************************************************************************

      subroutine TempIn(NMat,TPar,Ampl,tPeriod,kTopT,tTop,kBotT,tBot,
     !                  TopInf,BotInf,iCampbell,iVer,lSnow,SnowMF,
     !                  SnowSF,SnowLayer,lScreen,ierr,NVarH)

      logical TopInf,BotInf,lScreen,lInverse,lSnow,lSnowInput
      dimension TPar(10,NMat)

      lInverse=.false.
      if(NVarH.gt.0) lInverse=.true.
      if(lScreen) write(*,*) 'reading heat transport information'
      read(30,*,err=901)
      read(30,*,err=901)
      do 11 i=1,NMat
        if(lInverse) then
          read(30,*,err=901)
        else
          read(30,*,err=901)  (TPar(j,i),j=1,9)
        end if
11    continue
      read(30,*,err=901)
      if(iVer.le.2) then
        read(30,*,err=901) Ampl,tPeriod
        iCampbell=0
      else if(iVer.eq.3) then
        read(30,*,err=901) Ampl,tPeriod,iCampbell,SnowMF
      else
        read(30,*,err=901) Ampl,tPeriod,iCampbell,SnowMF,lSnowInput
        if(lSnow.and.lSnowInput) then
          read(30,*,err=901)
          read(30,*,err=901) SnowMF,SnowSF,SnowLayer
        end if
      end if
      write(50,110,err=902) Ampl
      do 12 i=1,NMat
        write(50,120,err=902) (TPar(j,i),j=1,9)
12    continue
      read(30,*,err=901)
      read(30,*,err=901) kTopT,tT,kBotT,tB
      if(.not.TopInf) tTop=tT
      if(.not.BotInf) tBot=tB
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(//' Heat transport information'/1x,26('=')//' ample = ',
     ! f10.3//'   Beta    Qn     Qo        B1         B2         B3
     !    Cn         Co         Cw')
120   format(3f7.3,6e11.3)
      end

************************************************************************

      subroutine Profil(N,NMat,x,MatNum,xSurf,Beta,Ah,AK,ATh,thr,thS,
     !                  ConS,hS,lScreen,ierr)

      dimension x(N),MatNum(N),Beta(N),Ah(N),AK(N),ATh(N),
     !          thr(NMat),thS(NMat),ConS(NMat),hS(NMat)
      logical lScreen

      if(lScreen) write(*,*)'printing profile information'
      write(78,110,err=901)
      ConSN=ConS(MatNum(N))*AK(N)
      do 11 i=N,1,-1
        M=MatNum(i)
        write(78,120,err=901) N-i+1,xSurf-x(i),thr(M),thr(M)+ATh(i)*
     !                        (thS(M)-thr(M)),hS(M)*Ah(i),ConS(M)*AK(i),
     !                        ConS(M)*AK(i)/ConSN,Beta(i),Ah(i),
     !                        AK(i),ATh(i)
11    continue
      write(78,'(''end'')',err=901)
      return

*     Error when writing into an output file
901   ierr=1
      return

110   format(//'    n      depth     THr       THs       hs       Ks',
     !'        Ks/KsTop     Beta      Ah        AK        ATh'/)
120   format(i5,f10.2,2f10.3,f10.1,e12.3,5f10.3)
      end

************************************************************************

      subroutine IrrigIn(iIrrig,hIrrig,rIrrig,tIrrig,iIrrigTr,nIrreg,
     !                   ierr)

      dimension tIrrig(4)

c      if(lIrrig) then
c        iLengthPath = Len_Trim(cDataPath)
c        cFileName = cDataPath(1:iLengthPath)//'/IRRIG.IN'
c        open(13,file=cFileName,status='old',err=901)
c        read(13,*,err=901)
c        read(13,*,err=901) lIrrig
      iIrrigTr=1
      nIrreg=0
      read(30,*,err=901)
      read(30,*,err=901) iIrrig,hIrrig,rIrrig,tIrrig(3),tIrrig(4)
      if(tIrrig(3).gt.0.) iIrrigTr=0
      write(14,110) iIrrig,hIrrig,rIrrig
      if(iIrrigTr.eq.0) write(14,120) tIrrig(3)
      if(iIrrigTr.eq.1) write(14,130)
      write(14,140) tIrrig(4)
      return

*     Error when reading from an input file 
901   ierr=1
      return

110   format(/' Irrigation Scheduling:'/
     !        ' ======================'/
     !        ' Observation node triggering irrigation: ', i12/
     !        ' Pressure head triggering irrigation [L]:',f12.5/
     !        ' Irrigation rate [L/T]:                  ',e12.5)
120   format( ' Duration of irrigation [T]:             ',f12.5)
130   format( ' Duration of irrigation calculated from transpiration')
140   format( ' Offset of irrigation [T]:               ',f12.5)
      end

************************************************************************

*     Input describing additional fumigant application
*     lAddF - additional fumigant application
*     tAddF - time of additional fumigant application
*     zAddFT - top of additional fumigant application
*     zAddFB - bottom of additional fumigant application
*     FumMass - mass of applied fumigant

      subroutine FumigIn(cDataPath,lFumig,lAddFum,tAddFum,FumMass,
     !                   zAddFT,zAddFB,tInit,tPrev,CumChPrev,lPrint,
     !                   NS,ierr)

      logical lFumig,lAddFum,lPrint
      double precision tAddFum,tInit
      character cFileName*260,cDataPath*260
      dimension CumChPrev(2)

      if(lFumig) then
        iLengthPath = Len_Trim(cDataPath)
        cFileName = cDataPath(1:iLengthPath)//'/FUMIGANTS.IN'
        open(15,file=cFileName,status='old',err=901)
        read(15,*,err=901)
        read(15,*,err=901)
        read(15,*,err=901) lFumig
        read(15,*,err=901)
        read(15,*,err=901) lAddFum
        if(lAddFum) then
          read(15,*,err=901)
          read(15,*,err=901) tAddFum,FumMass
          read(15,*,err=901)
          read(15,*,err=901) zAddFT,zAddFB
        end if
        close(15)
        if(lPrint) then
          cFileName = cDataPath(1:iLengthPath)//'/NOD_INF_C1.OUT'
          open(46,file=cFileName, status='unknown',err=902)
          cFileName = cDataPath(1:iLengthPath)//'/FUMIGANTS.OUT'
          open(48,file=cFileName, status='unknown',err=902)
          if(lAddFum) write(48,110) tAddFum,FumMass,zAddFT,zAddFB
          tPrev=sngl(tInit)
          CumChPrev(1)=0.
          CumChPrev(2)=0.
          write(48,120) 
          if(NS.ge.2) then
            cFileName = cDataPath(1:iLengthPath)//'/NOD_INF_C2.OUT'
            open(47,file=cFileName, status='unknown',err=902)
            cFileName = cDataPath(1:iLengthPath)//'/FUMIGANTS2.OUT'
            open(49,file=cFileName, status='unknown',err=902)
            write(49,120) 
          end if
        end if
      end if
      return

*     Error when reading from an input file 
901   lFumig=.false.
      return
*     Error when writing into an output file
902   ierr=2
      return

110   format(/' Fumigant Application:'/
     !        ' ===================='/
     !        ' Application Time [T]:           ',e12.5/
     !        ' Application Mass [M/L2]:        ',e12.5/
     !        ' Top of application zone [L]:    ',f12.4/
     !        ' Bottom of application zone [L]: ',f12.4/)
120   format(/' Fumigant Fluxes:'/
     !        ' ================'/
     !'      Time       Actual Flux    Cum. Flux  Average flux   ConcG(i
     !Obs)'/
     !'       [T]         [M/L2/T]       [M/L2]     [M/L2/T]     [M/La^3
     !]')
      end

************************************************************************

*     Read information about solute transport

      subroutine ChemIn(lUpW,lTDep,NMat,NS,NSD,MaxItC,ChPar,TDep,kTopCh,
     !                  cTop,kBotCh,cBot,epsi,tPulse,CumCh,cTolA,cTolR,
     !                  lLinear,lEquil,lArtD,PeCr,lScreen,dSurf,cAtm,
     !                  lTort,lMobIm,lBact,lFiltr,iMoistDep,WDep,NMatD,
     !                  iModel,Par,iVer,lDualNEq,lMassIni,lEqInit,iTort,
     !                  ierr,NVarS)

      logical lUpW,lTDep,lLinear(NSD),lEquil,lArtD,lScreen,lInverse,
     !        lTort,lMobIm(NMat),lBact,lFiltr,lMoistDep,lDualNEq,
     !        lMassIni,lEqInit,lVar
      dimension ChPar(NSD*16+4,NMat),TDep(NSD*16+4),cTop(NSD),cBot(NSD),
     !          CumCh(11,NSD),WDep(2+NMatD,NSD*9),Par(11,NMat)

      if(lScreen) write(*,*) 'reading solute transport information'
      lInverse=.false.
      if(NVarS.gt.0) lInverse=.true.
      
      write(50,110,err=902)
      read(30,*,err=901)
      read(30,*,err=901)
      if(iVer.le.2) then
        read(30,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,
     !                     PeCr,NS,lTort
      else
        read(30,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,
     !                     PeCr,NS,lTort,iBact,lFiltr
        lBact=.false.
        if(iBact.eq.1) lBact=.true.
      end if
      if(iVer.ge.4) then
        read(30,*,err=901)
        read(30,*,err=901) iNonEqul,lMoistDep,lDualNEq,lMassIni,lEqInit,
     !                     lVar
        if(lMoistDep) iMoistDep=1
        if(lVar) iTort=1
      end if
      PeCr=amax1(PeCr,0.1)
      if(lUpW) then
        write(50,120,err=902)
      else
        write(50,130,err=902)
        if(lArtD) write(50,140,err=902) PeCr
      end if
      write(50,150,err=902) lTDep,iMoistDep,cTolA,cTolR,MaxItC
      read(30,*,err=901)
      lEquil=.true.
      do 11 M=1,NMat
        if(lInverse) then
          read(30,*,err=901)
        else 
          read(30,*,err=901) (ChPar(j,M),j=1,4)
        end if
        write(50,160,err=902) M,(ChPar(j,M),j=1,4)
        if(ChPar(3,M).lt.1..or.ChPar(4,M).gt.0..or.lBact) lEquil=.false.
        lMobIm(M)=.false.
        if(.not.lBact.and.ChPar(4,M).gt.0.) lMobIm(M)=.true.
        if(.not.lEquil.and.ChPar(1,M).eq.0.) goto 903
11    continue
      do 13 jj=1,NS
        jjj=(jj-1)*16
        write(50,170,err=902) jj
        read(30,*,err=901)
        if(lInverse) then
          read(30,*,err=901)
        else 
          read(30,*,err=901) (ChPar(jjj+j,1),j=5,6)
c          Modification to have diffusion coef. different for different materials
c          Two lines (at the beginning of 12 cycle) need to be commented out.
c          nm=NMat
c          read(30,*,err=901)  ((ChPar(jjj+j,M),j=5,6),M=1,nm)
        end if
        write(50,180,err=902) (ChPar(jjj+j,1),j=5,6)
        read(30,*,err=901)
        lLinear(jj)=.true.
        do 12 M=1,NMat
          ChPar(jjj+5,M)=ChPar(jjj+5,1)
          ChPar(jjj+6,M)=ChPar(jjj+6,1)
          if(lInverse) then 
            read(30,*,err=901)
          else 
            read(30,*,err=901) (ChPar(jjj+j,M),j=7,20)
          end if
          write(50,190,err=902) M,(ChPar(jjj+j,M),j=7,20)
          if     (abs(ChPar(jjj+8,M)-0.0).gt.1.e-12) then
            write(50,200,err=902) M
          else if(abs(ChPar(jjj+9,M)-1.0).gt.0.001) then
            write(50,210,err=902) M
          else
            write(50,220,err=902) M
          end if
          if(.not.lEquil) then
            if(lMobIm(M)) then
              write(50,222,err=902)
            else
              write(50,224,err=902)
            end if
          end if
          if(abs(ChPar(jjj+8,M)-0.0).gt.1.e-12.or.
     !       abs(ChPar(jjj+9,M)-1.0).gt.0.001) lLinear(jj)=.false.
	    if(lBact.and.(ChPar(jjj+18,M).gt.0..or.ChPar(jjj+15,M).gt.0.))
     !                                         lLinear(jj)=.false.
12      continue
13    continue
      do 14 jj=1,NS*16+4
        TDep(jj)=0.
        if(jj.le.NS*9) WDep(1,jj)=1.
        if(jj.le.NS*9) WDep(2,jj)=0.
14    continue
      do 16 jj=1,NS
        do 15 i=1,11
          CumCh(i,jj)=0.
15      continue
        if(lTDep) then
          jjj=(jj-1)*16
          if(jj.eq.1) read(30,*,err=901)
          read(30,*,err=901)
          read(30,*,err=901) (TDep(jjj+j),j=5,6)
          read(30,*,err=901)
          read(30,*,err=901) (TDep(jjj+j),j=7,20)
        end if
16    continue
      if(iMoistDep.eq.1) then
        read(30,*,err=901)
        read(30,*,err=901)
        read(30,*,err=901) nPar2
        do 19 jj=1,NS
          jjj=(jj-1)*9
          read(30,*,err=901)
          read(30,*,err=901) (WDep(1,jjj+j),j=1,9)  ! Coeff. B
          read(30,*,err=901) (WDep(2,jjj+j),j=1,9)  ! References pressure head h_Ref
          do 18 M=1,NMat
            do 17 j=1,9
              WDep(2+M,jjj+j)=FQ(iModel,WDep(2,jjj+j),Par(1,M)) ! Theta_ref
17          continue
18        continue
19      continue
      end if

      read(30,*,err=901)
      read(30,*,err=901) kTopCh,(cTop(jj),jj=1,NS),kBotCh,
     !                   (cBot(jj),jj=1,NS)
      if(kTopCh.eq.-2) then 
        read(30,*,err=901)
        read(30,*,err=901) dSurf,cAtm
      end if
      write(50,230,err=902) kTopCh,(cTop(jj),jj=1,NS)
      write(50,240,err=902) kBotCh,(cBot(jj),jj=1,NS)
      read(30,*,err=901)
      read(30,*,err=901) tPulse
      write(50,250,err=902) tPulse
      return

*     Error when reading from an input file
901   ierr=1
      return
*     Error when writing into an output file
902   ierr=2
      return
*     Bulk Density is equal to zero 
903   ierr=3
      return

110   format(//' Solute transport information'/1X,28('='))
120   format(/' Upstream weighting finite-element method')
130   format(/' Galerkin finite-element method')
140   format (/' Artificial dispersion is added when Peclet number is',
     !         ' higher than',f10.3)
150   format(//' lTDep     lWDep     cTolA     cTolR   MaxItC'
     !        /l3,6x,i3,e13.3,f10.4,i7/
     !        //' Mat.     Bulk.D.    DispL    Fraction  Immobile WC')
160   format(i3,f13.4,3f10.4)
170   format(/'    Dif.w.      Dif.g.   ',50('-'),' (',i2,'.solute)')
180   format(2e12.4/' Mat.     KS         Nu         Beta      Henry
     !  SinkL1     SinkS1     SinkG1     SinkL1`    SinkS1`    SinkG1`
     !  SinkL0     SinkS0     SinkG0      Alfa')
190   format(i4,14e11.4)
200   format(/' Langmuir nonlinear adsorption isotherm for material ',
     !       i2)
210   format(/' Freundlich nonlinear adsorption isotherm for material ',
     !       i2)
220   format(/' No adsorption or linear adsorp. isotherm for material ',
     !       i2)
222   format(/' Physical non-equilibrium solute transport with mobile an
     !d imobile water.')
224   format(/' Chemical non-equilibrium solute transport with kinetic a
     !nd equilibrium sorption sites.')
230   format(/' kTopCh      cTop(1...NS)'/i4,7x,20e10.3)
240   format(/' kBotCh      cBot(1...NS)'/i4,7x,20e10.3)
250   format(/' tPulse =   ',f15.3)
      end

************************************************************************

      subroutine OpenSoluteFiles(NS,cDataPath,iLengthPath,cFileName,
     !                           ierr)

      character cFileName*200,cDataPath*200,cName*12,ch1*1,cName1*13,
     !          ch2*1
      
      do 11 i=1,NS
        if(i.le.9) then
          write(ch1,'(i1)') i
          cName = '/SOLUTEX.OUT'
          cName(8:8) = ch1
          cFileName = cDataPath(1:iLengthPath)//cName
        else
          write(ch1,'(i1)') 1
          write(ch2,'(i1)') i-10
          cName1 = '/SOLUTEXX.OUT'
          cName1(8:8) = ch1
          cName1(9:9) = ch2
          cFileName = cDataPath(1:iLengthPath)//cName1
        end if
        open(80+i,file=cFileName, status='unknown',err=901)
11    continue
      return

*     Error when reading from an input file 
901   ierr=1
      return
      end

*#######################################################################
*
*     iGetFileVersion - vrati verzi souboru
*
*#######################################################################
      integer*4 function iGetFileVersion(FileUnit,iText)
      integer FileUnit, iText
      character cLine*255, cVersion*10, cPcpStr*17
      integer*2 iVersion

      call EmptyStr(cVersion)   
      call EmptyStr(cLine)   

      iGetFileVersion = 0
      iVersion = 0
      rewind(FileUnit,err=1000)

      if(iText.eq.1) then
        ! Textovy soubor
        read(FileUnit, '(a)', err=1000) cLine
        i = index(cLine,'Pcp_File_Version=')
        if(i.eq.1) then
          iPcpLen = 17
          do i=iPcpLen+1,LEN_TRIM(cLine)
            j = i-iPcpLen
            cVersion(j:j) = cLine(i:i)
          end do
          read(cVersion,*,err=1000) iVersion
          iGetFileVersion = iVersion
        else
          rewind(FileUnit,err=1000)
        end if
      else
        ! Binarni soubor
        read(FileUnit,end=1000,err=1000) cPcpStr
        i = index(cPcpStr,'Pcp_File_Version=')
        if(i.eq.1) then
          read(FileUnit,end = 1000,err=1000) iVersion
          iGetFileVersion = iVersion
        else
          rewind(FileUnit,err=1000)
        end if
      end if   

1000  return
      end

*#######################################################################
*
*     EMPTYSTR - vycisteni stringu            
*
*#######################################################################
      subroutine EmptyStr(cString)
      character*(*) cString

      iMax = Len(cString)
      do i=1,iMax
        write(cString(i:i),100,err=200) ' '
      end do   
100   format(a1)
      return

200   write(*,'(1x,a)')  ' Internal error ! '
      end

************************************************************************

      subroutine Init(CosAlf,NTab,ItCum,TLevel,ALevel,PLevel,hRoot,
     !                vRoot,IterW,IterC,dtMaxC,wCumT,wCumA,err,lVarBC,
     !                NSD,cRoot,cCumT,cCumA,NumNPD,Sink,wc,CumQ,lMeteo,
     !                lBact,lVapor,lEnBal,lDayVar,lEnter,lFiltr,TauW,
     !                nPrStep,nTabMod,iDualPor,dtMaxT,lExtrap,lPrintD,
     !                lLAI,rExtinct,lDensity,ExcesInt,lMinstep,lPrint,
     !                lSnow,SnowMF,SnowSF,SnowLayer,lSolHyst,cTemp,
     !                iSunSh,iRelHum,xRoot,lCentrif,Radius,lSource,
     !                iNode,hSink,CumSource,GWL0L,tPrintInt,iTort,
     !                iEnhanc,hSeep,lEqInit,lSolub,lSinPrec,lIrrig,
     !                iIrrig,jIrrig,hIrrig,rIrrig,tIrrig,lThAds,
     !                iMoistDep,Sorb3,OmegaC,WTransf,lDualNEq,lFlux,
     !                iCrop,iRootIn,lMassIni,lMetDaily,PLevelI,iNOB,
     !                lSameMat,Sorb2,lActRSU,OmegaS,SPot,lOmegaW,
     !                OmegaW,lSurfact,lEnd,lVaporOut,lFluxOut,lIsotope,
     !                lFumig,lHanks,iCycle,nCycles,tCycle,iSolAU,
     !                lDblVeg,Sink2,lYield,Yield,PTLevel,crWat,crSal,
     !                iInterc,aInterc,lIsSolPres,lPrHeat,CumHeatFlB,
     !                CumHeatFlT,iDrKod,SinkDr,WDrain,lSurf,fET,lCond,
     !                lAPI,lAPI_ET,iInit)
 
      integer TLevel,ALevel,PLevel,PLevelI,PTLevel,err
      logical lVarBC,lMeteo,lBact,lVapor,lEnBal,lDayVar,lEnter,lFiltr,
     !        lVaporOut,lExtrap,lPrintD,lLAI,lDensity,lMinstep,lPrint,
     !        lSnow,lSameMat,lSolHyst,lCentrif,lSource,lEqInit,lSolub,
     !        lSeep,lSinPrec,lIrrig,lThAds,lOpen,lDualNEq,lAPI,lAPI_ET,
     !        lFlux,lSurfact,lEnd,lMassIni,lMetDaily,lActRSU,lOmegaW,
     !        lFluxOut,lIsotope,lFumig,lHanks,lDblVeg,lYield,lIsSolPres,
     !        lPrHeat,lSurf,lCond
      double precision tPrintInt,tCycle
      dimension cRoot(NSD),cCumT(NSD),cCumA(NSD),Sink(NumNPD),
     !          wc(NumNPD),cTemp(NumNPD),Sorb3(NSD,NumNPD),CumQ(17),
     !          NTab(1),tIrrig(4),Sorb2(NSD,NumNPD),Sink2(NumNPD),
     !          SinkDr(NumNPD)
      common /surfactant/ b1,b2,a1,a2,a3,a4,a5,a6,a7,a8

      lPrint  =.true.
      lMinstep=.false.
      lVarBC  =.false. ! Time-variable head/flux BC
      lVaporOut=.true. ! special Nod_inf_v.out file with thermal and isothermal conductivities and fluxes
      lFluxOut=.false. ! special T_Level1.out file with various boundary fluxes
      lCond=.false.    ! Decrease of conductivity with time

      CosAlf=1.
      NTab(1)=100
      ItCum=0
      TLevel=1
      ALevel=1
      PLevel=1
      PTLevel=1
      PLevelI=1
      hRoot=0.
      vRoot=0.
      IterW=0
      IterC=0
      dtMaxC=1.e+30
      dtMaxT=1.e+30
      wCumT=0.
      wCumA=0.
      WTransf=0.
      WDrain=0.
      iDrKod=1 ! =0: drainage flux through the bottom
               ! =1; drainage flux through sides
      err=0
      iNOB=1
      xRoot=0.
      GWL0L=0.
      lEnd=.false.
      do 11 i=1,NSD
        cRoot(i)=0.
        cCumT(i)=0.
        cCumA(i)=0.
11    continue
      do 12 i=1,NumNPD
        Sink(i)=0.
        Sink2(i)=0.
        SinkDr(i)=0.
        wc(i)=0.
        cTemp(i)=0.
        do 12 jS=1,NSD
          Sorb2(jS,i)=0.
          Sorb3(jS,i)=0.
12    continue
      do 13 i=1,16
        CumQ(i)=0.
13    continue

*     New Options - Supported by GUI
*     lBact   - Virus transport, ka,kd concept
*     lFiltr  - Filtration theory
*     lVapor  - Vapor flow
*     lEnter  - End the run with pushing Enter key
*     nPrStep - Print to the screen and T_Level files at each nPrStep
*     lPrintD - Print at a daily (given) interval
*     tPrintInt- Print interval
*     nTabMod - Kode for the input of the soil hydraulic properties tables
*     iDualPor- Dual porosity model; 
*               = 1: transfer proportional to difference in water contents
*               = 2: transfer proportional to difference in pressure heads
*     lSameMat- The same solute or heat transport parameters for all layers
*     lDualNEq- Both physical and chemical nonequilibrium are considered simultaneously, two-site sorption in the mobile zone
*               (fraction of equilibrium sites in mobile zone - ChPar(13), Rate of kinetic sorption - ChPar(16)
*     lMeteo  - Meterorological input to calculate ET
*     lLAI    - Distribution of pET based on LAI
*     lDayVar - Daily variations in root water uptake and evaporation
*     lSinPrec- Sinusoidal distribution of precipitation
*     lEnBal  - Evaporation and heat flux is calculated from energy balance
*     iSunSh  - =0 Sunshine hours; =1 Cloudeness; =2 Transmission coeff.
*     iRelHum - =0 Relative Humidity; =1 Vapor Pressure
*     lMetDaily - Daily variations of meteorological variables are generated from daily average, max, and min data
*     lSnow   - Snow accumulation at the soil surface
*     SnowMF  - Amount of snow melted per 1 degree [cm3/cm2/K/d]
*     SnowSF  - Sublimation constant
*     SnowLayer- Thickness of the snow layer
*     iTort   - tortuosity in solute transport, =0 for Millington and =1 for Moldrup
*     lSeep   - Seepage face initiated by different pressure head
*     hSeep   - seepage face with a different bottom pressure
*     lEqInit - initial noequilibrium phase is in equilibrium with liquid phase
*     lMassIni- initial condition is given in the total concentration [M_solute/M_soil]
*     iMoistDep- reaction rates are dependent on the water content (iMoistDep=2)
*     OmegaC  - Compensated root water uptake
*     lFlux   - Print fluxes in observation nodes instead of temperatures
*     lActRSU - Active root solute uptake
*     OmegaS  - Compensated root solute active uptake
*     SPot    - Potential root solute uptake
*     iSolAU  = Solute for which to consider active solute uptake
*     lOmegaW - Reduction of the potential root solute uptake due to reduction of the root water uptake
*     lIsotope- Solute transport for isotopes, i.e., they evaporate at the surface; change for Christine Stumpp
*     lFumig  - Special Options for fumigants
*     lHanks  - Hanks root water uptake model
*     lDblVeg - Double vegetation, Axx = Beta2, no dual-porosity, no heat transport
*     lYield  - Calculate crop yield
*     Yield   - Yield
*     crWat   - seasonal reduction of root water uptake due to pressure head stress
*     crSal   - seasonal reduction of root water uptake due to salinity stress
*     lIsSolPres - Is solute in the transport domain?
*     lSurf   - divide potential evapotranspiration into potential evaporation
*               and transpiration based on the root growth, i.e., fET=rR/rRMax

      lBact   =.false.
      lFiltr  =.false.
      lVapor  =.false.
      lEnter  =.true.
      lPrintD =.false.
      lSameMat=.false.
      lDualNEq=.false.
      lMeteo  =.false.
      lLAI    =.false.
      lDayVar =.false.
      lSinPrec=.false.
      lEnBal  =.false.
      lMetDaily=.false.
      lSnow   =.false.
      lSeep   =.false.
      lEqInit =.false.
      lMassIni=.false.
      lFlux   =.false.
      lActRSU =.false.
      lOmegaW =.false.
      lIsotope=.false.
      lFumig  =.true.   ! if input Fumigants.in does not exist, lFumig=.false.
      lHanks  =.false.
      nCycles=1
      iCycle=1
      tCycle=0.
      lDblVeg =.false.
      lYield=.false.
      Yield=1.
      crWat=0.
      crSal=0.
      iInterc=0         ! Interception
      aInterc=0.25      ! Interception [mm]
      lIsSolPres=.false.
      lPrHeat=.false.   ! Print heat fluxes
      CumHeatFlB=0.     ! Cumulative bottom heat flux
      CumHeatFlT=0.     ! Cumulative surface heat flux

      nPrStep =1
      tPrintInt=86400.
      nTabMod =10
      iDualPor=0
      iSunSh  =0
      iRelHum =0
      iCrop   =0
      SnowMF  =0.45  ! cm
      SnowSF  =0.4
      iTort   =0
      hSeep   =0.
      iMoistDep =0
      OmegaC  =1.
      OmegaS  =1.
      SPot    =0.
      iSolAU  =1
      OmegaW  =1.
      ExcesInt=0.
      rExtinct=0.39
      SnowLayer=0.
      lSurf   =.false.
      fET     =0.

*     New Options - Not supported by GUI
*     TauW    - Nonequilibrium water flow [Ross and Smettem, 2001]
*     lExtrap - Extrapolation of hNew from previous time step
*     lDensity- Density dependent flow and transport
*     lSolHyst- Hysteresis in solute parameters for mass transfer term for mobile-immobile model (20 in, and 16 out)
*     lCentrif- Gravitational acceleration in the cetrifuge.
*     Radius  - Distance of the sample from the centre of centrifuge
*     lSource - Internal pressure head sink/source [for Bridget]
*     hSink   - Pressure at the internal sink/source
*     iEnhanc - Enhancement factor for vapor flow, =1 OK, =0 no.
*     lSolub  - limit max concentration in the liquid phase (GamG1), enter initial conc. into Sorb,
*             - Solid phase concentration is then transfered into Sorb2 and stays there.
*             - Solubility limit is given in GamG1, Raoult's law is implemented.
*     lThAds  - Adsorption is a function of water content, implemented only with the equilibrium model
*     lAPI    - Modification for American Petroleum Institute
*     lSurfact- effect of surfactant

      lExtrap =.true.
      lDensity=.false.
      lSolHyst=.false.
      lSolub  =.false.
      lThAds  =.false.
      lAPI    =.false.
      lAPI_ET =.false.
      iRootIn =-1
      TauW    =0.
      iEnhanc=1

*     Irrigation Scheduling
      lIrrig=.false.
      iIrrig=1        ! Observation node triggering irrigation
      hIrrig=-200.    ! Pressure head triggering irrigation [L]
      jIrrig=0        ! =0: no irrigation
                      ! =1" irrigation is triggered
                      ! =2" irrigation is active
      rIrrig=2.0      ! Irrigation rate [L/T]
      tIrrig(1)=0.    ! Beginning of irrigation [T]
      tIrrig(2)=0.    ! End of irrigation [T]
      tIrrig(3)=0.3   ! Duration of irrigation [T]
      tIrrig(4)=0.5   ! Offset of irrigation [T]

*     Screen inputs -----------------------------------------------------
      iScreenInput=0

      lSurfact=.false.
      if(iInit.eq.0) then
        a1=1.
        a2=0.
        a3=0.
      	a4=1.
        a5=1.
        a6=0.
        a7=0.
      	a8=1.
        b1=0.
        b2=1.
        if(lSurfact.and.iScreenInput.eq.0) then
          inquire(unit=35,opened=lOpen)
          if(lOpen) then
            read(35,*)
            read(35,*)
            read(35,*) a1,a2,a3,a4
            read(35,*)
            read(35,*) a5,a6,a7,a8
            read(35,*) 
            read(35,*) b1,b2
          end if
        end if
      end if

      lCentrif=.false.
      if(iInit.eq.0) then
        Radius=0.
        if(lCentrif.and.iScreenInput.eq.1) then
          write(*,*) "Enter radius:"
          read(*,*) Radius
        end if
      end if

      lSource=.false.   ! needs lExtrap =.false. and iTemp=0
      CumSource=0.
      if(iInit.eq.0) then
        if(lSource) then
          if(iScreenInput.eq.1) then
            write(*,*) "Enter Node number:"
            read(*,*) iNode
            write(*,*) "Enter Pressure at Sink:"
            read(*,*) hSink
          else if(iScreenInput.eq.0) then
            inquire(unit=35,opened=lOpen)
            if(lOpen) then
              read(35,*)
              read(35,*)
              read(35,*) iNode,hSink
              write(*,*) "iNode=",iNode
              write(*,*) "hSink=",hSink
            end if
          end if
        end if
      end if
*     Screen inputs -----------------------------------------------------

      if(iScreenInput.eq.0) then ! Input file "Options.In"
        inquire(unit=35,opened=lOpen)
        if(lOpen) then
c          read(35,*)
c          read(35,*)
          if(iInit.eq.0) write(*,140) 
140       format(/' Information from the Options.in input file')  
        end if
      end if

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
