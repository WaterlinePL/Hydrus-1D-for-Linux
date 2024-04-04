* Source file SINK.FOR |||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine SetSnk(N,NMat,MatNum,x,hRoot,vRoot,Sink,TPot,hNew,
     !                  lMoSink,lSolRed,lSolAdd,P0,POptm,P2H,P2L,P3,r2H,
     !                  r2L,aOsm,c50,P3c,Beta,lChem,NS,NSD,Conc,cRoot,
     !                  lMsSink,ThNew,ParD,dt,OmegaC,iModel,Con,lOmegaW,
     !                  OmegaW,lHanks,t,lDblVeg,Sink2,Beta2,TPot2,crWat,
     !                  crSal)

      double precision t
      logical lChem,lMoSink,lSolRed,lSolAdd,lMsSink,lHanks,lOmegaW,
     !        lDblVeg
      dimension x(N),MatNum(N),hNew(N),POptm(NMat),Beta(N),Sink(N),
     !          Conc(NSD,N),cRoot(NS),aOsm(NS),ThNew(N),ParD(11,NMat),
     !          Con(N),Beta2(N),Sink2(N)
      common /DoubleVeg/ xP0,xP2H,xP2L,xP3,xr2H,xr2L,xPOptm

      Compen=1.
      nStep=1
      if(OmegaC.lt.1.) nStep=2
      Omega=0.
      Omega2=0.

      vRoot=0.
      hRoot=0.
      ARoot=0.
      do 11 ii=1,NS
        cRoot(ii)=0.
11    continue

      if(lHanks) then
        hR=P3
        do 22 iter=1,5
          Sum1=0.
          Sum2=0.
          xConst=1.  ! Penalty for distance
          do 21 i=2,N
            if(Beta(i).gt.0.) then
              if(i.eq.N) then
                dxM=(x(i)-x(i-1))/2.
              else
                dxM=(x(i+1)-x(i-1))/2.
              end if
              if(hNew(i).gt.hR-xConst*x(i)) then
                Sum1=Sum1+Con(i)*Beta(i)*(hNew(i)+xConst*x(i))*dxM
                Sum2=Sum2+Con(i)*Beta(i)*dxM
              end if
            end if
21        continue
          hR=max((Sum1-TPot)/Sum2,P3)
22      continue
        do 23 i=1,N
          if(Beta(i).gt.0.) then
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            Sink(i)=max(-Con(i)*Beta(i)*(hR-xConst*x(i)-hNew(i)),0.)
c            Sink(i)=-Con(i)*Beta(i)*(hR-xConst*x(i)-hNew(i))
            vRoot=vRoot+Sink(i)*dxM
            hRoot=hRoot+hNew(i)*dxM
            ARoot=ARoot+dxM
          end if
23      continue
        if(ARoot.gt.0.001) hRoot=hRoot/ARoot
        write(62,110) t,hR
110     format(e14.7,f12.4)
        return
      end if

      do 16 iStep=1,nStep
        rWat=0. ! profile reduction of root water uptake due to pressure head stress
        rSal=0. ! profile reduction of root water uptake due to salinity stress
        do 13 i=2,N
          if(Beta(i).gt.0.) then
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            M=MatNum(i)
            hRed=hNew(i)
            SAlfa=1.
            if(lChem.and.lSolRed) then
              cRed=0.
              do 15 j=1,NS
                cRed=cRed+aOsm(j)*Conc(j,i)
15            continue
              if(lSolAdd) then
                hRed=hRed+cRed
              else
                SAlfa=FSAlfa(lMsSink,cRed,c50,P3c)
              end if
            end if
            Alfa=FAlfa(lMoSink,TPot,hRed,P0,POptm(M),P2H,P2L,P3,r2H,r2L)
            if(iStep.ne.nStep) then
              Omega=Omega+Alfa*SAlfa*Beta(i)*dxM
              if(lDblVeg) goto 17
              goto 13
            else
              Compen=1.
              if(Omega.lt.OmegaC.and.Omega.gt.0.) Compen=OmegaC
              if(Omega.ge.OmegaC)                 Compen=Omega
            end if
            if(Compen.gt.0.) Sink(i)=Alfa*SAlfa*Beta(i)*TPot/Compen
            if(ThNew(i)-0.00025.lt.ParD(1,MatNum(i))) Sink(i)=0.
            if(lMoSink) PMin=P3
            if(.not.lMosink) PMin=10.*P0
            ThLimit=FQ(iModel,PMin,ParD(1,MatNum(i)))
c            Sink(i)=min(Sink(i),0.5*(ThNew(i)-ParD(1,MatNum(i)))/dt)
            Sink(i)=min(Sink(i),max(0.,0.5*(ThNew(i)-ThLimit)/dt))
            vRoot=vRoot+Sink(i)*dxM
            hRoot=hRoot+hNew(i)*dxM
            do 12 ii=1,NS
              if(lChem) cRoot(ii)=cRoot(ii)+Conc(ii,i)*dxM
12          continue
            ARoot=ARoot+dxM
            rWat=rWat+ Alfa*Beta(i)*dxM
            rSal=rSal+SAlfa*Beta(i)*dxM
          else
            Sink(i)=0.
          end if
17        if(lDblVeg) then
            if(Beta2(i).gt.0.) then
              if(i.eq.N) then
                dxM=(x(i)-x(i-1))/2.
              else
                dxM=(x(i+1)-x(i-1))/2.
              end if
              M=MatNum(i)
              hRed=hNew(i)
              SAlfa=1.
              if(lChem.and.lSolRed) then
                cRed=0.
                do 18 j=1,NS
                  cRed=cRed+aOsm(j)*Conc(j,i)
18              continue
                if(lSolAdd) then
                  hRed=hRed+cRed
                else
                  SAlfa=FSAlfa(lMsSink,cRed,c50,P3c)
                end if
              end if
              Alfa=FAlfa(lMoSink,TPot2,hRed,xP0,xPOptm,xP2H,xP2L,xP3,
     !                   xr2H,xr2L)
              if(iStep.ne.nStep) then
                Omega2=Omega2+Alfa*SAlfa*Beta2(i)*dxM
                goto 13
              else
                Compen=1.
                if(Omega2.lt.OmegaC.and.Omega2.gt.0.) Compen=OmegaC
                if(Omega2.ge.OmegaC)                  Compen=Omega2
              end if
              if(Compen.gt.0) Sink2(i)=Alfa*SAlfa*Beta2(i)*TPot2/Compen
              if(ThNew(i)-0.00025.lt.ParD(1,MatNum(i))) Sink2(i)=0.
              if(lMoSink) PMin=P3
              if(.not.lMosink) PMin=10.*P0
              ThLimit=FQ(iModel,PMin,ParD(1,MatNum(i)))
              Sink2(i)=min(Sink2(i),max(0.,0.5*(ThNew(i)-ThLimit)/dt))
              vRoot=vRoot+Sink2(i)*dxM
            else
              Sink2(i)=0.
            end if
          else
            Sink2(i)=0.
          end if
c          if(Beta(i).lt.0.) then ! Eddy Woehling's modification
c            if(i.eq.N) then
c              dxM=(x(i)-x(i-1))/2.
c            else
c              dxM=(x(i+1)-x(i-1))/2.
c            end if
c            Sink(i)=Beta(i)*rBot
c            Sink(i)=max(Sink(i),0.5*(ThNew(i)-ParD(2,MatNum(i)))/dt)
c          end if
13      continue
16    continue
      if(ARoot.gt.0.001) then
        hRoot=hRoot/ARoot
        do 14 ii=1,NS
          cRoot(ii)=cRoot(ii)/ARoot
14      continue
      end if
      if(lOmegaW.and.TPot.gt.0.) OmegaW=vRoot/TPot
      crWat=crWat+rWat*TPot*dt
      crSal=crSal+rSal*TPot*dt
      return
      end

************************************************************************

*     Subroutine calculating root solute uptake with and without compensation

      subroutine SetSSnk(jS,N,t,x,Beta,Sink,SinkS,NSD,Conc,OmegaW,
     !                   cRootMax,lActRSU,OmegaS,SPot,iSolAU,rKM,cMin,
     !                   TLevel,rA,rB,rC,rD)

      double precision t
      integer TLevel
      logical lActRSU
      dimension x(N),Beta(N),Sink(N),SinkS(N),Conc(NSD,N)

*     Inputs:
*     SPot      - potential root solute uptake
*     OmegaS    - solute stress index
*     rKM       - Michaelis-Menten constant
*     lActRSU   - consider active root solute uptake
*     iSolAU    - Solute for which to consider active solute uptake
*     cRootMax  - maximum concentration for the passive solute uptake

*     From Water Flow
*     Sink(i)   - Root water uptake
*     OmegaW    - ratio of actual and potential transpiration

*     SPUptake  - passive root solute uptake (step 1)
*     SAUptakeP - potential active solute uptake (step 1)
*     SAUptakeA - uncompensated actual active solute uptake (step 2)
*     SAUptakeAN- uncompensated actual active solute uptake (step 2)
*     SAUptakeA - compensated actual active solute uptake (step 3)
*     SinkS(i)  - local passive (and active; if considered) solute uptake

c      if(TLevel.eq.1.and.lActRSU.and.jS.eq.iSolAU) write(78,110)

*     Initialization
      Compen=1.
      nStep=1
      if(lActRSU)                  nStep=2
      if(lActRSU.and.OmegaS.lt.1.) nStep=3
*     step 1: Passive uptake
*     step 2: Active uptake without compensation
*     step 3: Active uptake with compensation
      if(jS.ne.iSolAU) nStep=1    ! Solute for which the active uptake is considered
      Omega=0.
      SPUptake=0.
      do 10 i=1,N
        SinkS(i)=0.
10    continue

      do 12 iStep=1,nStep
        SAUptakeA=0.
        do 11 i=1,N
          if(Beta(i).gt.0.) then
            UR=1.                                                  !Naftali
            if(jS.eq.20) UR=(1.-rC*Conc(jS,i))*(1.-rD*Conc(jS-1,i)) !Naftali (jS.eq.2)
            rKMx=rKM                                               !Naftali
            if(jS.eq.20) rKMx=rKM*(rA+rB*Conc(jS-1,i))             !Naftali (jS.eq.2)
            if(i.eq.N) then
              dxM=(x(i)-x(i-1))/2.
            else if(i.eq.1) then
              dxM=(x(i)-x(i+1))/2.
            else
              dxM=(x(i+1)-x(i-1))/2.
            end if
            cc=amax1(Conc(jS,i)-cMin,0.)
            if(iStep.eq.1) then
              SinkS(i)=Sink(i)*amax1(amin1(Conc(jS,i),cRootMax),0.)
              SPUptake=SPUptake+SinkS(i)*dxM
*             This is needed only for the last node, but that node may not have beta
              SAUptakeP=amax1(SPot*OmegaW-SPUptake,0.)
            else if(iStep.eq.2) then
              AUptakeA=0.
              if(rKMx+cc.gt.0.) 
     !          AUptakeA=cc/(rKMx+cc)*Beta(i)*SAUptakeP*UR
              Omega=Omega+AUptakeA*dxM
              if(nStep.eq.2) SinkS(i)=SinkS(i)+AUptakeA
*             This is needed only for the last node, but that node may not have beta
              SAUptakeA =Omega
              SAUptakeAN=Omega
              if(SAUptakeP.gt.0.) Omega1=Omega/SAUptakeP
            else if(iStep.eq.3) then
*             This is needed only for the first node, but that node may not have beta
              if(Omega1.lt.OmegaS.and.Omega1.gt.0.) Compen=OmegaS
              if(Omega1.ge.OmegaS)                  Compen=Omega1
              if(Compen.gt.0..and.rKMx+cc.gt.0.) 
     !          AUptakeA=cc/(rKMx+cc)*Beta(i)*SAUptakeP/Compen*UR
              SinkS(i)=SinkS(i)+AUptakeA
              SAUptakeA=SAUptakeA+AUptakeA*dxM
            end if
          else
            SinkS(i)=0.
          end if
11      continue
c        if(lActRSU.and.jS.eq.iSolAU.and.iStep.eq.nStep) 
c     !    write(78,120) t,SPUptake,SAUptakeP,SAUptakeA,SAUptakeAN ! the last is uncompensated
        tt=sngl(t)
12    continue
      return

110   format(/' SPUptake  - passive root solute uptake'/
     !        ' SAUptakeP - potential active solute uptake'/
     !        ' SAUptakeA - compensated actual active solute uptake'/
     !        ' SAUptakeAN- uncompensated actual active solute uptake'//
     !'         Time         SPUptake   SAUptakeP   SAUptakeA  SAUptakeA
     !N'/'          [T]         [M/L2/T]    [M/L2/T]    [M/L2/T]    [M/L
     !2/T]'/)
120   format(3x,e14.7,1x,4e12.4)
      end

************************************************************************

*     Stress response function for osmotic stress (for multiplicative model)

      real function FSAlfa(lMode,cRed,c50,P3c)

      logical lMode

      if(lMode) then
        FSAlfa=0.
        if(abs(c50).gt.0) FSAlfa=1./(1.+(cRed/c50)**P3c)
      else
        if(cRed.le.c50) then
          FSAlfa=1.
        else
          FSAlfa=max(0.,1.-(cRed-c50)*P3c*0.01)
        end if 
      end if
      return
      end

************************************************************************

*     Stress response function for water stress (and for additive model)

      real function FAlfa(lMoSink,TPot,h,P0,P1,P2H,P2L,P3,r2H,r2L)

      logical lMoSink

      if(lMoSink) then
        if(TPot.lt.r2L) P2=P2L
        if(TPot.gt.r2H) P2=P2H
        if((TPot.ge.r2L).and.(TPot.le.r2H).and.abs(r2H-r2L).gt.0.)
     !    P2=P2H+(r2H-TPot)/(r2H-r2L)*(P2L-P2H)
        FAlfa=0.0
        if((h.gt.P3).and.(h.lt.P2).and.abs(P2-P3).gt.0.) 
     !                              FAlfa=(h-P3)/(P2-P3)
        if((h.ge.P2).and.(h.le.P1)) FAlfa=1.0
        if((h.gt.P1).and.(h.lt.P0).and.abs(P1-P0).gt.0) 
     !                              FAlfa=(h-P0)/(P1-P0)
        if((P1.eq.0.).and.(P0.eq.0.)) FAlfa=1.0
*       Uptake even at full saturation, when both P1 and P2 are equal to zero
        if(h.ge.P2.and.P1.eq.0..and.P0.eq.0.) FAlfa=1.0
      else
        FAlfa=1./(1.+(h/P0)**P3)
      end if
      return
      end

************************************************************************

*     Stress response function for additive model (based on John Letey)

      real function FAlfa1(h,h50,c,c50,P3)

      FAlfa1=1./(1.+(h/h50+c/c50)**P3)
      return
      end

************************************************************************

      subroutine SetRG(NumNP,x,Beta,t,tRMin,tRHarv,xRMin,xRMax,RGR,
     !                 xRoot,lRoot,iRootIn,nGrowth,rGrowth,tRPeriod,fET)

      dimension x(NumNP),Beta(NumNP),rGrowth(1000,5)
      double precision t
      logical lRoot

      if(lRoot.and.iRootIn.eq.1) then
        i=1000
        j=5
        call Table(nGrowth,rGrowth,i,j,t,rDummy,rDummy,rDummy,xRoot)
      end if

      if(lRoot.and.iRootIn.eq.2) then
        tRoot=amod(sngl(t),tRPeriod)
        if(tRoot.lt.tRMin.or.tRoot.gt.tRHarv) then
          fET=0.
          do 11 i=1,NumNP
            Beta(i)=0.
11        continue
          return
        end if
        xR=xRMax
        if(xRMin.le.0.001) xRMin=0.001
        tt=tRoot-tRMin
        xR=(xRMax*xRMin)/(xRMin+(xRMax-xRMin)*exp(-RGR*tt))
        fET=amin1(amax1(0.,xR/xRMax),1.)
      else
        xR=xRoot
        fET=1.
      end if

      SBeta=0.
      Beta(NumNP)=0.
      Beta(1)=0.
      do 12 i=2,NumNP-1
        if(i.eq.NumNP-1.and.x(NumNP)-x(NumNP-1).gt.xR) then
          Beta(i)=1.        
        else if(x(i).lt.x(NumNP)-xR) then
          Beta(i)=0.
        else if(x(i).lt.x(NumNP)-0.2*xR) then
          Beta(i)=2.08333/xR*(1-(x(NumNP)-x(i))/xR)
        else
          Beta(i)=1.66667/xR
        end if
        if(i.ne.NumNP) then
          SBeta=SBeta+Beta(i)*(x(i+1)-x(i-1))/2.
        else
          SBeta=SBeta+Beta(i)*(x(i)-x(i-1))/2.
        end if
12    continue

      if(xR.gt.0.) then
        if(SBeta.lt.0.0001) then
          Beta(NumNP-1)=1./((x(NumNP)-x(NumNP-2))/2.)
        else
          do 13 i=2,NumNP-1
            Beta(i)=Beta(i)/SBeta
13        continue
        end if
      end if

      return
      end

************************************************************************

      subroutine Yield1(lChem,Yield,rRoot,vRoot,rSoil,cvChR,NS,SPot,CTp,
     !                  CRp,dt,t,rRootU,SPotU,TLevel,iKod)

      logical lChem,lAdjust
      dimension cvChR(NS)
      double precision t
      integer TLevel
  
*     Yield   - crop yield
*     rYield  - rate of yield loss
*     rYieldW - rate of yield loss due to water and salinity stress
*     rYieldS - rate of yield loss due to (lack of) nutrient stress (must be NS solute)
*     lAdjust - adjust potential fluxes due to past stress reduction
*     rRootU  - potential transpiration unadjusted to prior lost of yield
*     SPotU   - potential nutrient depand unadjusted to prior lost of yield
      lAdjust=.false.
      if(iKod.eq.0) then
        if(TLevel.eq.1) write(78,110) 
        rYieldW=0.
        if(CTp.gt.0.) rYieldW=(rRootU-vRoot)/CTp
        rYield=rYieldW
        if(lChem) then
          rYieldS=0.
          if(CRp.gt.0.) rYieldS=amax1(0.,(SPot-cvChR(NS))/CRp)
          rYield=max(rYieldW,rYieldS)
c          if(rYieldW.gt.0..and.rYieldS.gt.0.) rYield=(rYieldW+rYieldS)/2.
        end if
        Yield=max(0.,Yield-rYield*dt)
      else
        rRootU=rRoot   ! unadjusted
        SPotU=SPot
        if(lAdjust) then
          rRoot=rRoot*Yield
          rSoil=rSoil+(rRootU-rRoot)
          SPot=SPot*Yield
        end if
      end if

      if(iKod.eq.0) 
     !  write(78,120) t,rYieldW,rYieldS,rRoot,vRoot,SPot,cvChR(NS),Yield
      tt=sngl(t)

      return
110   format(/' rYieldW - rate of yield loss due to water stress'/
     !        ' rYieldS - rate of yield loss due to nutrient stress'/
     !        ' rRoot   - potential transpiration'/
     !        ' vRoot   - actual transpiration'/
     !        ' SPot    - potential nutrient demand'/
     !        ' cvChr   - actual nutrient root uptake'/
     !        ' Yield   - actual yield'//
     !'         Time         rYieldW     rYieldS      rRoot       vRoot 
     !      SPot       cvChr       Yield'/
     !'          [T]           [-]         [-]        [L/T]      [L/T]  
     !   [M/L2/T]    [M/L2/T]       [-]'/)
120   format(3x,e14.7,1x,10e12.4)
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

*     Jasper Vrugt function

      subroutine RootInN(NumNP,Beta,z)

      dimension Beta(NumNP),z(NumNP)

*     read input
      read(44,*)
      read(44,*) Zm,Z0,Za

*     coordinate of the surface
      ZMax =z(NumNP)

*     calculate non-normalized uptake intensity
      r1=0.
      r2=0.
      do 11 i=1,NumNP
        if(abs(Zm).gt.1.e-5) then
          r1=(Zm-(ZMax-z(i)))/Zm
          rRootA=Za
          if((ZMax-z(i)).gt.Z0) rRootA=1.
          r2=rRootA/(Zm)*abs(Z0-(ZMax-z(i)))
          r2=exp(-r2)
        end if
        Beta(i)=amax1(r1*r2,0.)
11    continue

*     normalize uptake intensity
      SBeta=Beta(NumNP)*(z(NumNP)-z(NumNP-1))/2.
      do 12 i=2,NumNP-1
        SBeta=SBeta+Beta(i)*(z(i+1)-z(i-1))/2.
12    continue
      do 13 i=2,NumNP
        if(SBeta.gt.0.) then
          Beta(i)=Beta(i)/SBeta
        else
          Beta(i)=0.
        end if
13    continue

      return
      end

************************************************************************

*     Reads from Options.in the parameters for the Couvreur et al. (2012) root water uptake model
*     and stress function
*     Initializes other variables needed by this model

      subroutine CouIn(StresFun,hx_min,Krs,Kcomp,NumNP,x,zSurf,BCtp,BCr,
     !                 lCouvreur,lDayVar)

      double precision Krs,Kcomp
      real hx_min
      integer StresFun,BCtp
      logical lCouvreur,lDayVar
      dimension x(NumNP)

      read(35,*,err=901)
      read(35,*,err=901)
      read(35,*,err=901)
      read(35,*,err=901) StresFun,BCtp,lDayVar
*               StresFun = 0: the transpiration rate is not sensitive to 
*                             the plant collar water potential (Tact = Tpot).
*               StresFun = 1: the transpiration rate can be adjusted depending on 
*                             the plant collar water potential.
*               BCtp - type of boundary condition at the plant collar 
*               = 1: for a pressure head BC (Collar water potential directly defined 
*                                            as root system boundary condition); 
*               = 2: for a flux BC.
*               lDayVar - Daily variations in root water uptake and evaporation
      read(35,*,err=901)
      read(35,*,err=901)
      read(35,*,err=901) hx_min
*               hx_min [L] is the threshold collar water potential that triggers 
*               a reduction of Tact, as compared to Tpot.
      read(35,*,err=901)
      read(35,*,err=901)
      read(35,*,err=901) Krs   
*               Equivalent conductance of the whole root system [cm³/hPa/day]
      read(35,*,err=901)
      read(35,*,err=901)
      read(35,*,err=901) Kcomp 
*               Compensatory RWU conductance of the root system [cm³/hPa/day]

      if(BCtp.eq.1) BCr=hx_min

      zSurf=x(NumNP)

      return

901   lCouvreur=.false.
      return
      end

************************************************************************

      subroutine SetSnkV(NumNP,SSF,zSurf,x,StresFun,Krs,Kcomp,hx_min,
     !                   BCr,BCtp,hNew,hTot,Sink,rRoot,hRoot,vRoot)

*     BCr - value of the head/flux at the plant collar.
*     VNodes - a vector containing the volume associated to each soil node.
*     zSurf - coordinate of the surface node [L]

      integer StresFun,BCtp
      double precision Krs,Kcomp
      real hx_min
      dimension hNew(NumNP),x(NumNP),hTot(NumNP),SSF(NumNP),Sink(NumNP)

      vRoot=0.
      hRoot=0.
      ARoot=0.

*     BCtp - type of boundary condition at the plant collar 
*          = 2: for a flux BC.
      if(BCtp.eq.2) BCr=rRoot    ! Potential transpiration

      hSeq=0.
      do 11 n=1,NumNP
        if(SSF(n).gt.0.) then
          if(n.eq.1) then
            dx=(x(2)-x(1))/2.
          else if(N.eq.NumNP) then
            dx=(x(NumNP)-x(NumNP-1))/2.
          else
            dx=(x(n+1)-x(n-1))/2.
          end if
          hTot(n)=hNew(n)+(x(n)-zSurf) ! The total soil water potential is the sum of 
*                                       the matric potential and the gravitational potential.
          hSeq=hSeq+SSF(n)*hTot(n)*dx  ! The "DOT_PRODUCT" of hTot by SSF.
        end if
11    continue

      if(BCtp.eq.2) then          ! Flux BC for root uptake
        Tpot=BCr                  ! [L/T]
        PHcollar=-abs(Tpot)/Krs+hSeq  ! First calculate the collar water potential considering no stress.
        if(StresFun.eq.1.and.PHcollar.lt.hx_min) then			!We check if we are under the threshold collar water potential.
c          write(*,910) PHcollar,hx_min
          PHcollar=hx_min         ! In case of stress, the collar water potential stays blocked at the treshold value, and will be used to calculate the actual transpiration rate.
          Tact=Krs*(hSeq-PHcollar)![L/T] Convert the actual collar water potential into the actual transpiration (Couvreur jul 2011)
        else
          Tact=Tpot               ! [L/T]
        end if
      else if(BCtp.eq.1) then     ! Head BC for root uptake  
        PHCollar=BCr
        Tact=Krs*(hSeq-PHcollar)  ! [L/T] Convert the actual collar water potential into the actual transpiration (Couvreur jul 2011)
      end if

      do 12 n=1,NumNP
        if(SSF(n).gt.0.) then
c          Sink(n)=abs(Tact)*SSF(n)/VNode(n)+
c     !            Kcomp*(hTot(n)-hSeq)*SSF(n)/VNode(n)
          Sink(n)=abs(Tact)*SSF(n)
     !            +Kcomp*(hTot(n)-hSeq)*SSF(n)
          vRoot=vRoot+Sink(n)*dx
          hRoot=hRoot+hNew(n)*dx
          ARoot=ARoot+dx
        end if
12    continue
      if(ARoot.gt.0.001) hRoot=hRoot/ARoot
      return

910   format('Stress in the collar xylem: PHcollar=',f12.3,' is lower th
     !an criterion ',f12.3)
      end

************************************************************************