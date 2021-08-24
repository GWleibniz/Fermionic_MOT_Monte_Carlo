function [deltaVLasers,atomicstate,ZeemanSublevel]   =  getVelocityChangev2(LaserInfoMat,AtomInfoMat,cteArray,CGTable,PolinomialFitMOTfield)%Initialize variabl

    UpdateStepsize = cteArray(7);
    vrecoil =cteArray(6);
    Pdecay = cteArray(8);
    
    LaserDirMat = LaserInfoMat(:,3:5);
    TransitionAdressedMat = LaserInfoMat(:,26);
    
    deltaVLasers = (zeros(1,3));
    zerosv = (zeros(1,3));
    auxVec = [-1,1,0];
    deltaJtransition = [1,-1,0];
    
    %Atom info
    atomicstate =  AtomInfoMat(:,7);
    ZeemanSublevel =  AtomInfoMat(:,8);
    %if we are in the ground state
    if (atomicstate == 0)

        [Psp,Psm,Spp] =  ScatteringRatev3(CGTable,cteArray,LaserInfoMat,AtomInfoMat,PolinomialFitMOTfield);
        ScatteringRateMatrix = [Psp,Psm,Spp];
%         nlasers = size(LaserInfoMat,1);
%         
%         [Psp,Psm,Pp] = arrayfun(@(n)ScatteringRatev3(CGTable,cteArray,LaserInfoMat(n,:),AtomInfoMat),1:nlasers);
        
%         ScatteringRateMatrix = [Psp',Psm',Pp'];
        ProbabilityVectorScattering = UpdateStepsize*reshape(ScatteringRateMatrix',[],prod(size(ScatteringRateMatrix)));

        %Use a random number to find wich laser will scatter
        ProbabilityVector = cumsum(ProbabilityVectorScattering/sum(ProbabilityVectorScattering));
        probabilityNumber = rand();  
        [~,selected] = max(probabilityNumber<ProbabilityVector);

        %Selected Laser
        Direction = LaserDirMat(ceil(selected/3),:);
        AdressedTransition = TransitionAdressedMat (ceil(selected/3));
        %Probability to scatter from the selected laser
        Psc = ProbabilityVectorScattering(selected);

        %Transition type
%         TransitionType2 = rem(selected+2,3)+1;% 1 = sp, 2 = sm , 0 = pi;
        TransitionType = selected + 3 - floor((selected+2)./3)*3;
%         gather(ScatteringRateMatrix)

        if (rand() < Psc ) %Atom scatter a photon and goes to the excited AtomicState
            if (AdressedTransition)
                atomicstate = 1;
            else
                atomicstate = 2;
            end

          ZeemanSublevel =ZeemanSublevel + deltaJtransition(TransitionType);
           deltaVLasers = vrecoil*Direction;
        else
           deltaVLasers = zerosv;
        end



    %if we are on one of the excited states
    elseif ( atomicstate ~=0 )
        if (rand() < Pdecay ) %Atom  emits a photon and goes to the ground AtomicStat

            %Calculate relative transition strenght of the possible,sp sm and pi decay paths
            if (atomicstate ==1)
                CGCoeffDecayTable = CGTable(21:32,:) ;
                n = (ZeemanSublevel*2-1)/2+7;
            elseif ( atomicstate == 2)
                CGCoeffDecayTable =   CGTable(33:end,:) ;
                n = (ZeemanSublevel*2-1)/2+6;
            end

            vectorDecay = CGCoeffDecayTable(n,:).^2;


            %Select decay path using a randon number ( same as above)
            ProbabilityVector = cumsum(vectorDecay/sum(vectorDecay));
            probabilityNumber = rand();
            [~,selected] = max(probabilityNumber<ProbabilityVector);


            %Calculate change of Zeeman sublevel based on decay path

           ZeemanSublevel =ZeemanSublevel + auxVec(selected);

            %Bring the atom back to the ground state
            atomicstate = 0;

            %Momentum kick in random direction
            v = randn(1,3);
            deltaVLasers = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)))*vrecoil;

        else
            deltaVLasers = zerosv;

        end

    end
end

  function [PSctsp,PSctsm,PSctpi] = ScatteringRatev3(CGTable,cteArray,LaserInfoMat,AtomInfoMat,PolinomialFitMOTfield)
    %Constants
    G =cteArray(9:11);
    kred =  cteArray(4);      
    CGtableTrapping =   CGTable(1:10,:);
    CGtableSitrring =  CGTable(11:20,:);
    Isat = cteArray(3);
    GammaRed =  cteArray(1);
    
    %Laser Info
    NumberOfcomblines =  LaserInfoMat(:,1);
    magneticMomentTransition = LaserInfoMat(:,2);
    BeamPropagation = LaserInfoMat(:,3:5);
    PeakIntensity =  LaserInfoMat(:,6);
    Offset = LaserInfoMat(:,7:9);
    DirWx = LaserInfoMat(:,10:12);
    DirWy = LaserInfoMat(:,13:15);
    Zr = LaserInfoMat(:,16:17);
    w = LaserInfoMat(:,18:19);
    PolarizationCR = LaserInfoMat(:,20:22);
    PolarizationCL = LaserInfoMat(:,23:25);
    AdressTransition = LaserInfoMat (:,26);
    DetuningMatrix = LaserInfoMat(:,27:end);
    
    %Atom info
    r =  AtomInfoMat(:,1:3);
    v =  AtomInfoMat(:,4:6);
    ZeemanSublevel =  AtomInfoMat(:,8);
    
    %ZeemanShifts 
%     auxVec = [1,-1,0];
%     MagneticField =G.*r;
%     ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);                
    [ZeemanDetunings ,MagneticField] =  getZeemanShift(G,r,magneticMomentTransition,ZeemanSublevel,PolinomialFitMOTfield);

    %DopllerShift

    DopplerShiftMatrix = kred*(BeamPropagation* v');

    %Detunings
    Deltasp = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,1);
    Deltasm = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,2);
    Deltapi = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,3) ;
    
    %Get polarization of laser
    PolComponentssq = GetLaserPolarizationComponentsVector(BeamPropagation,MagneticField,PolarizationCR,PolarizationCL);

    %Get average saturation parameter per combline
    Intensity = getlaserIntesities(r,Offset,DirWx,DirWy,BeamPropagation,PeakIntensity, Zr,w);  
    s0 =  Intensity./NumberOfcomblines/Isat;

    %Get CG coeficients
    n = (ZeemanSublevel*2-1)/2+6;
    TransitionStrenghtMAt = (AdressTransition == 1)* CGtableTrapping(n,:).^2  + (AdressTransition== 0)* CGtableSitrring(n,:).^2;
%     s0eff = s0*TransitionStrenght.*PolComponentssq;
    s0eff = bsxfun(@times,TransitionStrenghtMAt.*PolComponentssq,s0);


    %Calculate scattering rate.
    PSctsp = sum((s0eff(:,1)/2)./(1+s0eff(:,1)+ (2*Deltasp/GammaRed).^2),2);
    PSctsm = sum((s0eff(:,2)/2)./(1+s0eff(:,2)+ (2*Deltasm/GammaRed).^2),2);
    PSctpi = sum((s0eff(:,3)/2)./(1+s0eff(:,3)+ (2*Deltapi/GammaRed).^2),2);
    PSctsp(PSctsp>0.5) = 0.5;
    PSctsm(PSctsm>0.5) = 0.5;
    PSctpi(PSctpi>0.5) = 0.5;
    PSctsp = GammaRed*PSctsp;
    PSctsm = GammaRed*PSctsm;
    PSctpi = GammaRed*PSctpi;
  end
 
  function intensityBeamMatrix = getlaserIntesities(r,Offset,DirWx,DirWy,BeamPropagation,peakIntensity, Zr,w)   

    r = r-Offset;
    x = r.*DirWx;
    y = r.*DirWy;
    z = r.*BeamPropagation;          

    I0 =peakIntensity;

    zrx =Zr(:,1);
    zry =Zr(:,2);
    
    wx = w(:,1);
    wy = w(:,2);

    zsq = sum(z.*z,2);
    xsq = sum(x.*x,2);
    ysq = sum(y.*y,2);

    omegax = wx.*sqrt(1+(zsq./(zrx.^2)) );
    omegay = wx.*sqrt(1+(zsq./(zry.^2)) );

    intensityBeamMatrix =  I0.*(wx./omegax).*(wy./omegay).*exp(-2*(xsq)./(omegax.^2)).*exp(-2*(ysq)./(omegay.^2));


  end 
  function PolComponents = GetLaserPolarizationComponentsVector(BeamDirection,Bfield,PolarizationArrayCR,PolarizationArrayCL)

      Bdir = Bfield/norm(Bfield);
      LaserDir = BeamDirection;
      costeta = LaserDir*Bdir';
      sintetasq = 1-costeta.^2;

      PolComponents = (PolarizationArrayCR).*[(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2] + (PolarizationArrayCL).*[((1-costeta).^2)/4,((1+costeta).^2)/4,sintetasq/2];
%             

  end

  function [ZeemanDetunings,MagneticField] =  getZeemanShift(G,r,magneticMomentTransition,ZeemanSublevel,PolinomialFitMOTfield)
      %ZeemanShifts 
      if (r(2)<0.1)
          
        auxVec = [1,-1,0];
        MagneticField =G.*r;
        ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);  
      else
         auxVec = [1,-1,0];
         
         MagneticField =G.*r;
%          MagneticField(2) = 0;
%          MagneticField(2) = polyval(PolinomialFitMOTfield,r(2));
        ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);  
       end
  end
  
  function [ZeemanDetunings,MagneticField] =  getZeemanShiftv2(G,r,magneticMomentTransition,ZeemanSublevel,PolinomialFitMOTfield)
      %ZeemanShifts 
      if (r(2)>0)
          
        auxVec = [1,-1,0];
        MagneticField =G.*r;
        MagneticField(2) = G(2).*r(2).*exp(-r(2).^2/0.15^2);
        ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);  
      else
         auxVec = [1,-1,0];
         
         MagneticField =G.*r;
%          MagneticField(2) = 0;
%          MagneticField(2) = polyval(PolinomialFitMOTfield,r(2));
        ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);  
       end
  end
  