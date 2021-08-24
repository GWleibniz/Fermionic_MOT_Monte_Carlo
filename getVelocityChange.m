function [deltaVLasers,atom]   =  getVelocityChange(params,MOTLasersFormat,atom)%Initialize variables
    deltaVLasers = (zeros(1,3));
    zerosv = (zeros(1,3));
    auxVec = [-1,1,0];
    deltaJtransition = [1,-1,0];
    Pdecay = params.GammaRed*params.UpdateStepSize;
    
    %if we are in the ground state
    if (atom.AtomicState == 0)
        %Calculate scattering rate of all the lasers        
%          ScatteringRateMatrix  =cell2mat(arrayfun(@(Laser)ScatteringRate(atom,Laser),LaserArray','UniformOutput',false));
        [ScatteringRateMatrix] = ScatteringRatev2(atom,MOTLasersFormat);
        ProbabilityVectorScattering = params.UpdateStepSize*reshape(ScatteringRateMatrix',[],prod(size(ScatteringRateMatrix)));

        %Use a random number to find wich laser will scatter
        ProbabilityVector = cumsum(ProbabilityVectorScattering/sum(ProbabilityVectorScattering));
        probabilityNumber = rand();  
        [~,selected] = max(probabilityNumber<ProbabilityVector);

        %Selected Laser
        LaserScatered = MOTLasersFormat.LaserArray( ceil(selected/3));
        Direction = MOTLasersFormat.BeamPropagationMatrix(ceil(selected/3),:);

        %Probability to scatter from the selected laser
        Psc = ProbabilityVectorScattering(selected);

        %Transition type
%         TransitionType2 = rem(selected+2,3)+1;% 1 = sp, 2 = sm , 0 = pi;
        TransitionType = selected + 3 - floor((selected+2)./3)*3;
%         gather(ScatteringRateMatrix)

        if (rand() < Psc ) %Atom scatter a photon and goes to the excited AtomicState
            if (LaserScatered.IsTrappingLaser)
                atom.AtomicStateExcited = 11/2;
                atom.AtomicState = 1;
            else
                atom.AtomicStateExcited = 9/2;
                atom.AtomicState = 1;
            end

           atom.ZeemanSublevel = atom.ZeemanSublevel + deltaJtransition(TransitionType);
           deltaVLasers = params.vrecoil*Direction;
        else
           deltaVLasers = zerosv;
        end



    %if we are on one of the excited states
    elseif ( atom.AtomicState ~=0 )
        if (rand() < Pdecay ) %Atom  emits a photon and goes to the ground AtomicStat

            %Calculate relative transition strenght of the possible,sp sm and pi decay paths
            Fp = atom.AtomicStateExcited;
            if (Fp ==5.5)
                CGCoeffDecayTable = params.CGTable.TrappingTransitionDecay ;
                n = (atom.ZeemanSublevel*2-1)/2+7;
            elseif ( Fp == 4.5)
                CGCoeffDecayTable =  params.CGTable.StirringTransitionDecay ;
                n = (atom.ZeemanSublevel*2-1)/2+6;
            end

            vectorDecay = CGCoeffDecayTable(n,:).^2;


            %Select decay path using a randon number ( same as above)
            ProbabilityVector = cumsum(vectorDecay/sum(vectorDecay));
            probabilityNumber = rand();
            [~,selected] = max(probabilityNumber<ProbabilityVector);


            %Calculate change of Zeeman sublevel based on decay path

            atom.ZeemanSublevel = atom.ZeemanSublevel + auxVec(selected);

            %Bring the atom back to the ground state
            atom.AtomicState = 0;

            %Momentum kick in random direction
            v = randn(1,3);
            deltaVLasers = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)))*params.vrecoil;

        else
            deltaVLasers = zerosv;

        end

    end

end

  function [Psc] = ScatteringRatev2(particle,MOTlasers)
            
    G = particle.params.BGradient;
    kred =   particle.params.kred;      
    MagneticMoment3P1 =particle.params.MagneticMoment3P1;
    CGtableTrapping =   particle.params.CGTable.TrappingTransition;
    CGtableSitrring =   particle.params.CGTable.StirringTransition;

    Isat = particle.params.Isat;
    GammaRed = particle.params.GammaRed;
    MagneticMomentArray = MOTlasers.MagneticMomentAdresseTransitionArray;

    DetuningMatrix = MOTlasers.LaserDetuningMatrix;
    BeamPropagationDirectionMatrix = MOTlasers.BeamPropagationMatrix;
    PolarizationArray = MOTlasers.polarizationArray;

    %ZeemanShifts 
    auxVec = [1,-1,0];
    MagneticField =G.*particle.r;
    ZeemanDetunings =MagneticMomentArray *vecnorm(MagneticField,2,2)*(particle.ZeemanSublevel+auxVec);                


    %DopllerShift
    DopplerShiftMatrix = kred*(BeamPropagationDirectionMatrix*particle.v');


    %Detunings
    Deltasp = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,1);
    Deltasm = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,2);
    Deltapi = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,3) ;

    %Get polarization of laser
    PolComponentssq = GetLaserPolarizationComponentsVector(BeamPropagationDirectionMatrix,MagneticField,PolarizationArray);

    %Get average saturation parameter per combline
    Intensity =  getlaserIntesitiesv2(particle.r,MOTlasers)   ;
    s0 =  Intensity./MOTlasers.NumberoffCombLinesArray/Isat;

    %Get CG coeficients
    n = (particle.ZeemanSublevel*2-1)/2+6;
    TransitionStrenght = CGtableTrapping(n,:).^2;
    TransitionStrenghtMAt = (MOTlasers.adressedTransitionArray == 1)* CGtableTrapping(n,:).^2  + (MOTlasers.adressedTransitionArray == 0)* CGtableSitrring(n,:).^2;
%     s0eff = s0*TransitionStrenght.*PolComponentssq;
    s0eff = bsxfun(@times,TransitionStrenghtMAt.*PolComponentssq,s0);


    %Calculate scattering rate.
    PSctsp = sum(GammaRed*(s0eff(:,1)/2)./(1+s0eff(:,1)+ (2*Deltasp/GammaRed).^2),2);
    PSctsm = sum(GammaRed*(s0eff(:,2)/2)./(1+s0eff(:,2)+ (2*Deltasm/GammaRed).^2),2);
    PSctpi = sum(GammaRed*(s0eff(:,3)/2)./(1+s0eff(:,3)+ (2*Deltapi/GammaRed).^2),2);
    Psc = [PSctsp,PSctsm,PSctpi];

  end

  function intensityBeamMatrix = getlaserIntesitiesv2(r,MOTLasers)   

    r = r-MOTLasers.offsetarray;
    x = r.*MOTLasers.waistDirectionArrayx;
    y = r.*MOTLasers.waistDirectionArrayy;
    z = r.*MOTLasers.BeamPropagationMatrix;          

    I0 = MOTLasers.PeakIntesityArray;

    zrx = MOTLasers.RayleightLegnthArray(:,1);
    zry = MOTLasers.RayleightLegnthArray(:,2);

    wx = MOTLasers.WaistArray(:,1);
    wy = MOTLasers.WaistArray(:,2);


%             zsq = vecnorm (z,2,2).^2;
%             xsq = vecnorm (x,2,2).^2;
%             ysq = vecnorm (y,2,2).^2;
    zsq = sum(z.*z,2);
    xsq = sum(x.*x,2);
    ysq = sum(y.*y,2);

    omegax = wx.*sqrt(1+(zsq./(zrx.^2)) );
    omegay = wx.*sqrt(1+(zsq./(zry.^2)) );

    intensityBeamMatrix =  I0.*(wx./omegax).*(wy./omegay).*exp(-2*(xsq)./(omegax.^2)).*exp(-2*(ysq)./(omegay.^2));


  end  

  function PolComponents = GetLaserPolarizationComponentsVector(BeamPropagationDirectionMatrix,Bfield,PolarizationArray)

      Bdir = Bfield/norm(Bfield);
      LaserDir = BeamPropagationDirectionMatrix;
      costeta = LaserDir*Bdir';
      sintetasq = 1-costeta.^2;

      PolComponents = (PolarizationArray=="CR").*[(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2] + (PolarizationArray=="CL").*[((1-costeta).^2)/4,((1+costeta).^2)/4,sintetasq/2];
%             

end
