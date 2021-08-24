function [vx,vy,vz,atomicstate,ZeemanSublevel]   =  getVelocityChangev2ngpu(LaserInfoMat,AtomInfoMat,cteArray,CGTable)%Initialize variabl

    UpdateStepsize = cteArray(7);
    vrecoil =cteArray(6);
    
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
        %Calculate scattering rate of all the lasers        
%          ScatteringRateMatrix  =cell2mat(arrayfun(@(Laser)ScatteringRate(atom,Laser),LaserArray','UniformOutput',false));
%         [ScatteringRateMatrix] = ScatteringRatev2(params,MOTLasersFormat,MagneticMoment,ZeemanSublevel,r,v);
%         [Psp,Psm,Spp] =  ScatteringRatev3(params,LaserInfoMat,ZeemanSublevel,r,v)
%         ScatteringRateMatrix = [Psp,Psm,Spp];
        nlasers = size(LaserInfoMat,1);


        
        [Psp,Psm,Pp] = arrayfun(@(n)ScatteringRatev3(CGTable,cteArray,LaserInfoMat(n,:),AtomInfoMat),1:nlasers);
        
        ScatteringRateMatrix = [Psp',Psm',Pp'];
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

        if (rand() < Psc*100000 ) %Atom scatter a photon and goes to the excited AtomicState
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
    vx =  deltaVLasers(1);
    vy =  deltaVLasers(2); 
    vz =  deltaVLasers(3);
end

  function [PSctsp,PSctsm,PSctpi] = ScatteringRatev3(CGTable,cteArray,LaserInfoMat,AtomInfoMat)
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
    AtomicState =  AtomInfoMat(:,7);
    ZeemanSublevel =  AtomInfoMat(:,8);
    
    %ZeemanShifts 
    auxVec = [1,-1,0];
    MagneticField =G.*r;
    ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);                


    %DopllerShift

    DopplerShiftMatrix = kred*(BeamPropagation*v');


    %Detunings
    Deltasp = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,1);
    Deltasm = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,2);
    Deltapi = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,3) ;

      Bdir = MagneticField/norm(MagneticField);
      LaserDir = BeamPropagation;
      costeta = LaserDir*Bdir';
      sintetasq = 1-costeta.^2;

      PolComponents = (PolarizationCR).*[(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2] + (PolarizationCL).*[((1-costeta).^2)/4,((1+costeta).^2)/4,sintetasq/2];
%             

    %Get polarization of laser
    PolComponentssq = GetLaserPolarizationComponentsVector(BeamPropagation,MagneticField,PolarizationCR,PolarizationCL);

    %Get average saturation parameter per combline
    Intensity = getlaserIntesitiesv3(r,Offset,DirWx,DirWy,BeamPropagation,PeakIntensity, Zr,w);  
    s0 =  Intensity./NumberOfcomblines/Isat;

    %Get CG coeficients
    n = (ZeemanSublevel*2-1)/2+6;
    TransitionStrenght = CGtableTrapping(n,:).^2;
    TransitionStrenghtMAt = (AdressTransition == 1)* CGtableTrapping(n,:).^2  + (AdressTransition== 0)* CGtableSitrring(n,:).^2;
%     s0eff = s0*TransitionStrenght.*PolComponentssq;
    s0eff = bsxfun(@times,TransitionStrenghtMAt.*PolComponentssq,s0);


    %Calculate scattering rate.
    PSctsp = sum(GammaRed*(s0eff(:,1)/2)./(1+s0eff(:,1)+ (2*Deltasp/GammaRed).^2),2);
    PSctsm = sum(GammaRed*(s0eff(:,2)/2)./(1+s0eff(:,2)+ (2*Deltasm/GammaRed).^2),2);
    PSctpi = sum(GammaRed*(s0eff(:,3)/2)./(1+s0eff(:,3)+ (2*Deltapi/GammaRed).^2),2);


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
  function intensityBeamMatrix = getlaserIntesitiesv3(r,Offset,DirWx,DirWy,BeamPropagation,peakIntensity, Zr,w)   

    r = r-Offset;
    x = r.*DirWx;
    y = r.*DirWy;
    z = r.*BeamPropagation;          

    I0 =peakIntensity;

    zrx =Zr(:,1);
    zry =Zr(:,2);
    
    wx = w(:,1);
    wy = w(:,2);


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
  function PolComponents = GetLaserPolarizationComponentsVector(BeamDirection,Bfield,PolarizationArrayCR,PolarizationArrayCL)

      Bdir = Bfield/norm(Bfield);
      LaserDir = BeamDirection;
      costeta = LaserDir*Bdir';
      sintetasq = 1-costeta.^2;

      PolComponents = (PolarizationArrayCR).*[(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2] + (PolarizationArrayCL).*[((1-costeta).^2)/4,((1+costeta).^2)/4,sintetasq/2];
%             

end
