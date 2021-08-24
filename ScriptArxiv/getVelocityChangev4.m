function [vx,vy,vz,atomicstate,ZeemanSublevel]   =  getVelocityChangev4(LaserInfoMat,AtomInfoMat,cteArray,CGTable)%Initialize variabl

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
