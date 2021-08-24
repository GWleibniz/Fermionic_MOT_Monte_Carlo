

mfvector = -5/2:1:9/2


        numsims = 001000000;

        %initial conditions
        
        mf = -9/2
        r0 = [0,-13,0]*1e-3;
        v0 = [0,-0,0]*1e-3;

        %Create lasers
        YCaptureStirring = LaserClass("YcaptureStirring");
        XCaptureStirringPositive = LaserClass("XCaptureStirringPositive");
        XCaptureStirringNegative = LaserClass("XCaptureStirringNegative");

        YCapture = LaserClass("Ycapture");
        XCapturePositive = LaserClass("XcapturePositive");
        XCaptureNegative = LaserClass("XcaptureNegative");
        ZCapturePositive = LaserClass("ZcapturePositive");
        ZCaptureNegative = LaserClass("ZcaptureNegative");
        LaserArray = [YCapture,XCapturePositive,XCaptureNegative,ZCapturePositive,ZCaptureNegative]
        
        
%         LaserArray = [YCapture,XCapturePositive,XCaptureNegative]
        NumberOfLasers = length(LaserArray);
        MOTLasers = MOTClass(LaserArray);
        [a,maxnumberofComblines] = size(MOTLasers.LaserDetuningMatrix);
        params = ParamClass;
        params.SavePath = "C:\Users\PC de Rodrigo\Documents\MATLAB\save"
        params.LaserArray = LaserArray;



        %aux stuff
        auxVec = [-1,1,0];
        deltaJtransition = [1,-1,0];
        %Decay probability
        Pdecay = params.GammaRed*params.UpdateStepSize ;





        atomSlowedBoolean = false;
        atomLostBoolean   = false;
        counter =0;
        counterBad = 0;
        %Create an atom
        atom = ParticleClass(params,r0,v0,mf);
        r = []
for i = 1:numsims
    deltaVLasers = 0;
    %if we are in the ground state
    if (atom.AtomicState == 0)

        %Calculate scattering rate of all the lasers         
        [PSctMat2] = atom.ScatteringRatev2(MOTLasers.LaserDetuningMatrix,MOTLasers.BeamPropagationMatrix,MOTLasers.polarizationArray,MOTLasers);
        ProbabilityVectorScatteringv2 =  params.UpdateStepSize*reshape(PSctMat',[],prod(size(PSctMat)));
        ScatteringRateMatrix  =cell2mat(arrayfun(@(Laser)ScatteringRate(atom,Laser),LaserArray','UniformOutput',false));
        ProbabilityVectorScattering = params.UpdateStepSize*reshape(ScatteringRateMatrix',[],prod(size(ScatteringRateMatrix)));        

        %Use a random number to find wich laser will scatter
        probabilityNumber = rand();
        ProbabilityVector = cumsum(ProbabilityVectorScatteringv2/sum(ProbabilityVectorScatteringv2));
        [~,selected] = max(probabilityNumber<ProbabilityVector);

        %Selected Laser
        LaserScatered = LaserArray( ceil(selected/3));

        Direction = LaserScatered.BeamPropagationDirection;
        Psct = ProbabilityVectorScattering(selected);
        TransitionType = rem(selected+2,3)+1;% 1 = sp, 2 = sm , 0 = pi;

        %End vectorized code

               ScatteringRateMatrix  =cell2mat(arrayfun(@(Laser)ScatteringRate(atom,Laser),LaserArray','UniformOutput',false));
              ProbabilityVectorScattering = params.UpdateStepSize*reshape(ScatteringRateMatrix',[],prod(size(ScatteringRateMatrix)));
% % % % % %                 
% % % % % %                 %Use a random number to find wich laser will scatter
% % % % % %                 ProbabilityVector = cumsum(ProbabilityVectorScattering/sum(ProbabilityVectorScattering));
% % % % % %                 [~,selected] = max(probabilityNumber<ProbabilityVector);
% % % % % %                 
% % % % % %                 %Selected Laser
% % % % % %                 LaserScatered = LaserArray( ceil(selected/3));
% % % % % %                 
% % % % % %                 Direction = LaserScatered.BeamPropagationDirection;
% % % % % %                 Psct = ProbabilityVectorScattering(selected);
% % % % % %                 TransitionType = rem(selected+2,3)+1;% 1 = sp, 2 = sm , 0 = pi;

        if (rand() < Psct ) %Atom scatter a photon and goes to the excited AtomicState
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
           deltaVLasers = [0,0,0];
        end
    %if we are on one of the excited states
    elseif ( atom.AtomicState ~=0 )
        if (rand() < Pdecay ) %Atom  emits a photon and goes to the ground AtomicStat

            %Calculate relative transition strenght of the possible,sp sm and pi decay paths
            Fp = atom.AtomicStateExcited;
            mFp = atom.ZeemanSublevel;   

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
            deltaVLasers = [0,0,0];

        end
    end

       g = [0,-params.Gravity,0];
       vchange = deltaVLasers + g*params.UpdateStepSize  ;

       atom = UpdateVelocity (atom,vchange);
       atom = atom.UpdateLocation();  

       r(i,:) = atom.r;

%                atom = atom.storeTrayectoryData(atom.r,atom.v,i);
%                atom = atom.storeInternalStateData(atom.AtomicState,atom.AtomicStateExcited,atom.ZeemanSublevel,i);
       if (rem(i,10)==0)
%                atom.r*1000
       atom.r*1000
%                atom.ZeemanSublevel
       end

end



    



