
tic 
reptimes = 100;
stirringonBoleanArray= [false];
XCapPower = linspace(0,7,10)*1e-3;
simSpace = [size(XCapPower,2),reptimes];
numSim = prod(simSpace)
for indx = 1:length(XCapPower)
    for indx23 = 1:reptimes
        stirringonboolean = false
        numsims = 1e5;
        mf = -4.5
        repetitions =1;
        savePoints = 1;%Save data every 100 points;
        counter =0;

        %Get parameters
        params = ParamClass;
        params.StopTime = 10
        params.SavePath = 'Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\MOT\raw\run_data_two_changing_x_cap\';

        %initialize conditions
        r0 = [0,-13,0]*1e-3;
        v0 = [randn(),0,randn()*20]*1e-3;
        atom = ParticleClass(params,r0,v0,mf);
        atomarray = [atom];


        %Create lasers

        YCaptureStirring = LaserClass('YcaptureStirring');
        XCaptureStirringPositive = LaserClass('XCaptureStirringPositive');
        XCaptureStirringNegative = LaserClass('XCaptureStirringNegative');

        YCapture = LaserClass('Ycapture');
        XCapturePositive = LaserClass('XcapturePositive');
        XCaptureNegative = LaserClass('XcaptureNegative');
        XCapturePositive.Power = XCapPower;
        XCaptureNegative.Power = XCapPower;
        ZCapturePositive = LaserClass('ZcapturePositive');
        ZCaptureNegative = LaserClass('ZcaptureNegative');
        if(stirringonboolean)
            LaserArray = [YCapture,XCapturePositive,XCaptureNegative,ZCapturePositive,ZCaptureNegative,YCaptureStirring,XCaptureStirringPositive,XCaptureStirringNegative];
        else
            LaserArray = [YCapture,XCapturePositive,XCaptureNegative,ZCapturePositive,ZCaptureNegative];

        end
        %Format objects to make the code esailly vectoraizable.
        MOTLasersFormat = FormatInputsv2(LaserArray,atomarray);
        % MOTLasersFormat = FormatInputsgpu(LaserArray);
        laserInfo = (MOTLasersFormat.LaserInfoMat);
        atomInfoMat = (MOTLasersFormat.AtomInfoMat);
        cteArray = params.cteArray;
        CGTable =  (params.CGTableCombined);
        params.LaserArray = LaserArray;
        params.atomArray = atomarray;



        atomSlowedBoolean = 0;
        atomLostBolean = 0;
        i= 0;
        simInfo{indx,indx23} = {params,laserInfo,atomInfoMat,cteArray,CGTable,'MOTSim'};
        tic
    end
end

    parfor run_counter = 1:numSim  
        tic
        [indx1,indx2] = ind2sub(simSpace, run_counter);
        runSIM(simInfo{indx,indx2}{:},run_counter)
        toc
    end

