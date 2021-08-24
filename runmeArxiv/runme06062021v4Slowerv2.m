
tic 
mfvector =-9/2:9/2;
reptimes = 500;
stirringonBoleanArray= [false,true];
v_vector = linspace(1,4,20);
%Create progress bar
steps = length(v_vector);
h = waitbar(0,'Please wait runing simulation');
k= 0;
simSpace = [size(v_vector,2),size(mfvector,2),size(stirringonBoleanArray,2),reptimes];
numSim = prod(simSpace)
for indx1 = 1:length(v_vector)
    for indx2 =1:length( mfvector)
        for indx3 = 1:length(stirringonBoleanArray)
            for indx4 = 1:reptimes
                stirringonboolean = stirringonBoleanArray(indx3);
                velocity = v_vector(indx1);
                mf = mfvector(indx2);
                waitbar(indx1/steps);
                %Get parameters
                params = ParamClass;
                params.SavePath = 'Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\r10\';
                %initialize conditions
                r0 = [0,250,0]*1e-3;
                v0 = [0,-velocity,0];

                atom = ParticleClass(params,r0,v0,mf);
                atomarray = [atom];


                %Create lasers

                YCaptureStirring = LaserClass('YcaptureStirring');
                XCaptureStirringPositive = LaserClass('XCaptureStirringPositive');
                XCaptureStirringNegative = LaserClass('XCaptureStirringNegative');

                YCapture = LaserClass('Ycapture');
                XCapturePositive = LaserClass('XcapturePositive');
                XCaptureNegative = LaserClass('XcaptureNegative');
                ZCapturePositive = LaserClass('ZcapturePositive');
                ZCaptureNegative = LaserClass('ZcaptureNegative');
                if(stirringonboolean)
                    LaserArray = [YCapture,XCapturePositive,XCaptureNegative,ZCapturePositive,ZCaptureNegative,YCaptureStirring,XCaptureStirringPositive,XCaptureStirringNegative];
                else

                    LaserArray = [YCapture,XCapturePositive,XCaptureNegative,ZCapturePositive,ZCaptureNegative];
    %                 LaserArray = [YCaptureStirring,XCaptureStirringPositive,XCaptureStirringNegative];


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
                
                
                simInfo{indx1,indx2,indx3,indx4} = {params,laserInfo,atomInfoMat,cteArray,CGTable,'SlowerSIM'};
            end
        end

    end
end
close (h)
parfor run_counter = 1:numSim  
    [indx1,indx2 , indx3,indx4 ] = ind2sub(simSpace, run_counter);
    runSIM(simInfo{indx1,indx2,indx3,indx4}{:},run_counter)
    
end

Data_Slower = DataClass("Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\r10\","Data_2.mat")

% runme06062021v4Slower