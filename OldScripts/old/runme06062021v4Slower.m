
tic 
mfvector =-9/2:9/2;
reptimes = 50;
stirringonBoleanArray= [false,true];
v_vector = linspace(1,4,20);

%Create progress bar
steps = length(v_vector);
h = waitbar(0,'Please wait runing simulation');
k= 0;
for velocity = v_vector
    k = k + 1;
    waitbar(k/steps)
    for mf = mfvector
        for stirringonboolean = stirringonBoleanArray
                %Get parameters
                params = ParamClass;
                params.SavePath = 'Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\r8\';


                %initialize conditions
                r0 = [0,250,0]*1e-3;
    %             v0 = [0,-4500,0]*1e-3;
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
                tic
            parfor dummyvar = 1:reptimes
                tic
                runSIM(params,laserInfo,atomInfoMat,cteArray,CGTable,'SlowerSIM')
                toc
            end
        end

    end
end
% runme06062021v4Slower