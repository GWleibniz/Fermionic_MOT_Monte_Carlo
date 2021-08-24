
tic 
numsims = 1e5;
savePoints = 100;%Save data every 100 points;
counter =0;
%Get parameters
params = ParamClass;


%initialize conditions
r0 = [0,-12.6,0]*1e-3;
v0 = [0,0,0]*1e-3;
mf  = -9/2
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
LaserArray = [YCapture,XCapturePositive,XCaptureNegative,ZCapturePositive,ZCaptureNegative,YCaptureStirring,XCaptureStirringPositive,XCaptureStirringNegative]

%Format objects to make the code esailly vectoraizable.
MOTLasersFormat = FormatInputsv2(LaserArray,atomarray);
% MOTLasersFormat = FormatInputsgpu(LaserArray);
laserInfo = (MOTLasersFormat.LaserInfoMat);
atomInfoMat = (MOTLasersFormat.AtomInfoMat);
cteArray = params.cteArray;
CGTable =  (params.CGTableCombined);
params.LaserArray = LaserArray;




%Decay probability
atomInfoMat2Save = zeros(ceil(numsims/savePoints),8);
%Create an atom
tic
for i = 1:numsims
%        vchange = getVelocityChange(params,MOTLasersFormat,atom);
       [DeltaV,atomicstate,atomZeemanSublevel] = getVelocityChangev2(laserInfo,atomInfoMat,cteArray,CGTable);
       [DeltaV2,atom] = getVelocityChange(params,MOTLasersFormat,atom);       
       
       
       DeltaV = DeltaV + [0,-params.Gravity,0]*params.UpdateStepSize  ;
       atom = atom.UpdateLocation();
       atom = atom.UpdateVelocity(DeltaV);
       Deltar = atomInfoMat(4:6)*params.UpdateStepSize;
       
       atomInfoMat(1:3) = atomInfoMat(1:3) + Deltar;
       atomInfoMat(4:6) = atomInfoMat(4:6) + DeltaV;
       atomInfoMat(7) = atomicstate;
       atomInfoMat(8) = atomZeemanSublevel;

      assert(prod(DeltaV2 == DeltaV2)==1)

       if rem(i,savePoints)==0
         counter =counter+1;
         atomInfoMat2Save(counter,:) = atomInfoMat;
         atomInfoMat(1:3)*1000
       end

end
toc


