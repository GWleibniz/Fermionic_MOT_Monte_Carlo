
tic 
numsims = 1e3;
%Get parameters
params = ParamClass;



%initial conditions
r0 = [0,-13,0]*1e-3;
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

params.LaserArray = LaserArray;




%Decay probability
Pdecay = params.GammaRed*params.UpdateStepSize ;


%Create an atom
tic
for i = 1:numsims
%        vchange = getVelocityChange(params,MOTLasersFormat,atom);
       [vx,vy,vz] = getVelocityChangev2(params,MOTLasersFormat,atom);
       
       for i = 1:length(atomarray)
            atomarray(i) = UpdateVelocity (atomarray(i),[vx(i),vy(i),vz(i)]);
       end
%        atom = atom.storeTrayectoryData(atom.r,atom.v,i);
%        atom = atom.storeInternalStateData(atom.AtomicState,atom.AtomicStateExcited,atom.ZeemanSublevel,i);
       if (rem(i,10)==0)
       atom.r*1000
       end

end
toc



