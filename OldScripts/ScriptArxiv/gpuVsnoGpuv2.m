params = ParamClass;
iterations =  1:100000;
rr = rand(length(iterations),3)*1e-3-13*1e-3;
vv = rand(length(iterations),3)*1e-3;
for i = iterations


%initial conditions

r0 =rr(i,:);
v0 = vv(i,:);

mf  = -9/2
atom = ParticleClass(params,r0,v0,mf);
atomarray = [atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom];


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
MOTLasersFormat = FormatInputs(LaserArray,atomarray);

% MOTLasersFormat = FormatInputsgpu(LaserArray);

MOTLasersFormat2 = FormatInputsv2(LaserArray,atomarray);
MagneticMoment = MOTLasersFormat2.MagneticMomentAdresseTransitionArray;


atomicstate = atom.AtomicState;
atomZeemanSublevel = atom.ZeemanSublevel;
r = atom.r;
v = atom.v;
%% Particle Class: UpdateLocation
laserInfogpu = gpuArray(MOTLasersFormat2.LaserInfoMat);
atomInfoMatgpu = gpuArray(MOTLasersFormat2.AtomInfoMat);
cteArraygpu = gpuArray(params.cteArray);
CGTablegpu =  gpuArray(params.CGTableCombined);


laserInfo = (MOTLasersFormat2.LaserInfoMat);
atomInfoMat = (MOTLasersFormat2.AtomInfoMat);
cteArray = params.cteArray;
CGTable =  params.CGTableCombined;
 
% rng default
 tic
  [vx,vy,vz,~,~] =arrayfun(@(n) getVelocityChangev2(laserInfo,atomInfoMat(n,:),cteArray,CGTable),1:1);
 toc 
% % % 
%  rng default
tic
 [vx2,vy2,vz2,~,~] = arrayfun(@(n) getVelocityChangev2(laserInfogpu,atomInfoMatgpu(1,:),cteArraygpu,CGTablegpu),1:1);
toc
% assert(vx ==vx2)
% assert(vy==vy2)
% assert(vz == vz2)
end
