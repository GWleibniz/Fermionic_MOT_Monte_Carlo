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
MOTLasersFormat = FormatInputs(LaserArray,atomarray);

% MOTLasersFormat = FormatInputsgpu(LaserArray);

MOTLasersFormat2 = FormatInputsv2(LaserArray,atomarray);
MagneticMoment = MOTLasersFormat2.MagneticMomentAdresseTransitionArray;

r = atom.r
v = atom.v
%% Particle Class: UpdateLocation
rng default
laserInfo = (MOTLasersFormat2.LaserInfoMat);
atomInfoMat = (MOTLasersFormat2.AtomInfoMat);
cteArray = params.cteArray;
CGTable =  (params.CGTableCombined);


 [ DeltaV,atomicstate,atomZeemanSublevel] = getVelocityChangev2(laserInfo,atomInfoMat,cteArray,CGTable);
 vx= DeltaV(1);
 vy= DeltaV(2);
 vz= DeltaV(3);
 rng default
[vx2,vy2,vz2] = getVelocityChange(params,MOTLasersFormat2,atom);
assert(vx ==vx2)
assert(vy==vy2)
assert(vz == vz2)
atom.AtomicState = atomicstate
atomZeemanSublevel = atomZeemanSublevel
atom.r = [r+v*1e-6]
atom.v = v + [vx,vy,vz]
end
