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
atomarray = [atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom,atom];
atomarray = [atomarray,atomarray,atomarray,atomarray,atomarray,atomarray]
% atomarray = [atomarray,atomarray,atomarray,atomarray,atomarray,atomarray,atomarray,atomarray,atomarray,atomarray,atomarray]

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
laserInfogpu = 1*gpuArray(MOTLasersFormat2.LaserInfoMat);
atomInfoMatgpu = 1*gpuArray(MOTLasersFormat2.AtomInfoMat);
cteArraygpu = 1*gpuArray(params.cteArray);
CGTablegpu =  1*gpuArray(params.CGTableCombined);
auxarrayGPU = gpuArray([1,-1,0])

laserInfo = 1*(MOTLasersFormat2.LaserInfoMat);
atomInfoMat = 1*(MOTLasersFormat2.AtomInfoMat);
cteArray = params.cteArray;
CGTable =  params.CGTableCombined;
natom = length (atomarray)
nlasers = length(LaserArray)
% rng default
n = 1;
%  tic
%   [vx,vy,vz,~,~] = getVelocityChangev2ngpu(laserInfo,atomInfoMat(n,:),cteArray,CGTable);


    tic
%     [Psp,Psm,Pp] =arrayfun(@(n)ScatteringRatev3gpu(CGTablegpu,cteArraygpu,laserInfogpu(n,:),atomInfoMatgpu),1:nlasers);
    
    arrayfun(@(n)ScatteringRatev3gpu(CGTablegpu,cteArraygpu,laserInfogpu,atomInfoMatgpu(n,:),auxarrayGPU),1:natom);
    t2 = toc
    2*2
    tic
     arrayfun(@(n)ScatteringRatev3(CGTable,cteArray,laserInfo,atomInfoMat(n,:)),1:natom);
   t1 =  toc 
    

    t2/t1
% % % 
%  rng default
% % % tic
% % %  [vx2,vy2,vz2,~,~] = getVelocityChangev4(laserInfogpu,atomInfoMatgpu(1,:),cteArraygpu,CGTablegpu);
% % % toc
% assert(vx ==vx2)
% assert(vy==vy2)
% assert(vz == vz2)
end
