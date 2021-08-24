r0 = [0,0,10]*10^-3;
v0 = [0,0,0]*10^-3 ;
mf = 9/2;
params = ParamClass;
ZCaptureNegative = LaserClass('ZcaptureNegative');

%% Laser Class: get polarization components
Direction = [0,0,1];
expected = [0,1,0]
polComponets = ZCaptureNegative.getPolarizationComponents(Direction)
assert(isequal (polComponets,expected))

%% Laser Class: get polarization components
Direction = [0,0,-1];
expected = [1,0,0]
polComponets = ZCaptureNegative.getPolarizationComponents(Direction)
assert(isequal (polComponets,expected))

%% Laser Class: get polarization components
Direction = [0,0,-1];
ZCaptureNegative.Polarization = "CL"
expected = [0,1,0]
polComponets = ZCaptureNegative.getPolarizationComponents(Direction)
assert(isequal (polComponets,expected))

%% Laser Class: get polarization components
Direction = [0,1,0];
ZCaptureNegative.Polarization = "CL"
expected = [1/2,1/2,1/sqrt(2)]
polComponets = ZCaptureNegative.getPolarizationComponents(Direction)
assert(isequal (polComponets,expected))

%% Laser Class: get polarization components
Direction = [0,1,0];
ZCaptureNegative.Polarization = "CR"
expected = [1/2,1/2,1/sqrt(2)]
polComponets = ZCaptureNegative.getPolarizationComponents(Direction)
assert(isequal (polComponets,expected))

%% Laser Class: Intensity Profile
r = [0,ZCaptureNegative.wy,0]+ZCaptureNegative.offset;
intesity = ZCaptureNegative.LaserIntensityProfile(r)
expected = ZCaptureNegative.peakIntensity*exp(-2)
assert((intesity - expected) <10^-15)

%% Laser Class: Intensity Profile
r = [ZCaptureNegative.wx,ZCaptureNegative.wy,0]+ZCaptureNegative.offset;
intesity = ZCaptureNegative.LaserIntensityProfile(r)
expected = ZCaptureNegative.peakIntensity*exp(-4)
assert(abs(intesity - expected) <10^-15)
%% Laser Class: get polarization components
Direction = rand(1,3);
ZCaptureNegative.Polarization = "CL"
pol1 =  ZCaptureNegative.getPolarizationComponents(Direction)
pol2 = ZCaptureNegative.getPolarizationComponents(Direction*100)
assert(norm(abs(pol1 - pol2)) <10^-15)

%% Laser Class: get polarization componentssq
for i =  1 : 100
    Direction = rand(1,3);
    ZCaptureNegative.Polarization = "CL"
    pol1 =  ZCaptureNegative.getPolarizationComponents(Direction).^2
    pol2 = ZCaptureNegative.getPolarizationComponentssq(Direction)
    assert(norm(abs(pol1 - pol2)) <10^-15)
end