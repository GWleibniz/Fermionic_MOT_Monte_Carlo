params = ParamClass;


r0 = [0,0,10]*10^-3;
v0 = [0,0,0]*10^-3 ;
mf = 9/2;
particle = ParticleClass(params,r0,v0,mf);

ZCaptureNegative = LaserClass('ZcaptureNegative');
TestLaser = LaserClass("TestZPositive")

%% Particle Class: Initial conditions
expected =r0
assert(isequal (particle.r,expected))


%% Particle Class: Initial conditions
expected =v0
assert(isequal (particle.v,expected))


%% Particle Class: UpdateLocation
expected =r0 + v0*params.UpdateStepSize;
particle = particle.UpdateLocation()
assert(isequal (particle.r,expected))

%% Particle Class: Update UpdateVelocity
deltav = rand(1,3);
expected = v0 + deltav;
particle = particle.UpdateVelocity(deltav)
assert(isequal (particle.v,expected))


%% Particle Class: getZeemanShiftBosons 
G = params.BGradient;
particle.r = rand(1,3);
MagneticMoment3P1 =params.MagneticMoment3P1;
auxVec = [1,-1,0];
MagneticField =G.*particle.r;
Zs = -MagneticMoment3P1*norm(MagneticField);
assert(abs(Zs-particle.getZeemanShiftBosons())<1e-12)


%% Particle Class: getZeemanShiftFermions 
mf = -9/2:1:9/2
for i = mf
    G = params.BGradient;
    particle.ZeemanSublevel = i
    particle.r = rand(1,3);
    MagneticMoment3P1 =params.MagneticMoment3P1/(5.5);
    auxVec = [1,-1,0];
    MagneticField =G.*particle.r;
    Zs = -MagneticMoment3P1*norm(MagneticField)*(particle.ZeemanSublevel+[1,-1,0]);
    assert(norm(abs(Zs-particle.getZeemanShiftFermionsTrapping()))<1e-12)
end



%% Particle Class: getZeemanShiftFermionsSittirring 
mf = -9/2:1:9/2
for i = mf
    G = particle.params.BGradient;
    MagneticMoment3P1 =particle.params.MagneticMoment3P1/4.5/(4.5);
    auxVec = [1,-1,0];
    MagneticField =G.*particle.r;
    Zs = -MagneticMoment3P1*norm(MagneticField)*(particle.ZeemanSublevel+auxVec);
    
    assert(norm(abs(Zs-particle.getZeemanShiftFermionsSittirring()))<1e-12)
end

%% Particle Class: getDoppler2Shift 
mf = -9/2:1:9/2
for i = 1:10
    TestLaser = LaserClass("TestZPositive");
    v = rand(1,3);
    TestLaser.BeamPropagationDirection = v/norm(v)
    
    particle.v = rand(1,3)*100;
    dDopler = particle.getDoppler2Shift (TestLaser)
    expected = particle.params.kred*(TestLaser.BeamPropagationDirection*particle.v');
    assert(norm(abs(expected-dDopler))<1e-12)

end
%% Particle Class: ScatteringRate
%I have a floating point issue here!
TestLaser = LaserClass("TestZPositive")
particle.r = zeros(1,3)*1e-3;
particle.v = zeros(1,3);
zs = particle.getZeemanShiftFermionsTrapping

strenght(1) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,1,particle.ZeemanSublevel+1).^2;
strenght(2) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,-1,particle.ZeemanSublevel-1).^2;
strenght(3) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,0,particle.ZeemanSublevel).^2;


pol = TestLaser.getPolarizationComponentssq(particle.r.*params.BGradient)
s0 = TestLaser.LaserIntensityProfile(particle.r)/params.Isat;
s0 = s0*pol.*strenght;
Delta = zs;

expected = params.GammaRed*((s0/2)./(1+s0+(2*Delta/params.GammaRed).^2));
scatteringRateparticle = particle.ScatteringRate(TestLaser)
assert(norm(abs(scatteringRateparticle-expected))<1e-13)

%% Particle Class: ScatteringRate BField on

mf = -9/2:1:9/2
for i = mf
    TestLaser = LaserClass("TestZPositive")
    particle.r = rand(1,3)*1e-3;
    particle.v = zeros(1,3);
    particle.ZeemanSublevel = i;
    zs = particle.getZeemanShiftFermionsTrapping

    strenght(1) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,1,particle.ZeemanSublevel+1).^2;
    strenght(2) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,-1,particle.ZeemanSublevel-1).^2;
    strenght(3) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,0,particle.ZeemanSublevel).^2;


    pol = TestLaser.getPolarizationComponentssq(particle.r.*params.BGradient)
    s0 = TestLaser.LaserIntensityProfile(particle.r)/params.Isat;
    s0 = s0*pol.*strenght;
    Delta = zs;

    expected = params.GammaRed*((s0/2)./(1+s0+(2*Delta/params.GammaRed).^2));
    scatteringRateparticle = particle.ScatteringRate(TestLaser)
    assert(norm(abs(scatteringRateparticle-expected))<1e-13)
end
%% Particle Class: ScatteringRate BField on Doppler Shift on
for i = mf
TestLaser = LaserClass("TestZPositive")
particle.r = rand(1,3)*1e-3;
particle.v = rand(1,3);
particle.ZeemanSublevel = i;

zs = particle.getZeemanShiftFermionsTrapping
DeltaD = particle.getDoppler2Shift(TestLaser)

strenght(1) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,1,particle.ZeemanSublevel+1).^2;
strenght(2) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,-1,particle.ZeemanSublevel-1).^2;
strenght(3) = ClebschGordan(9/2,1,11/2,particle.ZeemanSublevel,0,particle.ZeemanSublevel).^2;


pol = TestLaser.getPolarizationComponentssq(particle.r.*params.BGradient)
s0 = TestLaser.LaserIntensityProfile(particle.r)/params.Isat;
s0 = s0*pol.*strenght;
Delta = zs-DeltaD;

expected = params.GammaRed*((s0/2)./(1+s0+(2*Delta/params.GammaRed).^2));
scatteringRateparticle = particle.ScatteringRate(TestLaser)
assert(norm(abs(scatteringRateparticle-expected))<1e-13)
end


%% Particle Class: MOT lasers intensity profile
for i = 1:100
    particle.r = rand();
    FormatObject = FormatInputsv2(TestLaser)
    expected =   particle.getlaserIntesities(particle.r,FormatObject);
    result  = FormatObject.LaserArray(1).LaserIntensityProfile(particle.r)
    assert(norm(abs(result-expected))<1e-16)
end

%% Particle Class: ScatteringRatev2
mf = -9/2:1:9/2
for j = 1:10
    for i = mf
        TestLaser = LaserClass("TestZPositive")
        particle.r =[0.001,0.12,-0.0132]*1e-3;
        particle.v = rand(1,3);
        particle.ZeemanSublevel = i;
        % 
        FormatObject = FormatInputsv2(TestLaser)
        scatteringRateparticle = particle.ScatteringRate(TestLaser)
        SctMat = particle.ScatteringRatev2(FormatObject)
        assert(norm(abs(scatteringRateparticle-SctMat))<1e-11)
    end
end


%% Particle Class: ScatteringRatev2
mf = -9/2:1:9/2
for j = 1:10
    for i = mf
        TestLaser = LaserClass("XcapturePositive")
        particle.r =[0.001,0.12,-0.0132]*1e-3;
        particle.v = rand(1,3);
        particle.ZeemanSublevel = i;
        % 
        FormatObject = FormatInputsv2(TestLaser)
        scatteringRateparticle = particle.ScatteringRate(TestLaser)
        SctMat = particle.ScatteringRatev2(FormatObject)
        assert(norm(abs(scatteringRateparticle-SctMat))<1e-11)
    end
end

%% Particle Class: ScatteringRatev2
mf = -9/2:1:9/2
for j = 1:10
    for i = mf
        TestLaser = LaserClass("ZcaptureNegative")
        particle.r =[0.001,0.12,-0.0132]*1e-3;
        particle.v = rand(1,3);
        particle.ZeemanSublevel = i;
        % 
        FormatObject = FormatInputsv2(TestLaser)
        scatteringRateparticle = particle.ScatteringRate(TestLaser)
        SctMat = particle.ScatteringRatev2(FormatObject)
        assert(norm(abs(scatteringRateparticle-SctMat))<1e-11)
    end
end



%% Particle Class: ScatteringRatev2
mf = -9/2:1:9/2
for j = 1:10
    for i = mf
        TestLaser = LaserClass("TestZPositive")
        LaserArray = [TestLaser,TestLaser]
        particle.r =[0.001,0.12,-0.0132]*1e-3;
        particle.v = rand(1,3);
        particle.ZeemanSublevel = i;
        % 
        FormatObject = FormatInputsv2(LaserArray)
        SctMat = particle.ScatteringRatev2(FormatObject)
        assert(norm(abs(SctMat(1,:)-SctMat(2,:)))<1e-11)
    end
end

%% Particle Class: ScatteringRatev2 : Stirring laser
mf = -9/2:1:9/2
for j = 1:10
    for i = mf
        TestLaser = LaserClass("XCaptureStirringPositive")
        
        particle.r =[0.001,0.12,-0.0132]*1e-3;
        particle.v = rand(1,3);
        particle.ZeemanSublevel = i;
        % 
        FormatObject = FormatInputsv2(TestLaser)
        scatteringRateparticle = particle.ScatteringRate(TestLaser)
        SctMat = particle.ScatteringRatev2(FormatObject)
        assert(norm(abs(scatteringRateparticle-SctMat))<1e-11)
    end
end
%% Particle Class: ScatteringRatev2 : Stirring laser
mf = -9/2:1:9/2
for j = 1:10
    for i = mf
        TestLaser = LaserClass("XCaptureStirringPositive")
        TestLaser2 = LaserClass("XcapturePositive")

        particle.r =[0.001,0.12,-0.0132]*1e-3;
        particle.v = rand(1,3);
        particle.ZeemanSublevel = i;
        % 
        FormatObject = FormatInputsv2([TestLaser,TestLaser2])
        scatteringRateparticle = particle.ScatteringRate(TestLaser)
        scatteringRateparticle2 = particle.ScatteringRate(TestLaser2)
        exp = [scatteringRateparticle;scatteringRateparticle2]
        SctMat = particle.ScatteringRatev2(FormatObject)
        assert(norm(abs(exp-SctMat))<1e-11)
    end
end

% %% Particle Class: ScatteringRate
% % TestLaser = LaserClass("TestZPositive")
% % particle.r = [0,0,1];
% % particle.v = [0,1000,0];
% % s0 = TestLaser.s0;
% % expected = params.GammaRed*((s0/2)./(1+s0 + ()));
% % scatteringRateparticle = particle.ScatteringRate(TestLaser)
% % assert(abs(scatteringRateparticle(1)-expected)<1e-12)
% % assert(abs(scatteringRateparticle(2))<1e-12)
% % assert(abs(scatteringRateparticle(3))<1e-12)

% % % %% Particle Class: ScatteringRate
% % % TestLaser = LaserClass("TestZPositive")
% % % particle.r = [0,0,1];
% % % particle.v = [0,1000,0];
% % % s0 = TestLaser.s0;
% % % expected = params.GammaRed*((s0/2)./(1+s0));
% % % scatteringRateparticle = particle.ScatteringRate(TestLaser)
% % % assert(abs(scatteringRateparticle(1)-expected)<1e-12)
% % % assert(abs(scatteringRateparticle(2))<1e-12)
% % % assert(abs(scatteringRateparticle(3))<1e-12)