
  function [PSC] = ScatteringRatev3(CGTable,cteArray,LaserInfoMat,AtomInfoMat)
      for jj = 1:100

    %Constants
    G =cteArray(9:11);
    kred =  cteArray(4);      
    CGtableTrapping =   CGTable(1:10,:);
    CGtableSitrring =  CGTable(11:20,:);
    Isat = cteArray(3);
    GammaRed =  cteArray(1);
    
    %Laser Info
    NumberOfcomblines =  LaserInfoMat(:,1);
    magneticMomentTransition = LaserInfoMat(:,2);
    BeamPropagation = LaserInfoMat(:,3:5);
    PeakIntensity =  LaserInfoMat(:,6);
    Offset = LaserInfoMat(:,7:9);
    DirWx = LaserInfoMat(:,10:12);
    DirWy = LaserInfoMat(:,13:15);
    Zr = LaserInfoMat(:,16:17);
    w = LaserInfoMat(:,18:19);
    PolarizationCR = LaserInfoMat(:,20:22);
    PolarizationCL = LaserInfoMat(:,23:25);
    AdressTransition = LaserInfoMat (:,26);
    DetuningMatrix = LaserInfoMat(:,27:end);
    
    %Atom info
    r =  AtomInfoMat(:,1:3);
    v =  AtomInfoMat(:,4:6);
    ZeemanSublevel =  AtomInfoMat(:,8);
    
    %ZeemanShifts 
    auxVec = [1,-1,0];
    MagneticField =G.*r;
    ZeemanDetunings =magneticMomentTransition *norm(MagneticField)*(ZeemanSublevel+auxVec);                


    %DopllerShift

    DopplerShiftMatrix = kred*(BeamPropagation*v');


    %Detunings
    Deltasp = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,1);
    Deltasm = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,2);
    Deltapi = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,3) ;

    Bdir = MagneticField/norm(MagneticField);
    LaserDir = BeamPropagation;
    costeta = LaserDir*Bdir';
    sintetasq = 1-costeta.^2;

    PolComponentssq = (PolarizationCR).*[(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2] + (PolarizationCL).*[((1-costeta).^2)/4,((1+costeta).^2)/4,sintetasq/2];
    %             

    %Get polarization of laser
%     PolComponentssq = GetLaserPolarizationComponentsVector(BeamPropagation,MagneticField,PolarizationCR,PolarizationCL);

    %Get average saturation parameter per combline
%     Intensity = getlaserIntesitiesv3(r,Offset,DirWx,DirWy,BeamPropagation,PeakIntensity, Zr,w);  
    
    %Get average intensity laser
    
    r = r-Offset;
    x = r.*DirWx;
    y = r.*DirWy;
    z = r.*BeamPropagation;          

    I0 =PeakIntensity;

    zrx =Zr(:,1);
    zry =Zr(:,2);
    
    wx = w(:,1);
    wy = w(:,2);


%             zsq = vecnorm (z,2,2).^2;
%             xsq = vecnorm (x,2,2).^2;
%             ysq = vecnorm (y,2,2).^2;
    zsq = sum(z.*z,2);
    xsq = sum(x.*x,2);
    ysq = sum(y.*y,2);

    omegax = wx.*sqrt(1+(zsq./(zrx.^2)) );
    omegay = wx.*sqrt(1+(zsq./(zry.^2)) );

    Intensity =  I0.*(wx./omegax).*(wy./omegay).*exp(-2*(xsq)./(omegax.^2)).*exp(-2*(ysq)./(omegay.^2));    
    
    
    s0 =  Intensity./NumberOfcomblines/Isat;

    %Get CG coeficients
    n = (ZeemanSublevel*2-1)/2+6;
    TransitionStrenght = CGtableTrapping(n,:).^2;
    TransitionStrenghtMAt = (AdressTransition == 1)* CGtableTrapping(n,:).^2  + (AdressTransition== 0)* CGtableSitrring(n,:).^2;
%     s0eff = s0*TransitionStrenght.*PolComponentssq;
    s0eff = bsxfun(@times,TransitionStrenghtMAt.*PolComponentssq,s0);


    %Calculate scattering rate.
    PSctsp = sum(GammaRed*(s0eff(:,1)/2)./(1+s0eff(:,1)+ (2*Deltasp/GammaRed).^2),2);
    PSctsm = sum(GammaRed*(s0eff(:,2)/2)./(1+s0eff(:,2)+ (2*Deltasm/GammaRed).^2),2);
    PSctpi = sum(GammaRed*(s0eff(:,3)/2)./(1+s0eff(:,3)+ (2*Deltapi/GammaRed).^2),2);
    PSC  = sum(sum([PSctsp,PSctsm,PSctpi]));
      end

  end
 