classdef ParticleClass
   properties
      
      %Initial location and velocity
      r0 = [0,0,0]*10^-3;
      v0 = [0,0,0]*10^-3 ;
      
      %Instantaneous velocity
      r = [];
      v = [];
      
      
      %Trayectory Matrix
      TrayectoryMat=[];
      IntenralStateMAt =[];
      
      %Atomic AtomicState ( 0 = ground AtomicState 1 = excited state)
      AtomicState = 0;
      % value of excited state hyperfine manifold. I dont label ground state 
      AtomicStateExcited = 9/2; 
      ZeemanSublevel = 9/2;
      
      %properties
      m = 87*1.66e-27;
      
    
      %other
      UpdateStepSize = 10e-6; %s
      
      params ;
   end
    methods
        function obj = ParticleClass(params,r0,v0,mf)
            if nargin > 1
                obj.r0 = r0;
                obj.v0 = v0;
                obj.ZeemanSublevel = mf;
            end
            obj.params = params;
            obj.r = obj.r0;
            obj.v = obj.v0;
        end

        function obj = UpdateLocation(obj)

            h =obj.params.UpdateStepSize;
            obj.r = obj.r + obj.v*h;
        end

        function obj = UpdateVelocity (obj,deltaV)
            obj.v = obj.v + deltaV;
        end
        
        function obj =storeInternalStateData(obj,state,mf,Fp,time)
            row  = [state,mf,Fp,time];
            obj.IntenralStateMAt(end+1,:) =row;
        end
        
        function obj =storeTrayectoryData(obj,r,v,time)
            row  = [r,v,time];
            obj.TrayectoryMat(end+1,:) =row;
        end
        
       function obj = SaveParticle(obj)
             SavePath = obj.params.SavePath ;
             filename = SavePath + "\runme1"  + "\atom.mat";
             atomnumber = 0;
             while isfile(filename)
                 atomnumber =atomnumber +1;
                 filename = SavePath + "\atom" + num2str(atomnumber) + ".mat";

             end
             save(filename,'obj')

        end

      
        function [Psc] = ScatteringRatev2(obj,MOTlasers)
            
                G = obj.params.BGradient;
                kred =   obj.params.kred;      
                MagneticMoment3P1 =obj.params.MagneticMoment3P1;
                CGtableTrapping =   obj.params.CGTable.TrappingTransition;
                CGtableSitrring =   obj.params.CGTable.StirringTransition;

                Isat = obj.params.Isat;
                GammaRed = obj.params.GammaRed;
                MagneticMomentArray = MOTlasers.MagneticMomentAdresseTransitionArray;
                
                DetuningMatrix = MOTlasers.LaserDetuningMatrix;
                BeamPropagationDirectionMatrix = MOTlasers.BeamPropagationMatrix;
                PolarizationArray = MOTlasers.polarizationArray;
                
                %ZeemanShifts 
                auxVec = [1,-1,0];
                MagneticField =G.*obj.r;
                ZeemanDetunings =MagneticMomentArray *norm(MagneticField)*(obj.ZeemanSublevel+auxVec);                
                
                
                %DopllerShift
                DopplerShiftMatrix = kred*(BeamPropagationDirectionMatrix*obj.v');

                
                %Detunings
                Deltasp = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,1);
                Deltasm = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,2);
                Deltapi = DetuningMatrix - DopplerShiftMatrix - ZeemanDetunings(:,3) ;
                
                %Get polarization of laser
                PolComponentssq = obj.GetLaserPolarizationComponentsVector(BeamPropagationDirectionMatrix,MagneticField,PolarizationArray);

                %Get average saturation parameter per combline
                Intensity =   obj.getlaserIntesities(obj.r,MOTlasers);
                s0 =  Intensity./MOTlasers.NumberoffCombLinesArray/Isat;
                
                %Get CG coeficients
                n = (obj.ZeemanSublevel*2-1)/2+6;
                TransitionStrenght = CGtableTrapping(n,:).^2;
                TransitionStrenghtMAt = (MOTlasers.adressedTransitionArray == 1)* CGtableTrapping(n,:).^2  + (MOTlasers.adressedTransitionArray == 0)* CGtableSitrring(n,:).^2;
                s0eff = s0*TransitionStrenght.*PolComponentssq;
                s0eff = bsxfun(@times,TransitionStrenghtMAt.*PolComponentssq,s0);

                
                %Calculate scattering rate.
                PSctsp = sum(GammaRed*(s0eff(:,1)/2)./(1+s0eff(:,1)+ (2*Deltasp/GammaRed).^2),2);
                PSctsm = sum(GammaRed*(s0eff(:,2)/2)./(1+s0eff(:,2)+ (2*Deltasm/GammaRed).^2),2);
                PSctpi = sum(GammaRed*(s0eff(:,3)/2)./(1+s0eff(:,3)+ (2*Deltapi/GammaRed).^2),2);
                Psc = [PSctsp,PSctsm,PSctpi];
                
        end
        
      function intensityBeamMatrix = getlaserIntesities(obj,r,MOTLasers)   
        
            r = r-MOTLasers.offsetarray;
            x = r.*MOTLasers.waistDirectionArrayx;
            y = r.*MOTLasers.waistDirectionArrayy;
            z = r.*MOTLasers.BeamPropagationMatrix;          

            I0 = MOTLasers.PeakIntesityArray;
            
            zrx = MOTLasers.RayleightLegnthArray(:,1);
            zry = MOTLasers.RayleightLegnthArray(:,2);
            
            wx = MOTLasers.WaistArray(:,1);
            wy = MOTLasers.WaistArray(:,2);
            
            
%             zsq = vecnorm (z,2,2).^2;
%             xsq = vecnorm (x,2,2).^2;
%             ysq = vecnorm (y,2,2).^2;
            zsq = sum(z.*z,2);
            xsq = sum(x.*x,2);
            ysq = sum(y.*y,2);

            omegax = wx.*sqrt(1+(zsq./(zrx.^2)) );
            omegay = wx.*sqrt(1+(zsq./(zry.^2)) );
            
            intensityBeamMatrix =  I0.*(wx./omegax).*(wy./omegay).*exp(-2*(xsq)./(omegax.^2)).*exp(-2*(ysq)./(omegay.^2));
        
        
      end   
        function intensityBeamMatrix = getlaserIntesitiesv2(obj,r,MOTLasers)   
          
            r = r-MOTLasers.offsetarray;
            x = r.*MOTLasers.waistDirectionArrayx;
            y = r.*MOTLasers.waistDirectionArrayy;
            z = r.*MOTLasers.BeamPropagationMatrix;          

            I0 = MOTLasers.PeakIntesityArray;
            
            zrx = MOTLasers.RayleightLegnthArray(:,1);
            zry = MOTLasers.RayleightLegnthArray(:,2);
            
            wx = MOTLasers.WaistArray(:,1);
            wy = MOTLasers.WaistArray(:,2);
            
            
%             zsq = vecnorm (z,2,2).^2;
%             xsq = vecnorm (x,2,2).^2;
%             ysq = vecnorm (y,2,2).^2;
            zsq = sum(z.*z,2);
            xsq = sum(x.*x,2);
            ysq = sum(y.*y,2);

            omegax = wx.*sqrt(1+(zsq./(zrx.^2)) );
            omegay = wx.*sqrt(1+(zsq./(zry.^2)) );
            
            intensityBeamMatrix =  I0.*(wx./omegax).*(wy./omegay).*exp(-2*(xsq)./(omegax.^2)).*exp(-2*(ysq)./(omegay.^2));
        
        
      end  
        function PolComponents = GetLaserPolarizationComponentsVector(obj,BeamPropagationDirectionMatrix,Bfield,PolarizationArray)
            
              Bdir = Bfield/norm(Bfield);
              LaserDir = BeamPropagationDirectionMatrix;
              costeta = LaserDir*Bdir';
              sintetasq = 1-costeta.^2;
              
              PolComponents = (PolarizationArray=="CR").*[(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2] + (PolarizationArray=="CL").*[((1-costeta).^2)/4,((1+costeta).^2)/4,sintetasq/2];
%             

        end
        
      
      function SctRate = ScatteringRate(obj,Laser)
          
        Params = obj.params;
        
        if (Laser.IsTrappingLaser)
            %Calculate Zeeman shift
            G = obj.params.BGradient;
            MagneticMoment3P1 =obj.params.MagneticMoment3P1/(5.5);
            auxVec = [1,-1,0];
            MagneticField =G.*obj.r;
            ZeemanDetunings = MagneticMoment3P1*norm(MagneticField)*(obj.ZeemanSublevel+auxVec);
            
            CGtable =   Params.CGTable.TrappingTransition;
        else
            %Calculate Zeeman shift
            G = obj.params.BGradient;
            MagneticMoment3P1 = Params.MagneticMoment3P1/4.5/(4.5);
            auxVec = [1,-1,0];
            MagneticField =G.*obj.r;
            ZeemanDetunings = MagneticMoment3P1*norm(MagneticField)*(obj.ZeemanSublevel+auxVec);
            
            CGtable =   Params.CGTable.StirringTransition;
        end
     
        kred =   Params.kred;
        vpr = Laser.BeamPropagationDirection*obj.v';
        DopplerShift = kred*vpr;
        %Get detunings of sp sm and sp transitions
%         DopplerShift = obj.getDopplerShift(Laser);
        Deltasp =Laser.CombLines - DopplerShift - ZeemanDetunings(1);
        Deltasm =Laser.CombLines - DopplerShift - ZeemanDetunings(2);
        Deltapi =Laser.CombLines - DopplerShift - ZeemanDetunings(3) ;
        
        %Get polarization of laser
        pol  = Laser.getPolarizationComponentssq(MagneticField);
        
        %Get CG coeficients
        
        n = (obj.ZeemanSublevel*2-1)/2+6;
        TransitionStrenght = CGtable(n,:).^2;
        s0 =  Laser.LaserIntensityProfile(obj.r)/length(Laser.CombLines)/Params.Isat;
        s0 = s0*pol.*TransitionStrenght;
        
        %Calculate scattering rate.
        PSctsp = sum((s0(1)/2)./(1+s0(1)+ (2*Deltasp/Params.GammaRed).^2));
        PSctsm = sum((s0(2)/2)./(1+s0(2)+ (2*Deltasm/Params.GammaRed).^2));
        PSctpi = sum((s0(3)/2)./(1+s0(3)+ (2*Deltapi/Params.GammaRed).^2));

        
        Psct =[PSctsp,PSctsm,PSctpi];
        if (max(Psct) >0.5)
            Psct(Psct>0.5) =  0.5;
%             warning('Excited AtomicState probability bigger thant 1/2')
        end
        SctRate = Params.GammaRed*Psct;
      end
      function Zs  = getZeemanShiftBosons (obj)
        values = obj.params;
        G = values.BGradient;
        
        MagneticMoment3P1 =values.MagneticMoment3P1;
        MagneticField =G.*obj.r;
        Zs = -MagneticMoment3P1*norm(MagneticField);
      end
    function Zs  = getZeemanShiftFermionsTrapping (obj)
        G = obj.params.BGradient;
        MagneticMoment3P1 =obj.params.MagneticMoment3P1/(5.5);
        auxVec = [1,-1,0];
        MagneticField =G.*obj.r;
        Zs = -MagneticMoment3P1*norm(MagneticField)*(obj.ZeemanSublevel+auxVec);
        
        
        
        %Zeeman shift of the [sp, sm pi] transitions
    end
    function Zs  = getZeemanShiftFermionsSittirring(obj)
        G = obj.params.BGradient;
        MagneticMoment3P1 =obj.params.MagneticMoment3P1/4.5/(4.5);
        auxVec = [1,-1,0];
        MagneticField =G.*obj.r;
        Zs = -MagneticMoment3P1*norm(MagneticField)*(obj.ZeemanSublevel+auxVec);
    end
      
     function deltaD  = getDoppler2Shift (obj,Laser)
          values = ParamClass;
          kred =   values.kred;
          
          v = Laser.BeamPropagationDirection*obj.v';
          deltaD = kred*v;
      end
   end
end