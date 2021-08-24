classdef LaserClass
   properties
   %
    LaserName;
       
   %Laser power in W
    Power = 10.8*10^-3; 
    peakIntensity;
    s0; 
    
    %Beam difrection
    BeamPropagationDirection = [0,1,0];
    offset = [0,0,0]
    directionwx =[1,0,0] ;
    directionwy = [0,0,1];
    
    %Polarization
    Polarization = "CL";%x pm iy sigma +- , x lineal 
    
    %Waist in m
    wx = 34*10^-3 ;
    wy = 34*10^-3; 
    RaleightLenghtx;
    RaleightLenghty;
    
    %Stirring or Trapping
    IsTrappingLaser = 1;
    MagneticMomentAdresseTransition;
    %detunings
    DeltaStart  = -2*pi*0.95*10^6;
    DeltaStop = -2*pi*5*10^6;
    modulation = 2*pi*0.015*10^6;
    
    CombLines ;
    powerPerCombline;
    IntensityPercombLine;
    IsatPerCombLine;
   end
   methods
       
      function obj = LaserClass(LaserName)
          params = ParamClass;
          
          obj = obj.setLaserParameters(LaserName);
          
          obj.peakIntensity = 2*obj.Power/pi/obj.wx/obj.wy;
          obj.s0 = obj.peakIntensity/params.Isat;
          
          obj.CombLines = obj.DeltaStop:obj.modulation:obj.DeltaStart;
          obj.powerPerCombline = obj.Power/length(obj.CombLines);
          obj.IntensityPercombLine = 2*obj.powerPerCombline/pi/obj.wx/obj.wy;
          obj.IsatPerCombLine =  obj.IntensityPercombLine/params.Isat;
          
          obj.RaleightLenghtx = params.kred*obj.wx^2/2;
          obj.RaleightLenghty = params.kred*obj.wy^2/2;
          
          obj.MagneticMomentAdresseTransition = getMomentAdressTransition(obj,params);
      end
      function polarizationComponents = getPolarizationComponents(obj,B)
          %Gives the polarization decomposition of the laser as 
          %[sp,sm,pi] = polarizationComponents
          if (B == zeros(1,3))
              polarizationComponents = [1,1,1]/sqrt(3);
          else
              Bdir = B/norm(B);
              LaserDir = obj.BeamPropagationDirection;
              costeta = Bdir*LaserDir';
              sinteta = norm(cross(Bdir,LaserDir));


              if (obj.Polarization =="CR")
                polarizationComponents = [(1+costeta)/2,(1-costeta)/2,sinteta/sqrt(2)];
              elseif(obj.Polarization =="CL")
                polarizationComponents = [(1-costeta)/2,(1+costeta)/2,sinteta/sqrt(2)];
              end
          end
          
      end
      function MagneticMomentAdresseTransition = getMomentAdressTransition(obj,params)
          switch  obj.IsTrappingLaser 
              case 1
                  moment = params.MagneticMoment3P1/5.5;
              case 0
                  moment = params.MagneticMoment3P1/4.5/4.5;
              otherwise
                  error('LaserClass : getMomentAdressTransition() , unknow transition')
          end
          
         MagneticMomentAdresseTransition = moment;
      end
    function polarizationComponentssq = getPolarizationComponentssq(obj,B)
          %Gives the polarization decomposition of the laser as 
          %[sp,sm,pi] = polarizationComponents
          if (B == zeros(1,3))
              polarizationComponentssq = [1,1,1]/3;
          else
              Bdir = B/norm(B);
              LaserDir = obj.BeamPropagationDirection;
              costeta = Bdir*LaserDir';
              sintetasq = 1-costeta.^2;


              if (obj.Polarization =="CR")
                polarizationComponentssq = [(1+costeta).^2/4,(1-costeta).^2/4,sintetasq/2];
              elseif(obj.Polarization =="CL")
                polarizationComponentssq =  [(1-costeta).^2/4,(1+costeta).^2/4,sintetasq/2];
              end
          end
          
      end
      function intensityBeam = LaserIntensityProfile(obj,r)   
        
        if obj.LaserName == "Ycapture"
            intensityBeam = yCaptureIntensityProfile(obj,r);  
        else
            r = r - obj.offset;
            x = r.*obj.directionwx;
            y = r.*obj.directionwy;
            z = r.*obj.BeamPropagationDirection;          

            I0 = obj.peakIntensity;
            zrx = obj.RaleightLenghtx;
            zry = obj.RaleightLenghty;
            omegax = obj.wx*sqrt(1+(z*z'/(zrx)^2));
            omegay = obj.wy*sqrt(1+(z*z'/(zry)^2));
            intensityBeam =  I0*(obj.wx/omegax)*(obj.wy/omegay)*exp(-2*(x*x')/omegax^2)*exp(-2*(y*y')/omegay^2);
        end
        
      end
    function intensityBeam = yCaptureIntensityProfile(obj,r)   
        
          
        r = r - obj.offset;
        x = r.*obj.directionwx;
        y = r.*obj.directionwy;
        z = r.*obj.BeamPropagationDirection;          
         
        I0 = obj.peakIntensity;
        zrx = obj.RaleightLenghtx;
        zry = obj.RaleightLenghty;
        omegax = ((1-z(2)/0.25))*obj.wx*sqrt(1+(z*z'/(zrx)^2));
        omegay = ((1-z(2)/0.25))*obj.wy*sqrt(1+(z*z'/(zry)^2));
        intensityBeam =  I0*(obj.wx/omegax)*(obj.wy/omegay)*exp(-2*(x*x')/omegax^2)*exp(-2*(y*y')/omegay^2);
        
        
      end
      
      
      

      function obj = setLaserParameters(obj,LaserName)
          obj.LaserName = LaserName;
          switch LaserName 
              case "Ycapture"          
                %Laser power in W
                obj.Power = 121.2*10^-3; 
                %Beam difrection
                obj.BeamPropagationDirection = [0,1,0];
                obj.directionwx =[1,0,0] ;
                obj.directionwy = [0,0,1];
                obj.offset = [0,440,0]*10^-3;
                
                %Waist in m (beam is focused 44 mm above quadrupole center)
                obj.wx = 2.8e-6 ;
                obj.wy = 2.8e-6; 
%                 obj.wx = 34*10^-3 ;
%                 obj.wy = 34*10^-3; 
%                 %detunings
                obj.DeltaStart  = -2*pi*0.95*10^6;
                obj.DeltaStop = -2*pi*3.5*10^6;
                obj.modulation = 2*pi*0.0150*10^6;
                obj.Polarization = "CR";%x pm iy sigma +- , x lineal 

              case "XcapturePositive"
                %Laser power in W
                obj.Power = 0.66*10^-3; 
                %Beam difrection
                obj.BeamPropagationDirection = [1,0,0];
                obj.directionwx =[0,0,1] ;
                obj.directionwy = [0,1,0];
                %Waist in m
                obj.wx = 23.5*10^-3 ;
                obj.wy = 23.5*10^-3; 
                obj.offset = [0,-13,0]*10^-3;
                %detunings
                obj.DeltaStart  = -2*pi*0.65*10^6;
                obj.DeltaStop = -2*pi*2.19*10^6;
                obj.modulation = 2*pi*0.015*10^6 ;               
                obj.Polarization = "CL";%x pm iy sigma +- , x lineal 
                
                obj.IsTrappingLaser = 1;

              case "XcaptureNegative"
                obj = obj.setLaserParameters("XcapturePositive");
                obj.LaserName = "XcaptureNegative";
                obj.BeamPropagationDirection = [-1,0,0];
              case "ZcapturePositive"

                  %Laser power in W
                obj.Power =0.1*10^-3; 
                %Beam difrection
                obj.BeamPropagationDirection = [0,0,1];
                obj.directionwx =[1,0,0] ;
                obj.directionwy = [0,1,0];
                %Waist in m
                obj.wx = 4*10^-3 ;
                obj.wy = 4*10^-3; 
%                 obj.offset = [0,-13,0]*10^-3;
                obj.offset = [0,-10,0]*10^-3;
                %detunings
                obj.DeltaStart  = -2*pi*0.82*10^6;
                obj.DeltaStop = -2*pi*1.24*10^6;
                obj.modulation = 2*pi*0.016*10^6;    
                obj.Polarization = "CR";%x pm iy sigma +- , x lineal 
                obj.IsTrappingLaser = 1;
              case "ZcaptureNegative"
                  obj = obj.setLaserParameters("ZcapturePositive");
                  obj.BeamPropagationDirection =[0,0,-1];
                  obj.LaserName = "ZcaptureNegative";
              case "YcaptureStirring"
                  obj = obj.setLaserParameters("Ycapture");
                  obj.LaserName = "YcaptureStirring";
                  obj.Power = 6.6*10^-3;
                  obj.IsTrappingLaser = 0;
              case "XCaptureStirringPositive"
                  obj = obj.setLaserParameters("XcapturePositive");
                  obj.LaserName = "XCaptureStirringPositive";
                  obj.Power = 0.3*10^-3;
                  obj.IsTrappingLaser = 0;  
                case "XCaptureStirringNegative"
                  obj = obj.setLaserParameters("XcaptureNegative");
                  obj.LaserName = "XCaptureStirringNegative";
                  obj.Power = 0.3*10^-3;
                  obj.IsTrappingLaser = 0;
              case "TestZPositive"
                  

                    %Laser power in W
                    obj.Power =0.1*10^-3; 
                    %Beam difrection
                    obj.BeamPropagationDirection = [0,0,1];
                    obj.directionwx =[1,0,0] ;
                    obj.directionwy = [0,1,0];
                    %Waist in m
                    obj.wx = 8*10^-3 ;
                    obj.wy = 8*10^-3; 
                    obj.offset = [0,0,0]*10^-3;

                    %detunings
                    obj.DeltaStart  = 6;
                    obj.DeltaStop = 0;
                    obj.modulation = 10;  
                    obj.Polarization = "CR";%x pm iy sigma +- , x lineal 
                    obj.IsTrappingLaser = 1;
              otherwise
                  error("obj.LaserName does not exist")
            end  
            
            
            
      end
   end
end