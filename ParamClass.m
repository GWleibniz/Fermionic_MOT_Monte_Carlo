classdef ParamClass
   properties
    %fundamental Constants
    kB =  1.3806503e-23;
    hbar = 1.054571817e-34;
     
    
    %Red transition properties
    GammaRed = 2*pi*7.5*10^3; %Hertz
    MagneticMoment3P1 = 2*pi*2.1*10^6; %2.1 MHz/Gauss 
    Isat = 3*10^-2; %W /m^2
    kred = 2*pi/(689e-9);
    vrecoil;
    Pdecay;
    
    
    cteArray;
    
    %ClebsGordan Coefficient table
    CGTable
    CGTableCombined;
    
    %Other constant
    Gravity = 9.81;
    m = 88*(1.66e-27); %88Sr mass

    %MOT magnetic gradient  
    BGradient = [-0.55,0.32,0.23]*10^(2); %G/m
    BfieldFitPolinomial;
        
    %Lasers
    LaserArray;
    atomArray;
    
    %
    StopTime;
    %SimData
    SimData;
    atomlostBolean;
    initialConditions;
    
    
    %save folder
    SavePath = 'Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\';
    %Update Step Size
    UpdateStepSize = 10e-6; %(s)
   end
   methods
       
      function obj = ParamClass()
          obj.vrecoil = obj.hbar*obj.kred/obj.m;
          obj.Pdecay = obj.UpdateStepSize*obj.GammaRed;
          obj = getCBTable(obj);
          obj.cteArray = obj.getCteArray();
%           obj.BfieldFitPolinomial = obj.getPol();
      end
   function pol = getPol(obj)
       path = [obj.SavePath , 'PolinomialFitMOTfield.mat'];
       var = load(path);
       pol = var.p;
       
      end
     function cteArray = getCteArray(obj)
            cteArray = [obj.GammaRed,obj.MagneticMoment3P1,obj.Isat,obj.kred,obj.hbar,obj.vrecoil,obj.UpdateStepSize,obj.Pdecay,obj.BGradient];
      end

      function obj = getCBTable(obj)
      [StrenghtTableTrapping ,StrenghtTableStirring,StrenghtTableTrappingDecay,StrenghtTableStirringDecay] = ConstructStrenghtTable ();
      obj.CGTable.TrappingTransition =StrenghtTableTrapping;
      obj.CGTable.StirringTransition =StrenghtTableStirring;
      obj.CGTable.TrappingTransitionDecay = StrenghtTableTrappingDecay;
      obj.CGTable.StirringTransitionDecay = StrenghtTableStirringDecay;
      obj.CGTableCombined = [obj.CGTable.TrappingTransition;obj.CGTable.StirringTransition;obj.CGTable.TrappingTransitionDecay;obj.CGTable.StirringTransitionDecay];



      end
      
      function SaveParams(obj,folder,Filename)
          try
              SavePath = obj.SavePath ;
              filename = SavePath + "\" + folder +  "\" + Filename + ".mat";             
              save(filename,'obj')
          catch
              display('unable to save file')
          end
      end
      
    function obj  = LoadParams(obj,path)

            try 
                data = load(path);
                var  = data.obj;
            catch
                var = obj;
                warning('file not found');
            end
            obj = var;
    end
          
   end
end