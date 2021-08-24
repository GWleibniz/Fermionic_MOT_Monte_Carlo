classdef PhysicsClass
   properties
    %fundamental Constants
     m = 88*(1.66e-27); %Kg
     
    %Red transition properties
    GammaRed = 2*pi*7.5*10^3; %Hertz
    MagneticMoment3P1 = 2*pi*2.1*10^6; %2.1 MHz/Gauss 
    Isat = 3*10^-2; %W /m^2
    kred = 2*pi/(689e-9);
    hbar = 1.054571817e-34;
    %recoil Velocity 
    vrecoil;
    
    Gravity = 9.81;
    
    %MOT magnetic gradient
    BGradient = [-0.55,0.32,0.23]*10^(2); %G/m
    
    %Lasers
    LaserArray;
   end
   methods
       
      function obj = PhysicsClass()
          obj.vrecoil = obj.hbar*obj.kred/obj.m;
      end
        
      function Pe = Pe(obj,Delta,s0)
            Pe = (s0/2)/(1+s0+ (2*Delta/obj.GammaRed)^2);
            
   end
   end
end