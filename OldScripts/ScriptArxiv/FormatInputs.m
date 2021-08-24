classdef FormatInputs
   properties
       %Lasers
       LaserArray
       LaserDetuningMatrix;
       BeamPropagationMatrix;
       polarizationArray;
       
       WaistArray;
       waistDirectionArrayx;
       waistDirectionArrayy;
       offsetarray;
       RayleightLegnthArray;
       PeakIntesityArray;
       NumberoffCombLinesArray;
       
       MagneticMomentAdresseTransitionArray;
       adressedTransitionArray;
       
       
       %Atoms
       velocityArray;
       positionArray;
       atomicStateArray;
       ZeemSublevelarray;
       atomicStateexcited;
       
      
       

   end
    methods
        function obj =  FormatInputs(LaserArray,atomarray)
            
            
            %MOT lasers ( Each row is the parameter of the laser)
            obj.LaserArray = LaserArray;
            obj.LaserDetuningMatrix = obj.MakeLaserDetuningMatrix(LaserArray);
            obj.BeamPropagationMatrix = obj.getBeamPropagationMatrix(LaserArray);
            obj.polarizationArray = obj.getPolarizationArray(LaserArray);
            obj.WaistArray = obj.getWaistArray(LaserArray);
            obj.RayleightLegnthArray =  obj.getRayleightLegnthArray(LaserArray)
            obj.waistDirectionArrayx = obj.getwaistdirectionxArray(LaserArray)
            obj.waistDirectionArrayy = obj.getwaistdirectionyArray(LaserArray)
            obj.offsetarray = obj.getoffsetArray(LaserArray);
            obj.PeakIntesityArray = obj.getPeakIntensityArray(LaserArray);
            obj.NumberoffCombLinesArray = obj.getNumberofComblinesArray(LaserArray);
            obj.MagneticMomentAdresseTransitionArray = obj.getTransitionAdressesMomentdArray(LaserArray);
            obj.adressedTransitionArray = obj.getAdressedTransitionArray(LaserArray);            
            
            if (nargin >1)
            obj.velocityArray = obj.getVelocityarray(atomarray);
            obj.positionArray = obj.getpositionArray(atomarray);
            obj.atomicStateArray = obj.getatomicStateArray(atomarray);
            obj.atomicStateexcited = obj.getatomicStateexcited(atomarray);
            obj.ZeemSublevelarray = obj.getZeemSublevelarray(atomarray);
            end
        end
        function ZeemSublevelarray = getZeemSublevelarray(obj,atomarray)
            for i = 1:length(atomarray)
                ZeemSublevelarray(i,:) = atomarray(i).ZeemanSublevel;   
            end
            
        end
        
        
        function atomicStateArray = getatomicStateArray(obj,atomarray)
            for i = 1:length(atomarray)
                atomicStateArray(i,:) = atomarray(i).AtomicState;   
            end
            
        end
        function atomicStateexcited = getatomicStateexcited(obj,atomarray)
            for i = 1:length(atomarray)
                atomicStateexcited(i,:) = atomarray(i).AtomicStateExcited;   
            end
            
        end

        function velocity = getVelocityarray(obj,atomarray)
            for i = 1:length(atomarray)
                velocity(i,:) = atomarray(i).v;   
            end
            velocity = velocity';
            velocity = velocity(:);
        end
        function positionArray = getpositionArray(obj,atomarray)
            for i = 1:length(atomarray)
                positionArray(i,:) = atomarray(i).r;   
            end
            positionArray = positionArray';
        
            positionArray = positionArray(:);
        end
        
        
        
        
        function LaserDetuningMatrix = MakeLaserDetuningMatrix(obj,LaserArray)
            for i = 1:length(LaserArray)
                combLines = LaserArray(i).CombLines;   
                LaserDetuningMatrix(i,1:length(combLines)) = LaserArray(i).CombLines;
%                 LaserDetuningMatrix(LaserDetuningMatrix == 0) = ;
            end
            sizeDetMAt = size(LaserDetuningMatrix,2)
            
            for i = 1:length(LaserArray)
                numberCombLines = length(LaserArray(i).CombLines);
                if (sizeDetMAt > numberCombLines)
                LaserDetuningMatrix(i,numberCombLines:end) =Inf;
%                 LaserDetuningMatrix(LaserDetuningMatrix == 0) = ;
                end
            end
        end
        function BeamPropagationMatrix = getBeamPropagationMatrix(obj,LaserArray)
            for i = 1:length(LaserArray)
                BeamPropagationMatrix(i,:) = LaserArray(i).BeamPropagationDirection;
            end
        end
        function PolArray = getPolarizationArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                PolArray(i,:) = LaserArray(i).Polarization;
            end
            PolArray = repmat(PolArray,1,3)
        end
        function WaistArray = getWaistArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                WaistArray(i,1) = LaserArray(i).wx;
                WaistArray(i,2) = LaserArray(i).wy;
            end
        end
        function RayleightLegnthArray = getRayleightLegnthArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                RayleightLegnthArray(i,1) = LaserArray(i).RaleightLenghtx;
                RayleightLegnthArray(i,2) = LaserArray(i).RaleightLenghty;
            end
        end
        function directionwxArray = getwaistdirectionxArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                directionwxArray(i,:) = LaserArray(i).directionwx;
            end
        end
        function directionwyArray = getwaistdirectionyArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                directionwyArray(i,:) = LaserArray(i).directionwy;
            end
        end
        function offsetArray = getoffsetArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                offsetArray(i,:) = LaserArray(i).offset;
            end
        end
        function peakIntensityarray = getPeakIntensityArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                peakIntensityarray(i,:) = LaserArray(i).peakIntensity;
            end
        end
        function numberofcomblinesArray = getNumberofComblinesArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                numberofcomblinesArray(i,:) = length(LaserArray(i).CombLines);
            end
        end
        function MagneticMomentAdresseTransitionArray = getTransitionAdressesMomentdArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                MagneticMomentAdresseTransitionArray(i,:) = LaserArray(i).MagneticMomentAdresseTransition;
            end
        end
        function adressedTransition = getAdressedTransitionArray(obj,LaserArray)
            for i = 1:length(LaserArray)
                adressedTransition(i,:) = LaserArray(i).IsTrappingLaser;
            end
        end
    end
end