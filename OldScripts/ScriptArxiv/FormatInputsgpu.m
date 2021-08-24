classdef FormatInputsgpu
   properties
       %Lasers
       LaserArray
       LaserDetuningMatrix;
       LaserDetuningMatrixgpu;
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

       
      
       

   end
    methods
        function obj =  FormatInputsgpu(LaserArray)
            obj.LaserArray = LaserArray;
            obj.LaserDetuningMatrix = gpuArray(obj.MakeLaserDetuningMatrix(LaserArray));
            obj.BeamPropagationMatrix = (obj.getBeamPropagationMatrix(LaserArray));
            obj.polarizationArray = (obj.getPolarizationArray(LaserArray));
            obj.WaistArray = (obj.getWaistArray(LaserArray));
            obj.RayleightLegnthArray = ( obj.getRayleightLegnthArray(LaserArray));
            obj.waistDirectionArrayx = (obj.getwaistdirectionxArray(LaserArray));
            obj.waistDirectionArrayy = (obj.getwaistdirectionyArray(LaserArray));
            obj.offsetarray = (obj.getoffsetArray(LaserArray));
            obj.PeakIntesityArray = (obj.getPeakIntensityArray(LaserArray));
            obj.NumberoffCombLinesArray = (obj.getNumberofComblinesArray(LaserArray));
            obj.MagneticMomentAdresseTransitionArray = (obj.getTransitionAdressesMomentdArray(LaserArray));
            obj.adressedTransitionArray =( obj.getAdressedTransitionArray(LaserArray));
            obj.LaserDetuningMatrixgpu = ( obj.LaserDetuningMatrix);
            
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