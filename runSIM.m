function runSIM(params,laserInfo,atomInfoMat,cteArray,CGTable,ExitMODE,run_number)
        SavePath = params.SavePath ;
        filename = num2str(run_number) ;
        path_fille = SavePath + "\" +filename + ".mat"; 
        if exist(path_fille, 'file' )== 2  
            display('file already exist');
            return            
        end
        initialConditions = atomInfoMat;
        i = 0;
        savePoints = 10;%Save data every 10 points;
        counter= 0;
        exitFlag = 1;
        atomInfoMat2Save = zeros(1e3,8);
        if (strcmp(ExitMODE,'MOTSim'))
            steps = params.StopTime/params.UpdateStepSize;
            h = waitbar(0,'Please wait...');
        end
        while exitFlag 
            if (rem(i,1000)==0) && (strcmp(ExitMODE,'MOTSim'))
                waitbar(i/steps);
            end
                i = i+1;        
               [DeltaV,atomicstate,atomZeemanSublevel] = getVelocityChangev2(laserInfo,atomInfoMat,cteArray,CGTable,params.BfieldFitPolinomial);


               DeltaV = DeltaV + [0,-params.Gravity,0]*params.UpdateStepSize  ;

               Deltar = atomInfoMat(4:6)*params.UpdateStepSize;

               atomInfoMat(1:3) = atomInfoMat(1:3) + Deltar;
               atomInfoMat(4:6) = atomInfoMat(4:6) + DeltaV;
               atomInfoMat(7) = atomicstate;
               atomInfoMat(8) = atomZeemanSublevel;
               
               
               if rem(i,savePoints)==0
                 counter =counter+1;
                 atomInfoMat2Save(counter,:) = atomInfoMat;
%                   atomInfoMat(1:3)*1000
               end
            runs =counter*savePoints;
            exitFlag = exitFUN(atomInfoMat,runs,ExitMODE,params);
        end
        if (strcmp(ExitMODE,'MOTSim'))
            close(h);
        end
        params.SimData = atomInfoMat2Save;
        params.atomlostBolean  = getExitReason(atomInfoMat,runs,ExitMODE,params);
        params.initialConditions = initialConditions;
        params.SaveParams("",filename);
        
end

function exitFlag = exitFUN(atomInfoMat,runs,ExitMODE,params)
    switch ExitMODE
        case "MOTSim"
            exitFlag =  ~((runs*params.UpdateStepSize > params.StopTime) || (atomInfoMat(2) < -100e-3));
        case "SlowerSIM"
            atomLostBoolean = (atomInfoMat(2) < -30e-3);
            atomSlowedBoolean = (norm(atomInfoMat(4:6))<0.1);
       
            exitFlag = ~(atomSlowedBoolean | atomLostBoolean);
        otherwise
            exitFlag = false;
    end


end

function exitReason = getExitReason(atomInfoMat,runs,ExitMODE,params)
    switch ExitMODE
        case "MOTSim"
            exitFlag =  (runs*params.UpdateStepSize > params.StopTime);
            atomLostBoolean = (atomInfoMat(2) < -100e-3);

        case "SlowerSIM"
            atomLostBoolean = (atomInfoMat(2) < -30e-3);
            atomSlowedBoolean = (norm(atomInfoMat(4:6))<0.1);
        otherwise
            exitFlag = false;
 
    end
    
    switch ExitMODE
        case "MOTSim"
            if (exitFlag)
                exitReason = "Sim Time finished"
            else
                exitReason = "Atom Lost"
            end
        case "SlowerSIM"
            if atomLostBoolean
                exitReason = "Atom Lost"
            elseif atomSlowedBoolean
                exitReason = "Atom Slowed"
            end
        otherwise
            exitReason = "exit mode undefined"
    end

end