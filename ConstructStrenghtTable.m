function [StrenghtTableTrapping ,StrenghtTableStirring,StrenghtTableTrappingDecay,StrenghtTableStirringDecay] = ConstructStrenghtTable ()
%Table for excitation probabilities

Fp = [9/2,11/2];
F = [9/2];
counter2 = 0;
for stateExcited = Fp
    counter2 = counter2 +1;
    for stateGround = F
        counter = 0;
        for mfG= -stateGround:1:stateGround
            mfe = mfG +1;
            counter = counter +1;
            vecsp(counter,counter2) = ClebschGordan(stateGround,1,stateExcited,mfG,mfe-mfG,mfe);
           
            mfe = mfG -1;
            vecsm(counter,counter2) = ClebschGordan(stateGround,1,stateExcited,mfG,mfe-mfG,mfe);
            
            mfe = mfG;
            vecpi(counter,counter2) = ClebschGordan(stateGround,1,stateExcited,mfG,mfe-mfG,mfe);
        end
    end
end
StrenghtTableTrapping = [vecsp(:,2),vecsm(:,2),vecpi(:,2)];
StrenghtTableStirring = [vecsp(:,1),vecsm(:,1),vecpi(:,1)];
counter2 = 0;
for stateExcited = Fp
    counter2 = counter2 +1;
    for stateGround = F
        counter = 0;
        for mfe= -stateExcited:1:stateExcited
            mfG = mfe - 1;
            counter = counter +1;
            vecspd(counter,counter2) = ClebschGordan(stateGround,1,stateExcited,mfG,1,mfe);
           
            mfG = mfe + 1;
            vecsmd(counter,counter2) = ClebschGordan(stateGround,1,stateExcited,mfG,-1,mfe);
            
            mfG = mfe;
            vecpid(counter,counter2) = ClebschGordan(stateGround,1,stateExcited,mfG,0,mfe);
        end
    end
end

StrenghtTableTrappingDecay = [vecspd(:,2),vecsmd(:,2),vecpid(:,2)];
StrenghtTableStirringDecay = [vecspd(:,1),vecsmd(:,1),vecpid(:,1)];