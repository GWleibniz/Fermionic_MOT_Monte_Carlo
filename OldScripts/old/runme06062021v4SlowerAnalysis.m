path = "Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\" +"Data1_2v2.mat";
data = load(path)
stirringOfflost = 0;
stirringOnLost =0;
totalstirring = 0 ;
totlaStirringoff  = 0;
 paramsArray = data.loadedData;
for k = 1:length(paramsArray)
    try
    mfinitial = paramsArray(k).initialConditions(end)
    if mfinitial == -9/2
        
            if length(paramsArray(k).LaserArray) >6 
                stirringOnLost = stirringOnLost + (paramsArray(k).atomlostBolean  == "Atom Lost");
                totalstirring =  totalstirring+1;

            else
                 stirringOfflost = stirringOfflost + (paramsArray(k).atomlostBolean == "Atom Lost");
                 totlaStirringoff = totlaStirringoff + 1;
            end

        nd
    end
    catch
    end
end

disp('Sittiring on % of lost atoms')
display(num2str(100*stirringOnLost/totalstirring))



disp('Sittiring off % of lost atoms')
display(num2str(100*stirringOfflost/totlaStirringoff))


% % path = "C:\Users\PC de Rodrigo\Documents\MATLAB\save\Slower\run1\" +"Data0_5.mat";
% % data = load(path)
% % stirringOfflost = 0;
% % stirringOnLost =0;
% % totalstirring = 0 ;
% % totlaStirringoff  = 0;
% %  paramsArray = data.loadedData;
% % for k = 1:length(paramsArray)
% %     try
% %     mfinitial = paramsArray(k).initialConditions(end)
% %     if mfinitial == 3/2
% %             if (isa(paramsArray(k).atomlostBolean,'logical'))
% %             if length(paramsArray(k).LaserArray) >6 
% %                 stirringOnLost = stirringOnLost + (paramsArray(k).atomlostBolean);
% %                 totalstirring =  totalstirring+1;
% % 
% %             else
% %                  stirringOfflost = stirringOfflost + (paramsArray(k).atomlostBolean );
% %                  totlaStirringoff = totlaStirringoff + 1;
% %             end
% %             end
% %         
% %     end
% %     catch
% %     end
% % end
% % 
% % disp('Sittiring on % of lost atoms')
% % display(num2str(100*stirringOnLost/totalstirring))
% % 
% % 
% % 
% % disp('Sittiring off % of lost atoms')
% % display(num2str(100*stirringOfflost/totlaStirringoff))
