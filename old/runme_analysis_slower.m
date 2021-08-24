% % % % Data_Slower = DataClass("Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\","Data_1.mat")
% % % % 
% % % % mf_arr = -9/2:1:9/2
% % % % for mf = mf_arr
% % % %     display(['For mf = ',num2str(mf) ])
% % % %     Data_Slower.GetNumberOfTrappedTrayectories(mf);
% % % % end
% % % % 
% % % % 
% % % % mf_arr = -9/2:1:9/2
% % % % for mf = mf_arr
% % % %     figure
% % % %     display(['For mf = ',num2str(mf) ])   
% % % %     Data_Slower.do_PSD_plot_Average_all(2,mf);
% % % % end


Data_Slower = DataClass("Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\r9\","Data_1.mat")
Data_Slower.Plot_MOT_Capture_efficiency(9/2)
for mf = -9/2:-3/2

    Initial_mf = zeros (size(Data_Slower.paramsArray));
    Initial_v = zeros (size(Data_Slower.paramsArray));
    TrappedBool = zeros (size(Data_Slower.paramsArray));
    stirringOnBool = zeros (size(Data_Slower.paramsArray));

    for indx = 1:length(Data_Slower.paramsArray)
        param = Data_Slower.paramsArray(indx);
        Initial_v(indx) = param.initialConditions(5);
        Initial_mf(indx) = param.initialConditions(8);
        TrappedBool(indx) =  (param.atomlostBolean == "Atom Slowed");
        stirringOnBool(indx) = length(param.LaserArray) >5;
    end


    xvalues = unique(Initial_v)
    yvalues_stirring= zeros(size(xvalues))
    yvalues_No_sitrring = zeros(size(xvalues))

    for indx = 1:length(xvalues)
        valueBool_arr = (Initial_v==xvalues(indx) & Initial_mf == mf);
        number_stirring_on = sum(valueBool_arr & stirringOnBool )
        number_stirring_off = sum(valueBool_arr & (~stirringOnBool ))

        if (number_stirring_on >30 && number_stirring_off>30)
            yvalues_stirring(indx) = sum(TrappedBool(valueBool_arr & stirringOnBool ))/  number_stirring_on;
            yvalues_No_sitrring(indx) = sum(TrappedBool(valueBool_arr & (~stirringOnBool )))/  number_stirring_off;
        end
    end
    figure
    hold on
    plot(xvalues,yvalues_No_sitrring)
    plot(xvalues,yvalues_stirring)
    legend('Stirring off','Stirring on')
    title(['mf = ', num2str(mf)])

end
Data_Slower2 = DataClass("Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\r5\","Data_1.mat")
mf_arr = -9/2:1:9/2
for mf = mf_arr
    display(['For mf = ',num2str(mf) ])
    Data_Slower2.GetNumberOfTrappedTrayectories(mf);
end
% % % mf_arr = -9/2:1:9/2
% % % for mf = mf_arr
% % %     figure
% % %     display(['For mf = ',num2str(mf) ])   
% % %     Data_Slower.do_PSD_plot_Average_all(2,mf);
% % % end


