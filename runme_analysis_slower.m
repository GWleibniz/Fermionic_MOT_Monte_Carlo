%%% Loading Data

% folder = "Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\r10\BatchData\"
% files = dir(folder+"*.mat");
% for k = 1:length(files)
%     path(k) =  folder + files(k).name;
% end
% for k = 1:length(files)
%     k
%     Data_Slower(k) = DataClass(folder,string(files(k).name))
% end
% close all

%%% Ploting average capture velocity
mf_vec = -9/2:9/2
for indx = 1:length(mf_vec)
    mf = mf_vec(indx)
    yvalues_No_sitrring_Mat = [];
    yvalues_stirring_Mat = [];
    xvalues_Mat = [];
    for k = 1:length(Data_Slower)
        [xvalues,yvalues_No_sitrring,yvalues_stirring] = Data_Slower(k).Plot_MOT_Capture_efficiency(mf,false);
        yvalues_No_sitrring_Mat(k,:) =yvalues_No_sitrring;
        yvalues_stirring_Mat(k,:) = yvalues_stirring;
        xvalues_Mat(k,:) = xvalues;
    end
    yvalues_No_stirring = mean(yvalues_No_sitrring_Mat)
    yvalues_stirring = mean(yvalues_stirring_Mat)
    xvalues = mean(xvalues_Mat)
    
    figure
    hold on
    plot(xvalues,yvalues_No_stirring)
    plot(xvalues,yvalues_stirring)
    legend('Stirring off','Stirring on')
    title(['mf = ', num2str(mf)])
    
    folder = "Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\PhD Thesis\Figures\FermiMOT\Fig6 - Phase-SpacePlot\Data\Capture_vel\"
    save(folder +"m_f_initial = " + string(mf_vec(indx)) +  ".mat",'yvalues_No_stirring','yvalues_stirring','xvalues')
    
    y_stirring_on(indx,:)   = yvalues_stirring;
    y_stirring_off(indx,:)  =  yvalues_No_stirring;
    
end
figure
hold on
mean_stir_on = mean(y_stirring_on)
mean_stir_off = mean(y_stirring_off)

plot(-xvalues,mean(y_stirring_off))
plot(-xvalues,mean(y_stirring_on))
legend('Stirring off','Stirring on')
folder = "Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\PhD Thesis\Figures\FermiMOT\Fig6 - Phase-SpacePlot\Data\Capture_vel\"
save(folder +"m_f_initial =average " +  ".mat",'mean_stir_on','mean_stir_off','xvalues')


%%% Ploting all calculated Phase Space trayectories
close all
figure
mf_arr = -9/2:-7/2
v_vector = -linspace(1,4,20)
str_Vec = [false,true]
mf_vec = [-9/2:9/2]
for bool_str_on = str_Vec
for mf = mf_vec
    for v_initial_val = v_vector
    stirring_on = bool_str_on
    initial_mf =mf

    v = {}
    r = {}
    t = {}

    for data = Data_Slower
        paramsArray = data.paramsArray;
        paramsArray = data.remove_captured_particles();
        paramsArray = Filter_mf(data,initial_mf);

        h = paramsArray(1).UpdateStepSize;

        for indx = 1:length(paramsArray)
            ic = paramsArray(indx).initialConditions;
            stirringON = (length(paramsArray(indx).LaserArray)> 5)
            if (ic(5) ==v_initial_val && stirringON == stirring_on  )
                SimData = paramsArray(indx).SimData;

                v{end + 1 } = abs(SimData(:,2 + 3));
                r{end + 1} = SimData(:,2)*100;



                l = length(r{end});
                h = paramsArray(indx).UpdateStepSize;
                t{end + 1} =  linspace(0,l*h,l);
            end
        end
    end
    vi = [] 
    ri =  []
    tTarget = linspace(0,t{1}(end)*1.5,1000);
    for indx = 1:length(v)
        vi(:,indx) = interp1(t{indx},v{indx},tTarget);
        ri(:,indx) = interp1(t{indx},r{indx},tTarget);
    end


    figure(1)
    hold on
    meanr = mean(ri,2);
    meanv = mean(vi,2);
%     folder = "Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\PhD Thesis\Figures\FermiMOT\Fig6 - Phase-SpacePlot\Data\Phase_space_plot_Monte_carlo\"
%     save(folder +"OnlyLost mf=" + string(mf) + "v=" +string(v_initial_val)+ "Stirring_on =" + string(bool_str_on)+ ".mat",'meanr','meanv')
    plot(meanv,meanr)
    end
    end
end
figure
for mf = -9/2:9/2
    yvalues_No_sitrring_Mat = [];
    yvalues_stirring_Mat = [];
    xvalues_Mat = [];
    for k = 1:length(Data_Slower)
        [xvalues,yvalues_No_sitrring,yvalues_stirring] = Data_Slower(k).Plot_MOT_Capture_efficiency(mf,false);
        yvalues_No_sitrring_Mat(k,:) =yvalues_No_sitrring;
        yvalues_stirring_Mat(k,:) = yvalues_stirring;
        xvalues_Mat(k,:) = xvalues
    end
    yvalues_No_stirring = mean(yvalues_No_sitrring_Mat)
    yvalues_stirring = mean(yvalues_stirring_Mat)
    xvalues = mean(xvalues_Mat)
    
    hold on
    plot(xvalues,yvalues_No_stirring)
    plot(xvalues,yvalues_stirring)
    legend('Stirring off','Stirring on')
    title(['mf = ', num2str(mf)])
end

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



