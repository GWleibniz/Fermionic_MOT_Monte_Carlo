% % % % % Data = DataClass('Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\Data_1.mat')
Data_MOT = DataClass("Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\MOT\raw\run_data_1\","Data_0.mat")
% % % % close all
% % %  axis_arr = ['x','y','z','x','y','z']
% % % %% Get spatial width
% % % for axis = 1:3
% % %     Data_MOT.Histogram_x_axis =linspace(-20,20,100) ;
% % %     [x,y] = Data_MOT.PlotHistogram(axis);
% % %     [ymax,loc] = max(y);
% % %     y_norm = y/ymax;
% % %     fun = @(p,x) Data_MOT.gaussian1Dfn(p,x);
% % %     x0 = [x(loc),2];
% % %     x_values = x(1:end-1)+ diff(x(1:2))/2;
% % %     fit_res = lsqcurvefit(fun,x0,x_values,y_norm);
% % %     yfit = fun(fit_res,x);
% % %     hold on 
% % %     plot(x,yfit*max(y))
% % %     legend(['\sigma', axis_arr(axis),' = ',num2str(fit_res(2)), 'mm' ])
% % % 
% % % end
% % % % % % % % %% Get temperature width
% % % for axis = 4:6
% % %     cte = ParamClass;
% % %     Data_MOT.Histogram_x_axis =linspace(-0.4,0.4,100)*1000 ;
% % %     [x,y] = Data_MOT.PlotHistogram(axis);
% % %     [ymax,loc] = max(y);
% % %     y_norm = y/ymax;
% % %     fun = @(p,x) Data_MOT.gaussian1Dfn(p,x);
% % %     x0 = [x(loc),100];
% % %     x_values = x(1:end-1)+ diff(x(1:2))/2;
% % %     fit_res = lsqcurvefit(fun,x0,x_values,y_norm);
% % %     yfit = fun(fit_res,x);
% % %     hold on 
% % %     plot(x,yfit*max(y))
% % %     T = (fit_res(2)/1000)^2*cte.m/cte.kB*1e6;
% % % 
% % %     legend(['\sigma', axis_arr(axis),' = ',num2str(fit_res(2)),' mm/s temperature is',num2str(T),'\mu k' ])
% % % 
% % % %     display(['Mot width on slected axis = ',axis_arr(axis) ,' ' ,num2str(fit_res(2)), 'mm' ])
% % % 
% % % end
% % % 
% % % 

%% Get temperature
axis_arr = 'xyzxyz'
for axis = 1:3
    [T,Delta_t] = Data_MOT.get_avg_temp(axis) ;
    
    display(['Mot Temperature on axis = ',axis_arr(axis) ,' ' ,num2str(T*1e6),'(',num2str(Delta_t*1e6),')', 'uk' ])
end

%% Get size
for axis = 1:3
    [T,Delta_t] = Data_MOT.get_avg_size(axis) ;
    
    display(['Mot size on axis = ',axis_arr(axis) ,' ' ,num2str(T*1e3),'(',num2str(Delta_t*1e3),')', 'mm' ])
end

% % % % 
% % % % %%ZMOT offsett
% % % % Data_MOT_offset = DataClass("Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\MOT\raw\Old_runs\r3\","Data_0.mat")
% % % % 
% % % % %% Get temperature
% % % % for axis = 1:3
% % % %     [T,Delta_t] = Data_MOT_offset.get_avg_temp(axis) ;
% % % %     
% % % %     display(['Mot Temperature on axis = ',axis_arr(axis) ,' ' ,num2str(T*1e6),'(',num2str(Delta_t*1e6),')', 'uk' ])
% % % % end
% % % % 
% % % % %% Get size
% % % % for axis = 1:3
% % % %     [T,Delta_t] = Data_MOT_offset.get_avg_size(axis) ;
% % % %     
% % % %     display(['Mot size on axis = ',axis_arr(axis) ,' ' ,num2str(T*1e3),'(',num2str(Delta_t*1e3),')', 'mm' ])
% % % % end


%%xMOT power
Data_MOT_x_power = DataClass("Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\MOT\raw\Old_runs\r7\","Data_0.mat")

%% Get temperature
for axis = 1:3
    [T,Delta_t] = Data_MOT_x_power.get_avg_temp(axis) ;
    
    display(['Mot Temperature on axis = ',axis_arr(axis) ,' ' ,num2str(T*1e6),'(',num2str(Delta_t*1e6),')', 'uk' ])
end

%% Get size
for axis = 1:3
    [T,Delta_t] = Data_MOT_x_power.get_avg_size(axis) ;
    
    display(['Mot size on axis = ',axis_arr(axis) ,' ' ,num2str(T*1e3),'(',num2str(Delta_t*1e3),')', 'mm' ])
end
% % % 
% % % %% Plot average MF ocupation
% Data_MOT.Plot_avg_mf_ocupation()
[mean_mf_occupation,std_mf_occupation,mf_val] = Data_MOT.Get_avg_mf_ocupation()
folder = 'Y:\Groups\strontium\People\Rodrigo Gonzalez Escudero\PhD Thesis\Figures\FermiMOT\Figures Python scripts\Data\Mean_mf_occ\'
filename = 'data.mat'
save([folder,filename],'mean_mf_occupation','std_mf_occupation','mf_val')