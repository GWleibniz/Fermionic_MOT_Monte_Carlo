classdef DataClass
   properties
    path;
    paramsArray;
    Histogram_x_axis = linspace(-10,10,50);
    batchSize = 1
   end
   methods
       
    function obj = DataClass(path,filename,batch_size)
        
      if nargin >2
          obj.batchSize = batch_size
      end
          
      obj.path = path;
      if(isfile(obj.path + filename))
          obj.paramsArray = obj.loadData(filename);
      else
          obj.paramsArray = obj.Construct_params_array(filename);
      end
      obj.paramsArray =  filter_bad_params(obj)
    end
    function newParamsArray =  filter_bad_params(obj)
        paramsArray = obj.paramsArray;
        newParamsArray=ParamClass.empty
        for indx = 1:length(paramsArray)
            if ~ isempty(paramsArray(indx).SimData)
                newParamsArray(end+1) = paramsArray(indx);
            end
        end
    end
    
    function [xvalues,yvalues_No_sitrring,yvalues_stirring] =  Plot_MOT_Capture_efficiency(obj,mf,do_plot)

    Initial_mf = zeros (size(obj.paramsArray));
    Initial_v = zeros (size(obj.paramsArray));
    TrappedBool = zeros (size(obj.paramsArray));
    stirringOnBool = zeros (size(obj.paramsArray));

    for indx = 1:length(obj.paramsArray)
        param = obj.paramsArray(indx);
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

        yvalues_stirring(indx) = sum(TrappedBool(valueBool_arr & stirringOnBool ))/  number_stirring_on;
        yvalues_No_sitrring(indx) = sum(TrappedBool(valueBool_arr & (~stirringOnBool )))/  number_stirring_off;
        
    end
    
    if do_plot
        figure
        hold on
        plot(xvalues,yvalues_No_sitrring)
        plot(xvalues,yvalues_stirring)
        legend('Stirring off','Stirring on')
        title(['mf = ', num2str(mf)])
    end
    
    end
 function [mean_mf_occupation,std_mf_occupation,mf_val] = Get_avg_mf_ocupation(obj)
    paramsArray = obj.paramsArray;
    mf_arr = [];
    excited_arr = [];
        counter = 1
        for param = paramsArray
        Atomic_state = param.SimData(1000:end,7);
        mf_arr = param.SimData(1:end,8);
        mf_arr = mf_arr(Atomic_state==0);
        counter =counter + 1;       
        number_of_values = prod(size(mf_arr));
        mf_val = -9/2:1:9/2;
            for i = 1:length(mf_val)
                n_mf(i) = length(mf_arr(mf_arr ==mf_val(i)))/number_of_values;
                
            end
            mat_mf (counter,:) = n_mf;
        end
        
        mean_mf_occupation = mean(mat_mf)
        std_mf_occupation =  std(mat_mf)
        for i = 1:length(mf_val)
            display("Percentage of atoms that occupy the mf = " + num2str(mf_val(i)) + " is" + num2str(mean_mf_occupation(i)) + "(" + num2str(std_mf_occupation(i)) + ")" )
        end

        

    end   
    function Plot_avg_mf_ocupation(obj)
    paramsArray = obj.paramsArray;
    mf_arr = []
    excited_arr = []

        for param = paramsArray
            
            mf = param.SimData(1000:end,8);
            Atomic_state = param.SimData(1000:end,7);
            mf = mf(Atomic_state == 0);
            mf_arr = [mf_arr;mf];

        end
        mf_edges = (-9/2:1:9/2) - 0.5
        [Bin_counts,Bin_edges] = histcounts(mf_arr,mf_edges, 'Normalization', 'probability')
        figure
         histogram('BinCounts', Bin_counts, 'BinEdges', Bin_edges);

    end
    function [Bin_edges,Bin_counts_total]  = PlotHistogram(obj,axis)
        paramsArray = obj.paramsArray;
        Bin_edges = obj.Histogram_x_axis;
        Bin_counts_total = zeros(1,length(Bin_edges)-1);
        for param = paramsArray
            
            x=  param.SimData(10000:end,axis)*1000; 
            [Bin_counts,Bin_edges,bin] = histcounts(x,Bin_edges,'Normalization', 'probability');
            Bin_counts_total  =Bin_counts_total + Bin_counts;

        end
        figure
        histogram('BinCounts', Bin_counts_total, 'BinEdges', Bin_edges);
        
    end
    function [dev,Delta_sigma] = get_avg_size(obj,axis)
        paramsArray = obj.paramsArray;
        cte = ParamClass;
        counter = 0;
        for param = paramsArray
            counter = counter + 1;
            r=  param.SimData(1:end,axis );
            vy = param.SimData(1:end,5 );
            r = r( abs(vy)<0.15 ) ;

            r_std(counter) = std(r);
            
            
        end
            dev = mean(r_std);
            Delta_sigma = std(r_std);
    end
    function [T,Delta_T]  = get_avg_temp(obj,axis)
        paramsArray = obj.paramsArray;
        cte = ParamClass;

        v_total = [];
        counter = 0;
        for param = paramsArray
            counter = counter + 1;
            v=  param.SimData(1:end,axis +3); 
            vy = param.SimData(1:end,5 );

            v = v( abs(vy)<0.15 ) ;

            
            v_std(counter) = std(v);
            
        end
            
            T = mean(v_std.^2)*cte.m/cte.kB;
            Delta_T = std(v_std.^2)*cte.m/cte.kB;
        
    end

    function [z] = gaussian1Dfn(obj,p,x);


        z = exp(-0.5*(x-p(1)).^2./(p(2)^2));
        z = double(z);
    end
    function paramsArray = Construct_params_array(obj,filename)
        
        params = ParamClass
        files = dir(obj.path+"*.mat");
        for k = 1:length(files)
            path(k) =  obj.path+ files(k).name;
        end
        if length(files)==0
            error('No files in target directory')
        end

        a = params.LoadParams(path(1));
        paramsArray = repmat(a,1,length(files));
        steps = length(files);
        h = waitbar(0,'Please wait...');
        for k = 1:length(files)
            waitbar(k/steps)
            paramsArray (k) = params.LoadParams(path(k));
        end
        
        close(h)
        % paramsArray = arrayfun(@(k)params.LoadParams(path(k)),1:length(path));
        path_save = obj.path + filename;
        loadedData = paramsArray;
        if obj.batchSize ==1
         save(path_save,'loadedData')
        else
         obj.SaveDataInBatches(paramsArray)  
        end

    end
    function SaveDataInBatches(obj,paramsArray)  
        a = length(paramsArray)
        batchSize = obj.batchSize;
        batch =floor( linspace(1,a,batchSize+1))
        for i = 1:(batchSize)
            loadedData = paramsArray(batch(i):batch(i+1)-1);
            filename2Save = ['BatchN',num2str(i),'.mat']
            save([obj.path,filename2Save],'loadedData')
        end
    end
    function SaveData(obj,path)
        data = obj;
        save(path,data)
    end
    function paramsArray= loadData(obj,filename) 
        path = obj.path + filename ;
        DataMat = load(path);
        try
        paramsArray =  DataMat.loadedData;
        catch
        paramsArray =  DataMat.loadedData2Save;
        end

    end

    function GetNumberOfTrappedTrayectories(obj,mf)
        paramsArray = obj.paramsArray;
        stirringOfflost = 0;
        stirringOnLost = 0;
        totlaStirringoff = 0;
        totalstirring = 0;
        paramsArray = Filter_mf(obj,mf);

        for k = 1:length(paramsArray)
             mfinitial = paramsArray(k).initialConditions(end);
                if length(paramsArray(k).LaserArray) >6 
                    stirringOnLost = stirringOnLost + (paramsArray(k).atomlostBolean  == "Atom Lost");
                    totalstirring =  totalstirring+1;

                else
                     stirringOfflost = stirringOfflost + (paramsArray(k).atomlostBolean == "Atom Lost");
                     totlaStirringoff = totlaStirringoff + 1;
                end  
        end
        

        disp('Sittiring on % of lost atoms')
        display(num2str(100*stirringOnLost/totalstirring))



        disp('Sittiring off % of lost atoms')
        display(num2str(100*stirringOfflost/totlaStirringoff))

    end
    function do_PS_plot(obj,i,axis)
        close all
        SimData = obj.paramsArray(i).SimData;
        
        v = abs(SimData(:,axis + 3));
        r = SimData(:,axis)*100;
        figure
        plot(v,r)
        xlim([0 6])
        ylim([-2.5 10])
    end
    function newParamsArray = remove_captured_particles(obj)
        paramsArray = obj.paramsArray;
        newParamsArray=ParamClass.empty
        for indx = 1:length(paramsArray)
            if paramsArray(indx).atomlostBolean == "Atom Lost"
                newParamsArray(end+1) = paramsArray(indx);
            end
        end
    end
    function newParamsArray = remove_lost_particles(obj)
        paramsArray = obj.paramsArray;
        newParamsArray=ParamClass.empty
        for indx = 1:length(paramsArray)
            if paramsArray(indx).atomlostBolean == "Atom Slowed"
                newParamsArray(end+1) = paramsArray(indx);
            end
        end
    end
    function newParamsArray = Filter_mf(obj,mf)
        paramsArray = obj.paramsArray;
        newParamsArray=ParamClass.empty;
        for indx = 1:length(paramsArray)
            
            if paramsArray(indx).initialConditions(end) ==  mf
                newParamsArray(end+1) = paramsArray(indx);
            end

        end
    end
    function do_PSD_plot_Average_all(obj,axis,mf)
        paramsArray = obj.paramsArray;
        paramsArray = remove_lost_particles(obj);
        paramsArray = Filter_mf(obj,mf);

        h = paramsArray(1).UpdateStepSize;
        tTarget = linspace(0,900*h,1000);
        
        for indx = 1:length(paramsArray)
            
            SimData = paramsArray(indx).SimData;
            v{indx} = abs(SimData(:,axis + 3));
            r{indx} = SimData(:,axis)*100;
            l = length(r{indx});
            h = paramsArray(indx).UpdateStepSize;
            t{indx} = linspace(0,l*h,l);
            vi(:,indx) = interp1(t{indx},v{indx},tTarget);
            ri(:,indx) = interp1(t{indx},r{indx},tTarget);
        end
        
% %         for indx = 1:length(paramsArray)
% %             hold on
% %             figure(1)
% %             plot(ri(:,indx))
% %             plot(vi(:,indx))
% %             figure(2)
% %             plot(v{indx},r{indx})
% % 
% %         end
% %         
        figure
        meanr = mean(ri,2);
        meanv = mean(vi,2);
        stdv = abs(std(vi,0,2));
        stdr = std(ri,0,2);
        hold on
        plot(meanv+stdv/2,meanr-stdr/2);
        plot(meanv-stdv/2,meanr+stdr/2);
        plot(meanv-stdv/2,meanr-stdr/2);
        plot(meanv+stdv/2,meanr+stdr/2);
        
        plot(meanv,meanr)

        

%         for indx = 1:length(paramsArray)
%             figure
%             hold on
%             plot(v{indx},r{indx})
%             plot(vi{indx},ri{indx})
%             xlim([0 6])
%             ylim([-2.5 10])
%         end
    end
   end
end