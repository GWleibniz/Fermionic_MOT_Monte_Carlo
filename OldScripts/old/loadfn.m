params = ParamClass;
files = dir('Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\*.mat');
for k = 1:length(files)

    path(k) = "Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\" + files(k).name;

    

end
a = params.LoadParams(path(1));
paramsArray = repmat(a,1,length(files));
for k = 1:length(files)
    k
    paramsArray (k) = params.LoadParams(path(k));
end
% paramsArray = arrayfun(@(k)params.LoadParams(path(k)),1:length(path));
path = "Z:\Groups\strontium\People\Rodrigo Gonzalez Escudero\MATLAB\MOT_MOntecarlo\save\Slower\raw\" + "Data1_2v2.mat";
loadedData = paramsArray;
save(path,'loadedData')
