function [w]   = beamwaist(w0,z)%Initialize variables
params = ParamClass;
zr = params.kred*w0.^2/2;
w = w0.*sqrt(1+(z./zr).^2);
end