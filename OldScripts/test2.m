a  = rand(1000,1000);
agpu = gpuArray(a);

tic 
a.^2;
toc

tic 
agpu.^2;
toc

% 
tic 
arrayfun(@testfn,agpu);
toc
