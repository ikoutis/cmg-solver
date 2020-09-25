
addpath(pwd); 
cd mex;
makeSolverMex;
eval('mex mx_splitforest_.c');
addpath(pwd);
cd ..;
savepath;

