function InstallGurobi(gurobiPath)
if nargin < 1
gurobiPath = '/opt/gurobi751/linux64/'; % Please specify another path if you would like to use a different Gurobi installation
end
disp([ 'Gurobi path is set to ' gurobiPath]);
addpath([gurobiPath 'matlab/']);
gurobi_setup;
end