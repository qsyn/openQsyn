function oqsyn_ini()
% OQSYN initiates OpenQsyn  
%
% Open Qsyn is a modern open source toolbox for QFT control synthesis
% The project is hosted on GitHub: https://github.com/rubindan/openQsyn

% Author: Daniel R. 31-Jan-2019


fprintf('Open Qsyn toolbox is now initiated \n\n')

% all it does is to add utilities and doc to the path
fullpath=mfilename('fullpath');
homepath=fullpath(1:end-9);
addpath([homepath,'utilities'])
addpath(genpath([homepath,'doc']))

builddocsearchdb('C:\Users\RubinDan\Documents\MATLAB\Qsyn\oo\doc\html')

end
