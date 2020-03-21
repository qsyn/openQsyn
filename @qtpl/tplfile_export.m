function tplfile_export(obj,filename,varargin)
%TPLFILE_EXPORT exports templates into a Qsyn *.tpl file
%
%   TPLFILE_EXPORT(T,filename)   exports tpl T into <filename.tpl>
%
%   TPLFILE_EXPORT(T,filename,w)   specifies which frequencies to export
%
%
%   Inputs:
%   T           qtpl object array
%   filename    char array of the tpl file name w/o the *.tpl extension
%   w           vector of frequencies to import (def=all)
%

p = inputParser;
addRequired(p,'obj',...
    @(x) validateattributes(x,{'qtpl'},{'nonempty'}));
addRequired(p,'filename',...
    @(x) validateattributes(x,{'char'},{'nonempty'}));
addOptional(p,'w',[],...
    @(x) validateattributes(x,{'double'},{'positive','vector'}))
parse(p,obj,filename,varargin{:})
w = p.Results.w;

w_tpl = [obj.frequency];

if isempty(p.Results.w)
    N = length(w_tpl);
    idx = 1:N;
elseif all(ismember(w,w_tpl(:,1)))
    N = length(p.Results.w);
    [~,idx] = ismember(w,w_tpl(:,1));
else
    error('at least one of specified frequencies not found in tplfile')
end

w_tpl = [w_tpl' (1:N)'];
parnames = obj(1).parNames;
info = 'exported from OpenQsyn';
data = date;

save(filename,'w_tpl','parnames','info','data');


for k=1:N
    
    I = idx(k);
    eval(sprintf('par_%i = obj(I).parameters;',k));
    eval(sprintf('t_w%i = obj(I).template;',k));
    save(filename,sprintf('par_%i',k),sprintf('t_w%i',k),'-append');
    
end

movefile([filename,'.mat'],[filename,'.tpl']);





end

