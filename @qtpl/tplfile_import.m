function [ T ] = tplfile_import(filename,varargin)
%TPLFILE_IMPORT imports templates from existing Qsyn *.tpl file
%
%   [ T ] = tplfile_import(filename,w)   imports tpl from FILENAME.tpl into 
%   qtpl object T
%
%   Inputs:
%   filename    char array of the tpl file name w/o the *.tpl extension
%   w           vector of frequencies to import (def=all)
%

p = inputParser;
addRequired(p,'filename',...
    @(x) validateattributes(x,{'char'},{'nonempty'}));
addOptional(p,'w',[],...
    @(x) validateattributes(x,{'double'},{'positive','vector'}))
parse(p,filename,varargin{:})
w = p.Results.w;

S=load([p.Results.filename,'.tpl'],'-mat');
w_tpl = sortrows(S.w_tpl,2); 

if isempty(p.Results.w)
    N = size(w_tpl,1);
    idx = 1:N;
elseif all(ismember(w,w_tpl(:,1)))
    N = length(p.Results.w);
    [~,idx] = ismember(w,w_tpl(:,1));
else
    error('at least one of specified frequencies not found in tplfile')
end

T = qtpl(N); % pre-allocation

nom = [];
if isfield(S, 'w_nom')
    %inom = ismember(S.w_nom,S.w_tpl(:,1)); % get nominal point at each frequency
    %nom = S.nom(inom);
    [isnom,inom]=ismember(w_tpl(idx,1),S.w_nom);
else
    isnom = zeros(N,1);
end

par_nom = [];
if isfield(S, 'par_nom')
    par_nom = S.par_nom.';
end

for n=1:N    
    k = w_tpl(idx(n),2);
    T(n).frequency = w_tpl(k,1);
    if isnom(n)
        nom = S.nom(inom(n));
    else
        nom=[];
    end
    %fprintf('t_w%i \n',k);
    tpl = [nom ; S.(sprintf('t_w%i',k))];
    T(n).template = unwrap(real(tpl)*pi/180)*180/pi + 1i*imag(tpl); %unwrap
    if isfield(S,sprintf('par_%i',k))
        T(n).parameters = [par_nom S.(sprintf('par_%i',k))];
    else
        T(n).parameters = [par_nom nan(size(S.(sprintf('t_w%i',k))))];
    end
end


end

