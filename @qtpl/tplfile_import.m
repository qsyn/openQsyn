function [ T ] = tplfile_import(filename)
%TPLFILE_IMPORT imports templates from existing Qsyn *.tpl file
%
%   [ T ] = tplfile_import(FILENAME)    import tpl from FILENAME.tpl into a 
%   qtpl object T

S=load([filename,'.tpl'],'-mat');

N = size(S.w_tpl,1);
T = qtpl(N);

inom = ismember(S.w_nom,S.w_tpl(:,1)); % get nominal point at each frequency
nom = S.nom(inom);
par_nom = S.par_nom.';

for n=1:N
    
    k = S.w_tpl(n,2);
    T(n).frequency = S.w_tpl(k,1);
    T(n).parameters = [par_nom S.(sprintf('par_%i',k))];
    T(n).template = [nom(n) ; S.(sprintf('t_w%i',k))];
    
end


end

