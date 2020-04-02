function tplfile_export(obj,filename)
%TPLFILE_EXPORT exports templates to an old qsyn *.tpl file

N = length(obj.templates);

w_tpl = [obj.templates.frequency ; 1:N].';
w_nom = obj.nominal.w_nom;
nom = obj.nominal.nic;
par_nom = [obj.pars.nominal].';
Par_name = sprintf('%s ',obj.templates(1).parNames{:});
info = 'File exported from open Qsyn';

for k=1:N

    par_k = obj.templates(k).parameters;
    tpl_k = obj.templates(k).template;
    eval(sprintf('par_%i = %g',k,par_k));
    eval(sprintf('t_w%i = %g',k,tpl_k));
     
end

clear nargin nargout 
save('filename')

end

