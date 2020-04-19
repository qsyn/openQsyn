function obj = unwrap(obj)
%UNWRAP unwarps a qplant nominal response and templates

if isempty(obj.nominal)
    error('nominal plant must be computed before unwrapping')
end
if isempty(obj.templates)
    error('no template sto unwrap')
end

nom = obj.nominal.unwrap;
tpl = obj.templates;

for k = 1:length(tpl)
    
    [~,I] = min(abs(nom.frequency - tpl(k).frequency));
    n_r=ceil((real(tpl(k).template(1))-real(nom.response(I)))/360);
    tpl(k).template = tpl(k).template-(n_r)*360;

end

obj.nominal = nom;
obj.templates = tpl;

end

