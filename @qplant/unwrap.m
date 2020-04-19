function obj = unwrap(obj)
%UNWRAP unwarps a qplant nominal response and templates

if isempty(obj.nominal)
    error('nominal plant must be computed before unwrapping')
end
if isempty(obj.templates)
    error('no template sto unwrap')
end

nom = obj.nominal.unwrap;







end

