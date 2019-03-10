function T = union(T1,T2)
%UNION union of qtpl objects
% 
%   it is assumed that T1 and T2 share the same parameters.

if ~isa(T2,'qtpl') 
    error('inputs to tpl.union must be qtpl objects')
end

w1 = [T1.frequency];
w2 = [T2.frequency];
[~,idx] = ismember(w2,w1);
T = qtpl(length(idx));
for k=1:length(idx)
    w = w2(k);
    tpl = [T1(w1==w).template ; T2(k).template];
    par = [T1(w1==w).parameters T2(k).parameters];
    T(k).frequency = w;
    T(k).template = tpl;
    T(k).parameters = par;
    T(k).parNames = T1.parNames;
end

end

