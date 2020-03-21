function newobj = remove(obj,I,varargin)
%REMOVE removes a point from a qtpl object
%
%   REMOVE(T,I)   remove the point of index I from qtpl object T
%
%   REMOVE(T,I,w)   specifies which at which frequencies to remove
%
%
%   Inputs:
%   T           qtpl object array
%   I           index of the point 
%   w           vector of frequencies to import (def=all)
%

p = inputParser;
addRequired(p,'obj',...
    @(x) validateattributes(x,{'qtpl'},{'nonempty'}));
addRequired(p,'I',...
    @(x) validateattributes(x,{'numeric'},{'positive','vector'}));
addOptional(p,'w',[],...
    @(x) validateattributes(x,{'double'},{'positive','vector'}))
parse(p,obj,I,varargin{:})
w = p.Results.w;

w_tpl = [obj.frequency];

if isempty(p.Results.w)
    N = length(w_tpl);
    tplidx = 1:N;
elseif all(ismember(w,w_tpl))
    N = length(p.Results.w);
    [~,tplidx] = ismember(w,w_tpl);
else
    error('at least one of specified frequencies not found in tplfile')
end

newobj = obj;

for k=1:N
    
    T = obj(tplidx(k)).template;
    par = obj(tplidx(k)).parameters;
    
    n = length(T);
    idx = 1:n;
    Tnew = T(idx~=I);
    parnew = par(:,idx~=I);
    
    newobj(tplidx(k)).template = Tnew;
    newobj(tplidx(k)).parameters = parnew;

end

