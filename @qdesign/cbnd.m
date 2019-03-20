function obj = cbnd(obj,spcname,w,spcfunc)
%CBND Summary of this function goes here
%   Detailed explanation goes here

if nargin<3, w=[]; end
if nargin<4, spcfunc=[]; end

switch spcname
    case 'odsrs', spcfunc = @fodsrs; 
    case 'rsrs',  spcfunc = @frsrs; 
    %otherwise, error('unrecognized specification %s',spcname);
end
     
if isempty(spcfunc)
    error('for non standard spc names an spc. function must be specified')
end

if isempty(w), w = [obj.tpl.frequency]; end

is = find(strcmp({obj.spc.name},spcname));
spec = obj.spc(is);

it = ismember([obj.tpl.frequency],w);
tpl = obj.tpl(it);

bnd.name = spcname;
bnd.w = w;  

fprintf('Calculating bounds for %s \n',spcname);
bnd.c = cbnd1(tpl,spcfunc,spec);

if isempty(obj.bnd)
    obj.bnd = bnd; 
elseif any(strcmp({obj.bnd.name},spcname))
    obj.bnd(strcmp({obj.bnd.name},spcname)) = bnd;
else
    obj.bnd(end+1) = bnd;
end


end


function c = cbnd1(tpls,spcfunc,spec)
%CBND1 computes a single bound at all frequencies

%default initial grid
gphase0=-360:10:0;
gmag0=-50:5:50;

w = [tpls.frequency]';
upper = interp1(spec.frequency,spec.upper,w);
if isempty(spec.lower)
    lower = [];
else
    lower = interp1(spec.frequency,spec.lower,w);
end
spcval = [upper lower];

c = cell(length(tpls),1);
for k=1:length(tpls)
    fprintf('--> w(%i) = %g [rad/s]\n',k,w(k));
    tpl = tpls(k).template;   
    c0 = qdesign.makebnd(tpl,spcfunc,spcval(k,:),gphase0,gmag0);   
    
    % refined grid
    gphase = floor(( max(min(real(c0))-2,-360) ) : 2 : ( min(max(real(c0))+2,0) ));
    gmag = floor(( min(imag(c0))-2 ) : 2 : ( max(imag(c0))+2 ));
    %bphase  = floor([max(min(real(c0))-2,-360) min(max(real(c0))+2,0)]);
    %bmag = floor([min(imag(c0))-2 max(imag(c0))+2]);
    
    % refined grid around the initial points
    %gphase = repmat(real(c0),1,20) + repmat(repmat([-2 -1  0  1  2],1,4),1,length(c0));
    %gmag = repmat(imag(c0),1,20) + repmat([-2*ones(1,4) -1*ones(1,4) 0*ones(1,4) 1*ones(1,4) 2*ones(1,4)],1,length(c0));
    
    c{k} = qdesign.makebnd(tpl,spcfunc,spcval(k,:),gphase,gmag);   
end

end

