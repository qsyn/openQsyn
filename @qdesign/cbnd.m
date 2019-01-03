function obj = cbnd(obj,spcname,w)
%CBND Summary of this function goes here
%   Detailed explanation goes here

if nargin<3, w=[]; end

switch spcname
    case 'odsrs', spcfunc = @fodsrs; 
    case 'rsrs',  spcfunc = @frsrs; 
    otherwise, error('unrecognized specification %s',spcname);
end
     
if isempty(w), w = [obj.tpl.frequency]; end

is = find(strcmp({obj.spc.name},spcname));
spec = obj.spc(is);

it = ismember([obj.tpl.frequency],w);
tpl = obj.tpl(it);

bnd.name = spcname;

fprintf('Calculating bounds for %s \n',spcname);
bnd.c = cbnd1(tpl,spcfunc,spec);

bnd.w = w;  
obj.bnd = bnd;    % single bound! in future must be changed to bounds from multiple specs.


end


function c = cbnd1(tpls,spcfunc,spec)
%CBND1 computes a single bound at all frequencies

%default initial grid
gphase0=-360:10:0;
gmag0=-50:5:50;

spcval = [spec.upper spec.lower];

c = cell(length(tpls),1);
for k=1:length(tpls)
    fprintf('--> w(%i) \n',k);
    tpl = tpls(k).template;   
    c0 = qdesign.makebnd(tpl,spcfunc,spcval(k),gphase0,gmag0);   
    
    % refined grid
    %gphase = floor(( max(min(real(c0))-2,-360) ) : 1 : ( min(max(real(c0))+2,0) ));
    %gmag = floor(( min(imag(c0))-2 ) : 1 : ( max(imag(c0))+2 ));
    bphase  = floor([max(min(real(c0))-2,-360) min(max(real(c0))+2,0)]);
    bmag = floor([min(imag(c0))-2 max(imag(c0))+2]);
    for g = length(c0)+2
        
        gphase = 
        
    end
    
    c{k} = qdesign.makebnd(tpl,spcfunc,spcval(k),gphase,gmag);   
end

end

