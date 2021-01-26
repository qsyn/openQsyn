function newObj = reduce( obj,dist,w )
%TPLREDUCE reduce template size by removing points
% 
%Usage: 
%
%newObj = TPLREDUCE(obj,dist)   return a reduced qtpl array by
%removing points with distance less than 'dist' in each object (freq)
% 
%newObj = TPLREDUCE(obj,dist,w)   specifies frequecnies to reduce, other
%qtpl objects in array are returned as they were
%
%
%Inputs:  
%
%	dist         [phase mag] demanded for minimal distance 
%                between adjacent points normelized by the templkate width 
%                and height
%
%	w            vector of frequencies [rad/s], default = all frequencies
%
%Example:
%   
%   T.reduce([0.1 0.05])   reduces points by 10% along phase and by 5%
%   along magnitude

plot_on=1;

if nargin==0, disp('newObj = reduce(obj,dist)');return;end
if ~exist('w','var'), w=[]; end

if dist(1)>=1 || dist(2)>=1 || dist(1)<=0 || dist(2)<=0
    error('2nd argument must be a vector of length 2 with values in range (0,1)');
end

w_tpl = [obj.frequency];
if isempty(w)
    w = w_tpl; 
end
[~,w_idx] = ismember(w,w_tpl);

newObj = obj;

for k = w_idx
    
    T = obj(k).template;
    par = obj(k).parameters;
    phi=real(T); mag=imag(T);
    dx=abs(max(phi)-min(phi))*dist(1);
    dy=abs(max(mag)-min(mag))*dist(2);
    X=min(phi):dx:max(phi);
    Y=min(mag):dy:max(mag);
    inx=zeros((length(X)-1)*(length(Y)-1),1);
    kk=1;
    for kx=1:length(X)-1
        for ky=1:length(Y)-1
            in=inpolygon(phi,mag,X([kx,kx,kx+1,kx+1]),Y([ky,ky+1,ky+1,ky]));
            if sum(in)>0
                inx(kk)=find(in,1,'first');               
            else
                inx(kk)=nan;
            end
            % plotting
            if plot_on && sum(in)>1
                fill(X([kx,kx,kx+1,kx+1]),Y([ky,ky+1,ky+1,ky]),[1 0 0]); hold on
                scatter(phi(in),mag(in),'b*'); 
                scatter(phi(inx(kk)),mag(inx(kk)),'k*')
            end
            kk=kk+1;
        end
    end
    inx(isnan(inx))=[]; % remove nans     
    
    newObj(k).parameters = par(:,inx);
     newObj(k).template = T(inx);

end


end