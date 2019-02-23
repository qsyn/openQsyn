function dombnd = dombnd(obj)
%DOMBND computes the dominante HS bounds 
%   Used to compute dominante HS bounds  for qdes object with bounds from 
%   multiple specifications. 
%
%   Called by: qdes.sheobnd('dom',...)
%
%   

bnd = obj.bnd;
w = unique([bnd.w]);

dombnd = struct('name','dombnd','w',w,'c',[]);

for kw=1:length(w)
    a = zeros(length(bnd),360);
    theta = linspace(-2*pi+0.01,-0.01,360);
    for k = 1:length(bnd)
        [~,idx ] = ismember(w(kw),bnd(k).w);
        if isempty(idx) 
            continue
        end
        b = n2c(bnd(k).c{idx})+1; % move origing to (-1)
        if all(isnan(b))
            continue
        end
        inan = [find(isnan(b)) length(b)+1]; % nans indicate holes in bounds
        angB = angle(b(1:inan(1)-1)); 
        angB(angB>0)=angB(angB>0)-2*pi;
        [~,I] = unique(angB);
        a(k,:) = interp1(angB(I),abs(b(I)),theta,'linear','extrap');
    end
    amax=max(a);
    C = amax.*exp(1i*theta); % complex form w.r.t to (-1)
    BND = c2n(C-1);          % shift to (-1) an trasfrom into Nichols form
    %plot(real(BND),imag(BND),'color',obj.col(kw,:)); hold on    
    %dombnd.c{kw} = BND;
    dombnd.c{kw} = closecontour(BND);
end
%xlim([-360 0])
%ngrid

end


function CC = closecontour(CO)
% make sure that bounds that are supposed to be closed contours are indeed 
% closed 
if (max(real(CO))>-0.5) || (min(real(CO))<-355.5)
    CC=CO;
else
    CC = [CO CO(1)];
end

end
