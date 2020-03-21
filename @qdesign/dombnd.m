function dombnd = dombnd(obj)
%DOMBND computes the dominante HS bounds 
%   Used to compute dominante HS bounds  for qdesign object with bounds from 
%   multiple specifications. 
%
%   Called by: qdes.showbnd('dom',...)
%
%   

bnd = obj.bnd;
w = unique([bnd.w]);

dombnd = struct('name','dombnd','w',w,'c',[]);
N = 360; % number of points along contour

for kw=1:length(w)
    a = zeros(length(bnd),N);
    theta = linspace(-2*pi+0.01,-0.01,N);
    for k = 1:length(bnd)
        [~,idx ] = ismember(w(kw),bnd(k).w);
        if isempty(idx) 
            continue
        end
        b = n2c(bnd(k).c{idx})+1;   % translate bounds to complex form with 
                                    % origin at (-1)
        if all(isnan(b))
            continue
        end
        if any([min(real(b))>=0, max(real(b))<=0, min(imag(b))>=0, max(imag(b))<=0])
            warning('point -1 is not inside the tpl for w=%g, skipping.',w(kw))
            continue
        end
        
        %figure(2), clf
        %scatter(real(b),imag(b)); hold on
        
        %inan = [find(isnan(b)) length(b)+1]; % nans indicate holes in bounds
        %angB = angle(b(1:inan(1)-1)); 
        if any(isnan(b))
            b = b(~isnan(b)); % b without nans
            %K = convhull(real(b),imag(b));  % nans indicate holes in bounds
            %angB = angle(b(K));
            %b = b(K);
            angB = angle(b);
        else
            angB = angle(b); 
        end
        angB(angB>0)=angB(angB>0)-2*pi; % force angles to be in range [-2*pi 0]
        [~,I] = unique(angB);
        a(k,:) = interp1(angB(I),abs(b(I)),theta,'linear','extrap');
        
        %b_inter = a(k,:).*exp(1i*theta);  
        %plot(real(b_inter),imag(b_inter))
    end
    amax=max(a);
    C = amax.*exp(1i*theta); % complex form w.r.t to (-1)
    BND = c2n(C-1);          % shift to (-1) an trasfrom into Nichols form
    figure(1), plot(real(BND),imag(BND),'color',obj.col(kw,:)); hold on    
    dombnd.c{kw} = BND;
    dombnd.c{kw} = closecontour(BND);
end
%xlim([-360 0])
%ngrid

end


function CC = closecontour(CO)
% make sure that bounds that are supposed to be closed contours are indeed 
% closed 

if abs(CO(1)-CO(end)) < 20
    CC = [CO CO(1)];
else
    CC = CO;
end

end
