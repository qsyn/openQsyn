function tpl = crff(obj,w_tpl,options)
%CRFF computes rff template
%   Detailed explanation goes here


disp(['Calculating templates using the Real Factored Form method']);

ndist = options.Tacc; 
angDist = ndist(1);

if (abs(360-fix(360/angDist)*angDist)>eps)
    error('360 must be an integer multiple of the angular resolution, dist(1) [deg]');
end

tpl = qtpl(length(w_tpl));
for ii=1:length(w_tpl)
    w = w_tpl(ii);
    disp(['--> for w=',num2str(w),' [rad/sec]']);
    
    
    % compute numerator
    t_=[];
    for k = 1:length(obj.num)
        rffElement = obj.num(k);
        if any(strcmp(rffElement.type,{'gain','delay','uns','int'}))
            t = rffElement.rffel(w,angDist); 
        elseif strcmp(rffElement.type,'poly')
            f = rffElement.qrff2func;
            c = f(w*1j);
            t = c2n(c,-180);
        else
            if isempty(rffElement.par2)
                t = rffElement.rffpz(w,'z',angDist);
            else
                t = rffElement.rffcpz(w,'z',angDist);
            end
        end
        t_ = qrff.rffmul(t,t_,angDist);
    end
    
    if ~isempty(t_)
        num_r=t_;
    else
        num_r=0;
    end
    
    % compute denumerator
    t_ = [];
    for k = 1:length(obj.den)
        rffElement = obj.den(k);
        if any(strcmp(rffElement.type,{'gain','delay','uns','int'}))
            t = rffElement.rffel(w,angDist); 
        elseif strcmp(rffElement.type,'poly')
            f = rffElement.qrff2func;
            c = f(w*1j);
            t =  c2n(c,-180);
        else
            if isempty(rffElement.par2)
                t = rffElement.rffpz(w,'z',angDist); % note that z is used anyways! 
            else
                t = rffElement.rffcpz(w,'z',angDist); % note that z is used anyways!
            end
        end
        t_ = qrff.rffmul(t,t_,angDist);
    end
    
    [n,~]=size(t_);
    if n == 0
        den_r = 0;        
    elseif n>2
        den_r = [t_((n/2+1:n),:) ; t_((1:n/2),:)];         
    else
        den_r = t_;
    end
    
    t_ = qrff.rffmul(num_r,-den_r,angDist);
       
    [n,nn]=size(t_);
    if ((n*nn) > 1)
        min_t=t_(1:n/2);
        max_t=t_(n:-1:(n/2+1));
        [bb_,index_min]=sort(real(min_t));
        [bb_,index_max]=sort(real(max_t));
        min_t = min_t(index_min);
        max_t = max_t(index_max);
        max_t = max_t(length(max_t):-1:1);
        t_ = [min_t ; max_t];
        t_ = qrff.cltmp(t_,ndist); 
    end
    %if length(n_dif)>1  % add uncertain differentiators
    %    t_=tplop('A+B',t_,c2n((1j*w_tpl(ii)).^n_dif(1:(length(n_dif)-1)),'unwrap'));
    %end
    %tnom=real(tab_le1([w_nom(:),nom(:)],w_tpl(ii))); %Mattias 960903
    %t_=qunwrap(t_);%mattias 961122
    %t_=t_+round((-max(real(t_))/2-min(real(t_))/2+real(tnom))/360)*360;%mattias 961122
    %add2tpl(tplf,w_tpl(ii),t_,[],option);
    
    T = t_;
    nanpars = NaN(length(obj.pars),length(T)); % until RFF outputs its parameter grid 
    tpl(ii) = qtpl(w,T.',nanpars);
    
end





end

