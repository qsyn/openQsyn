function tpl = crff(w_tpl,options)
%CRFF computes rff template
%   Detailed explanation goes here

disp(['Calculating templates using the Real Factored Form method']);

% if ~isempty(Par_range)
%     index2=[1 find(Par_name==',') ,  length(Par_name) ];
%     for eval_p=1:length(index2)-1
%         P_name=Par_name(index2(eval_p)+1:index2(eval_p+1)-1);
%         P_range=Par_range(eval_p,[1 2]);
%         eval([P_name,'=P_range;'])
%     end
% end


ndist=[0.25 1]; % to replace with options.accuracy?
dist=ndist(1);

if (abs(360-fix(360/dist)*dist)>eps)
    error('360 must be an integer multiple of the angular resolution, dist(1) [deg]');
end

%[sn1,sn2]=size(Un_t1);
%[sd1,sd2]=size(Un_t2);


tpl = qtpl(length(w_tpl));
for ii=1:length(w_tpl)
    disp(['--> for w=',num2str(w_tpl(ii)),' [rad/sec]']);
    %s = 1j*w_tpl(ii);
    t_=[];
    
    % compute numerator
    for k = 1:length(obj.num)
        rffElement = obj.num(k);
        if strcmp(rffElement.type,{'gain','delay','uns','int'})
            t = rffElement.rffel(w,dist); 
        else
            t = rffElement.rffpz('z',dist);
        end
        t_ = rffmul(t,t_,dist);
    end
    
    if ~isempty(t_)
        num_r=t_;
    else
        num_r=0;
    end
    
    t_ = [];
    for k = 1:length(obj.den)
        rffElement = obj.den(k);
        if strcmp(rffElement.type,{'gain','delay','uns','int'})
            t = rffElement.rffel(w,dist); 
        else
            t = rffElement.rffpz('p',dist);
        end
        t_ = rffmul(t,t_,dist);
    end
    
    [n,nn]=size(t_);
    if n > 0
        den_r = [t_((n/2+1:n),:) ; t_((1:n/2),:)];
    else
        den_r = [0];
    end
    t_=rffmul(num_r,-den_r,dist);
    
    %Jens: .'
    %if length(Kn1) > 0 ,    t_=t_+c2n(eval(Kn1),'unwrap').';      end;
    %if length(Kn2) > 0,     t_=(t_)-c2n(eval(Kn2),'unwrap').';   end;
    
    %%% certain RFF element must be implamented inside their functions!!!
    
    [n,nn]=size(t_);
    if ((n*nn) > 1)
        min_t=t_(1:n/2);
        max_t=t_(n:-1:(n/2+1));
        [bb_,index_min]=sort(real(min_t));
        [bb_,index_max]=sort(real(max_t));
        min_t = min_t(index_min);
        max_t = max_t(index_max);
        max_t = max_t(length(max_t):-1:1);
        t_=[min_t ; max_t];
        %t_=cltmp(t_,ndist); % to add as a qtpl method
    end
    %if length(n_dif)>1  % add uncertain differentiators
    %    t_=tplop('A+B',t_,c2n((1j*w_tpl(ii)).^n_dif(1:(length(n_dif)-1)),'unwrap'));
    %end
    %tnom=real(tab_le1([w_nom(:),nom(:)],w_tpl(ii))); %Mattias 960903
    %t_=qunwrap(t_(:));%mattias 961122
    %t_=t_+round((-max(real(t_))/2-min(real(t_))/2+real(tnom))/360)*360;%mattias 961122
    %add2tpl(tplf,w_tpl(ii),t_,[],option);
    
    T = t_;
    tpl(ii) = qtpl(w_tpl(k),T.',[]);
    
end





end

