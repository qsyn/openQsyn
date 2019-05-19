function tpl = recgrid(obj,w,options)
%ADGRID computes template by the recursive grid method

fprintf('Calculating templates by recurcive grid\n')

dist=options.Tacc;  % to be given as option in future

fprintf('Accuracy: %g [deg], %g [dB] \n',dist(1),dist(2));
%
idx = ~strcmp({obj.pars.name},'uncint_par');
Pars = obj.pars(idx);

if length(Pars)==1
    disp('RECGRID: Only one uncertain parameter. Using Recursive Edge Grid instead');
    tpl = recedge(obj,w,options);
    return
end
    
nomPars = Pars.nom;
rangePars = [[Pars.lower ].' [Pars.upper].'];
% npar=length(Pars);
% qdist = ([Pars.upper]' - [Pars.lower]')./double([Pars.cases]');
% indgrid = obj.idxgrid(npar);
% qg = grid(Pars,0,2);

f = obj.qplant2func();
%c = qplant.pack(qg);
N = length(w);
tpl = qtpl(N); % pre-allocating
for k=1:N
    fprintf('--> for w=%g [rad/s] \n',w(k));
    [t,par]= obj.adgrid(f,w(k)*1j,rangePars,dist,[]); %the last 1 turns plotting on
    
    % unwrap according to nominal point
    tnom = cases(obj,nomPars,w(k));
    t_ = unwrap(real(t)*pi/180,180)*180/pi+1i*imag(t);
    T=t_+round((-max(real(t_))/2-min(real(t_))/2+real(tnom))/360)*360;%mattias 961122
    
    %add2tpl(tplf,w_tpl(k),t_,[],option,par);
    tpl(k) = qtpl(w(k),T.',par);
    
    if options.plot_on
        if k==1
            h=figure('Name','Template Computation');
            col = lines(N);
        end
        tpl(k).show(h,'color',col(k,:));
        drawnow
    end
end


end


