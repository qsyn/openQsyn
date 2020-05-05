function tpl = recedge(obj,w,options)
%ADEGDE copmputes templates via the recurcive edge method

%fixed: the parameters retuerned do no match the actual response
%fixed: set prune method by OPTIONS

fprintf('Calculating templates by recurcive edge grid\n')

prune_on = 1;
Tacc=options.Tacc;

idx = ~strcmp({obj.pars.name},'uncint_par');
Pars = obj.pars(idx);
npar=length(Pars);
qdist = ([Pars.upper]' - [Pars.lower]')./double([Pars.cases]');
indgrid = obj.idxgrid(npar);
if length(Pars)>1
    qg = grid(Pars,2,0);    % reordered inputs!
else
    qg = linspace(Pars.lower,Pars.upper,Pars.cases); % single uncartain parameter
end

f = obj.qplant2func();

%C = num2cell(qg,2);
%C{end+1} = 0;
c = qplant.pack(qg);

tpl = qtpl(length(w)); % pre-allocating
for kw = 1:length(w)
    fprintf('--> for w=%g [rad/s] \n',w(kw));
    s = 1j*w(kw);
    Tg = qplant.funcval(f,c,s);
    T = [];
    Qpar = []; %Qpar = qg;
    for k=1:npar
        ind0=find(~indgrid(k,:));
        ind1=find(indgrid(k,:));
        for l=1:length(ind0)
            Td=Tg(ind1(l))-Tg(ind0(l));
            Td=rem(real(Td)+540,360)-180+1j*imag(Td);
            qfix=max(qg(:,ind1(l))-qg(:,ind0(l)))/qdist(k);
            if abs(real(Td)/Tacc(1)+1j*imag(Td)/Tacc(2))>=1 || qfix>1
                [T1,Q1] = recedge1(f,s,qg(:,ind0(l)),qg(:,ind1(l)),Tg(ind0(l)),Tg(ind1(l)),Tacc,qfix);
                T=[T T1];
                Qpar=[Qpar Q1];
            end
        end
    end
    T = unwrap(real(T)*pi/180)*180/pi + 1i*imag(T); % unwrap again
    
    if prune_on && numel(T)>3
        %idx=boundary(real(T)',imag(T)',0.4); % replaces PRUNE (introduced in R2014b)
        [~,idx] = prune(T,options.Tacc);
    else
        idx=1:length(T);
    end
    if options.plot_on
        col = lines(length(w));
        scatter(real(T),imag(T),5,col(kw,:)); hold on
        scatter(real(T(idx)),imag(T(idx)),10,col(kw,:),'marker','o');
    end
    tpl(kw) = qtpl(w(kw),T(idx).',Qpar(:,idx));
end

end


function [Tnew,Qpar] = recedge1(trf,s,qmin,qmax,Tmin,Tmax,Tacc,qdist)
%RECEDGE subroutine used by ADEDGE
%   [Tnew,Qpar]=recedge(trf,s,qmin,qmax,Tmin,Tmax,Tacc,qdist);
%
%   Adapted from Qsyn. Original Author: M Nordin

%global prune_on
qbreak=abs(qmin-qmax);
qbreak=min(qbreak((qbreak~=0)));
if qbreak<1e-8
    %prune_on=0 % return value???
    return
end
qnew=0.5*(qmin+qmax);
qdist=0.5*qdist;

%avoid NaN from ~plant.m = xxxplant.m  peo 970318
c=qplant.pack(qnew);
Tcurr = qplant.funcval(trf,c,s);
if isnan(Tcurr)
    %c=obj.pack(qnew+eps*ones(size(qnew)));
    qnew = qnew+eps*ones(size(qnew));
    c=obj.pack(qnew);
    Tcurr = qplant.funcval(trf,c,s);
end

Qpar=qnew;
Tnew=Tcurr;

Td=Tcurr-Tmin;
Td=rem(real(Td)+540,360)-180+1j*imag(Td);
%The wrapped distance
if abs(real(Td)/Tacc(1)+1j*imag(Td)/Tacc(2))>1 || qdist>=1
    [T1,Q1]=recedge1(trf,s,qmin,qnew,Tmin,Tcurr,Tacc,qdist);
    Tnew=[Tnew T1];
    Qpar=[Qpar Q1];
end
Td=Tmax-Tcurr;
Td=rem(real(Td)+540,360)-180+1j*imag(Td);
if abs(real(Td)/Tacc(1)+1j*imag(Td)/Tacc(2))>1 || qdist>=1
    [T1,Q1]=recedge1(trf,s,qnew,qmax,Tcurr,Tmax,Tacc,qdist);
    Tnew=[T1 Tnew];
    Qpar=[Q1 Qpar];
end
end
