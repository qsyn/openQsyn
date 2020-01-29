function [T,Qpar]=adgrid(trf,s,qpar,Tacc,n)
%ADGRID   template computation with the recursive grid method
%
%   function[T,Qpar]=adgrid(trfun,s,qpar,Tacc,n,plot_on);
%
%   Subroutine to RECGRID
%
% Template generation for transfer functions with parametric uncertainty.
% The method recursively makes a grid over the parameter space finer
% and finer until the prescribed accuracy is achieved. The final result (in Nichols
% form) is a set of points representing the border of the value set, after pruning
% (see the m-function T=prune(t,Tacc).
%
% Outputs:
%
% T      -  final template vector, with elements in Nichols form (dB + j*dB)
%
% Qpar   -  parameter matrix for the final template points. Each column represents
%           one plant case. The parameters are in the same order as in qpar.
%
% Inputs:
%
% trf    -  user defined function [t]=trfun(q,s) that gives the frequency
%           response (in Nichols form) for all parameter combinations
%           in a column of q. s is j*omega
% s      -  j*omega
% qpar  -   [qmin,qmax] The range of the uncertain parameters. The first
%           column is the minimum values, and second column is the maximum values.
%           if qmin(i)=qmax(i) the algorithm will detect that this is a constant
%           parameter.
% Tacc   -  [degree_accuracy , dB_accuracy], prescribed 2-norm accuracy in Nichols form
% n      -  (Optional) Initial grid of the parameterspace. If n(i)=1 then parameter
%           i will be considered to be constant and the mean value of qmin(i)
%           and qmax(i) will be used. Default: n=3*(number of uncertain parameters)
%
% Reference:
% [B. Cohen, M. Nordin and P.-O. Gutman, Recursive Grid Methods to Compute Value
% Sets for Uncertain Transfer Functions, ACC 1995]

% Adopted from original Qsyn toolbox
% Author: M Nordin
% Version Upgrade: A. & Y. Greenhut


qmin=qpar(:,1);
qmax=qpar(:,2);

if ~(exist('Tacc')==1); Tacc=[]; end
if ~(exist('n')==1); n=[]; end %J
if ~(exist('plot_on')==1); plot_on=[]; end

if isempty(Tacc), Tacc=[5 5]; end  % default accuracy
if isempty(n)
    n=3*ones(size(qmin));
else
    n=n(:);
end
nconst=(qmin==qmax) | n==1;  % The constant parameters!  --> deprecated

% if sum(1-nconst)==0
%     disp('All parameters constant');
%     Qpar=0.5*qmin+0.5*qmax;
%     T=c2n(feval(trf,Qpar,s));
% elseif sum(1-nconst)==1     
%     disp('ADGRID:Only one uncertain parameter. Using adedge(rfun,s,qpar,Tacc,n) instead');
%     [T,Qpar]=adedge(trf,s,qpar,Tacc,n);
% else
    disp(['ADGRID: ',num2str(sum(1-nconst)),' uncertain and ',num2str(sum(nconst)),' constant parameter(s). Accuracy [',num2str(Tacc(1)),' deg ',num2str(Tacc(2)),'dB]']);
    qconst=0.5*(qmin+qmax).*nconst;
    qmin=qmin(logical(1-nconst));
    qmax=qmax(logical(1-nconst));
    n=n(logical(1-nconst));
    Nmul=zeros(length(nconst),length(n));
    ind=find(1-nconst);
    for k=1:length(n)
        Nmul(ind(k),k)=1;
    end
    
    % Create Global Help Variables --> taken off! (Daniel R.)
    %global Nfun prune_on
    phandle=[];
    N=length(n);
    %Nfun=0;
    clear regrid1 
    %prune_on=1;
    j1=0:2^length(n)-1;
    indgrid=zeros(size(n))*zeros(size(j1));
    for k=1:N
        indgrid(k,:)=rem(fix(j1/2^(k-1)),2);
    end
    % Skapa grid och plotta
    qg=Nmul*qgrid(2*ones(size(qmin)),qmin,qmax)+qconst*ones(1,2^length(n));
    %T=c2n(feval(trfun,qg,s));
    c = qplant.pack(qg);
    T = qplant.funcval(trf,c,s);
    Qpar=qg;
    Told=T;
    [T1,Q1,Nfun,prune_on]=recgrid1(trf,s,n,qmin,qmax,T,Tacc,phandle,indgrid,Nmul,qconst);
    Qpar=[Qpar Q1];
    if prune_on
        [T, tindex]=prune([T T1],Tacc);
        Qpar=Qpar(:,tindex);
        if plot_on
            set(phandle,'XData',real([T,T(1)]),'YData',imag([T,T(1)]),'color','b','linestyle','-');
        end
    else
        T=wrap([T T1]);
        disp('    ADGRID: A singularity occured (imaginary pole or zero?), no pruning performed')
    end
    disp(['# function evaluations = ',num2str(Nfun)]);
    disp(['Final Border Size = ',num2str(length(T))]);
% end
%if plot_on;  close; end;
%pause


end

function [T,Qpar,nfun,doPrune] = recgrid1(trfun,s,n,qmin,qmax,Tgrid,Tacc,phandle,indgrid,Nmul,qconst)
%RECGRID    Subroutine used by ADGRID
%
% [T,Qpar]=recgrid(trfun,s,n,qmin,qmax,Tgrid,Tacc,phandle,indgrid,Nmul,qconst);
%
% Remark:   This code is difficult to read, due to the complicated handling of
%           n-dimensional grids, without support for them in Matlab.
%

% Author: M Nordin
% Version Upgrade:A. &  Y. Greenhut

%global Nfun prune_on
persistent Nfun % conuter for function evaluations
persistent prune_on;
if isempty(Nfun), Nfun=0; end 
if isempty(prune_on), prune_on=1; end 

N=length(n);
j=0:prod(n)-1;
qg=zeros(N,1)*zeros(size(j));
nn=max([n-1 ones(size(n))]')';
jj=0:prod(nn)-1;
ind0=zeros(N,1)*zeros(size(jj));
ind1=zeros(N,1)*zeros(size(jj));
Igrid=indgrid(N,:)*(n(N)-1);
for k=1:N
    qg(k,:)=rem(fix(j/prod(n(1:k-1))),n(k))+1;
    qg(k,:)=((n(k)-qg(k,:))*qmin(k)+(qg(k,:)-1)*qmax(k))*(1/(n(k)-1));
    ind0(k,:)=rem(fix(jj/prod(nn(1:k-1))),nn(k))+1;
    ind1(k,:)=ind0(k,:)+(n(k)>1);
end
n0=zeros(size(jj));
n1=zeros(size(jj));
n0=ind0(N,:)-1;
n1=ind1(N,:)-1;
for k=N-1:-1:1
    Igrid=n(k)*Igrid+indgrid(k,:)*(n(k)-1);
    n0=n(k)*n0+ind0(k,:)-1;
    n1=n(k)*n1+ind1(k,:)-1;
end
Igrid=Igrid+1;
Irest=ones(1,length(qg));
Irest(Igrid)=zeros(1,length(Igrid));
T=zeros(1,N);
T(Igrid)=Tgrid;
Qgrid=Nmul*qg(:,logical(Irest))+qconst*ones(1,sum(logical(Irest)));
Qpar=Nmul*qg(:,logical(Irest))+qconst*ones(1,sum(logical(Irest)));  %New parameter combinations
%Tnew=c2n(feval(trfun,Qpar,s));
c = qplant.pack(Qpar);
Tnew = qplant.funcval(trfun,c,s);
T(logical(Irest))=Tnew;   % New values needed
Nfun=Nfun+length(Irest)-length(Igrid);
n0=n0+1;
n1=n1+1;
ncorn0=ones(size(n))*n0;
ncorn1=ones(size(n))*n1;
ncorn0(1,:)=ncorn0(1,:)+1;
ncorn1(1,:)=ncorn1(1,:)-1;
for j1=2:N
    ncorn0(j1,:)=ncorn0(j1,:)+prod(n(1:j1-1));
    ncorn1(j1,:)=ncorn1(j1,:)-prod(n(1:j1-1));
end
T0=T(n0);
T1=T(n1);
Tcorn0=zeros(size(ncorn0));
Tcorn1=zeros(size(ncorn1));
for k=1:N
    % calculate the wrapped distance (factors n*360 are removed!)
    Tcorn0(k,:)=rem(real(T(ncorn0(k,:)))-real(T0)+540,360)-180+sqrt(-1)*imag(T(ncorn0(k,:)))-sqrt(-1)*imag(T0);
    Tcorn1(k,:)=rem(real(T(ncorn1(k,:)))-real(T1)+540,360)-180+sqrt(-1)*imag(T(ncorn1(k,:)))-sqrt(-1)*imag(T1);
end
nrec0=fix(abs(real(Tcorn0)*(1.1/Tacc(1))+sqrt(-1)*imag(Tcorn0)*(1.1/Tacc(2))));
nrec1=fix(abs(real(Tcorn1)*(1.1/Tacc(1))+sqrt(-1)*imag(Tcorn1)*(1.1/Tacc(2))));

if ~isempty(phandle)
    set(phandle,'XData',real(Tnew),'YData',imag(Tnew),'marker','o');
end
for j=jj(max(nrec0)>0 | max(nrec1)>0)+1
    qminrec=qg(:,n0(j));
    qmaxrec=qg(:,n1(j));
    Igrid=zeros(1,length(indgrid));
    for k=1:N
        Igrid=Igrid+indgrid(k,:)*prod(n(1:k-1));
    end
    Igrid=Igrid+n0(j);
    Trec=T(Igrid);
    if min( abs(qmaxrec-qminrec)./max(abs([qmaxrec,qminrec])')')<1e-8
        %disp('singularity in template, breaking the recursion');
        prune_on=0;
    else
        [T1,Q1,Nfun]=recgrid1(trfun,s,...
            min([max([nrec0(:,j)'+2;nrec1(:,j)'+2]);5*ones(1,N)])',...
            qminrec,qmaxrec,Trec,Tacc,phandle,indgrid,Nmul,qconst);
        Tnew = [Tnew, T1];
        Qpar = [Qpar, Q1];
    end
end
T=Tnew;
if length(Tnew)>1000 && prune_on
    Qpar=[Qgrid,Qpar];
    %	insert('tmp.mat',[Tgrid Tnew],'T','r')
    [T, tindex]=prune([Tgrid Tnew],Tacc);
    Qpar=Qpar(:,tindex);
    if ~isempty(phandle)
        set(phandle,'XData',real([T T(1)]),'YData',imag([T T(1)]),'linestyle','-');
    end
    %disp(['Still gridding. Function evaluations:',num2str(Nfun)]); % !!! DR
end
nfun = Nfun;
doPrune = prune_on;

end